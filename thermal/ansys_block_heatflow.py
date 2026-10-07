# Ansys Mechanical 2025 R1 scripting (IronPython) - run from Automation > Scripting.
# Replaces volumetric Internal Heat Generation on the 22 component blocks with a
# Heat Flow on each block's bottom face (the face touching its PCB). Same watts.
# Board IHG loads (ENIAC, C6, G6) are left alone.
#
# Why: blocks use Chip_Mold (k ~0.7 W/m-K, mold compound). Volumetric heat in a
# low-k block traps heat inside it and gives unphysical hot spots. Real packages
# dump most heat into the board through leads/thermal pad, so a face load at the
# block-board interface is the closer simplification.
#
# Block wattage is read from the body name, e.g. "C6_U8_TLV1117_1.659W".
# Block -> board by name prefix (C6_, G6_, ENIAC_) matched to "<prefix>_PCB".
import re

model = ExtAPI.DataModel.Project.Model
analysis = model.Analyses[0]
bodies = [b for a in ExtAPI.DataModel.GeoData.Assemblies for p in a.Parts for b in p.Bodies]

def vmean(body):
    vs = body.Vertices
    n = float(len(vs))
    return [sum(v.X for v in vs) / n, sum(v.Y for v in vs) / n, sum(v.Z for v in vs) / n]

def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

def sub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]

boards = {}
for b in bodies:
    m = re.match(r'^(C6|G6|ENIAC)_PCB', b.Name)
    if m:
        boards[m.group(1)] = b

# 1) suppress block IHG loads, keep board ones
keep = ('(ENIAC)', '(C6)', '(G6)')
nsup = 0
for ihg in analysis.GetChildren(DataModelObjectCategory.InternalHeatGeneration, True):
    if not ihg.Name.strip().endswith(keep) and not ihg.Suppressed:
        ihg.Suppressed = True
        nsup += 1

# 2) Heat Flow on bottom face of each block
existing = [h.Name for h in analysis.GetChildren(DataModelObjectCategory.HeatFlow, True)]
total = 0.0
made = []
for b in bodies:
    m = re.match(r'^(C6|G6|ENIAC)_(.+)_([0-9.]+)W', b.Name)
    if not m:
        continue
    brd, tag, watts = m.group(1), m.group(2), float(m.group(3))
    name = 'Heat Flow (%s_%s %.3fW)' % (brd, tag, watts)
    total += watts
    if name in existing:
        continue
    cb = vmean(b)
    toward = sub(vmean(boards[brd]), cb)  # direction from block toward its board
    # bottom face = face whose centroid is furthest toward the board
    face = max(b.Faces, key=lambda f: dot(sub(list(f.Centroid), cb), toward))
    sel = ExtAPI.SelectionManager.CreateSelectionInfo(SelectionTypeEnum.GeometryEntities)
    sel.Ids = [face.Id]
    hf = analysis.AddHeatFlow()
    hf.Location = sel
    hf.Magnitude.Output.DiscreteValues = [Quantity(watts, 'W')]
    hf.Name = name
    made.append(name)

print('Suppressed block IHG loads: %d' % nsup)
print('Heat Flow loads created: %d' % len(made))
print('Total block watts: %.3f W' % total)
for n in made:
    print('  ' + n)
