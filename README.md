# IREC Cooling 2026-2027

Passive avionics-bay thermal model for UGA AIAA's IREC rocket, Church.

The solver represents **one isothermal node per PCB**, one sealed bay-air node,
and one inner-wall node. Individual component heat estimates are retained in
`COMPONENT_HEAT_W`; their sums plus the small-parts allowance become the full
board loads. Component blocks are not heat sources or thermal nodes in this
baseline. Results are average board temperatures, not chip junction temperatures.

## Run

Python 3.10 or newer:

```bash
python -m pip install -r requirements.txt
python thermal/Avbay_Steady.py
python -m unittest discover -s tests -v
```

Optional sensitivity cases:

```bash
python thermal/Avbay_Steady.py --h-inside 2.1
python thermal/Avbay_Steady.py --contact-h 100
python thermal/Avbay_Steady.py --length-mm 508
```

`--contact-h` is tube/coupler interface conductance in W/(m² K). Its default is
ideal contact, not a measured fit. No suggested finite value is a calibrated input.
`--length-mm 508` treats the protruding coupler as exposed to the same external
environment as the body tube. Do not use that case if a neighboring rocket section
covers the extension without modeling that section.

## Geometry and heat paths

Diameters from the user and lengths/PCB body data from
`Avbay_Floating_PCBs_v2.step`:

| Item | Value |
|---|---:|
| Body tube OD / ID | 156.718 / 152.400 mm |
| Coupler OD / ID | 152.400 / 147.828 mm |
| Body tube / coupler length | 400.05 / 508 mm |
| **Default provisional modeled length** | **400.05 mm, common overlap only** |
| Body tube / coupler radial thickness | 2.159 / 2.286 mm |
| Default outside cylindrical area | 0.196962 m² |
| Default inside cylindrical area | 0.185789 m² |
| Body tube / coupler conductivity | 0.29 W/(m K), assumed |

The actual sealed bay length and bulkhead positions are not confirmed. The
108 mm coupler extension is excluded from the default analysis domain. Change
that domain after confirming the assembly; these dimensions are not a claim that
the real sealed bay is 400.05 mm long. Endcap heat transfer is omitted, and the
modeled gas is assumed sealed. Mounting paths default to zero because they are
absent from the supplied STEP.

The two concentric layers use cylindrical resistance in series:

```text
R_body = ln(D_body_OD / D_body_ID) / (2*pi*k_body*L)
R_coupler = ln(D_coupler_OD / D_coupler_ID) / (2*pi*k_coupler*L)
R_contact = 1 / (h_contact*pi*D_interface*L), if supplied
```

Heat enters the outer skin from sunlight and exchanges with outside air, sky,
and ground. The skin connects to the inner wall through the tube layers/contact.
PCBs exchange convection through the solved bay-air node and radiation directly
with the inner-wall node. Optional `Settings.mount_conductance` adds PCB-to-wall
conduction in W/K. No imposed board or inner-wall temperature is used.

## Power bookkeeping

All component watts are **preliminary user-provided estimates**, not a completed
datasheet or bench audit. Edit individual estimates in `COMPONENT_HEAT_W`; do
not independently maintain another set of full-board totals.

| Board | Named components | Small parts | Total heat | Full-board IHG for supplied STEP body |
|---|---:|---:|---:|---:|
| C6 | 4.533 W | 0.194 W | **4.727 W** | **216,270 W/m³** |
| G6 | 0.600 W | 0.091 W | **0.691 W** | **31,587 W/m³** |
| ENIAC | 6.510 W | 0.133 W | **6.643 W** | **707,176 W/m³** |

Total: **12.061 W**. Full-board IHG is total heat divided by the actual bare-board
body volume, and must be recomputed if geometry changes. Surface areas include
both PCB faces and edges. Component block surfaces are excluded.

Verify regulator losses from their actual rail currents; do not count regulator
output energy again as regulator heat. The buck allocation may contain inductor
and other converter losses, so it is a board-level allowance, not validated IC
self-heating. Verify RF amplifier DC input and net RF output separately. The
BMT321 value remains uncertain; the manufacturer datasheet is now available at
https://d2lxe0fofddnat.cloudfront.net/BMT321-V7.4.pdf.
The top-level schematic uploads omit referenced child sheets; PCB files provide
nets, footprints and stackups, but firmware activity and transmit duty cycle are
still unknown. A measured board input-power budget is an independent check,
allowing for exported electrical/RF power and inter-board power transfers.

All three environmental scenarios use the **same electronics operating-case
budget**. The previous 7.4 V idle-power estimate is not silently mixed into the
moderate/cool cases; a separate verified idle component budget is needed.

## Environmental scenarios and results

Inputs retained from the previous solver, not independently certified as site
extremes: direct normal irradiance 900 W/m², diffuse horizontal 100 W/m², sun
elevation 45 degrees, ground albedo 0.40, solar absorptivity 0.60, infrared
emissivity 0.90. Solar loading assumes a vertical cylinder, circumferentially
averaged direct beam, isotropic diffuse sky and ground reflection. Exterior long-
wave radiation uses equal sky/ground view fractions and a clear-sky Swinbank
sky-temperature approximation. Solar absorptivity and IR emissivity are separate.

| Scenario | Air / ground | Exterior h | Solar | Bay air | C6 | G6 | ENIAC |
|---|---|---:|---|---:|---:|---:|---:|
| Hot pad | 42 / 60°C | 10 W/(m² K) | on | 67.56°C | 78.65°C | 66.97°C | 107.87°C |
| Moderate pad | 34 / 52°C | 21.6 W/(m² K) | on | 51.70°C | 64.06°C | 51.19°C | 95.84°C |
| Cool / no sun | 20.5 / 20.5°C | 22 W/(m² K) | off | 27.99°C | 42.55°C | 27.66°C | 78.52°C |

These are conditional model outputs for the default overlap-only geometry,
ideal tube contact, internal h=3 W/(m² K), no mounts and the same 12.061 W budget.
The hot case is **not a verified worst-case bound**, and the moderate case is not
an established statistical median. Calm wind and orientation require separate
sensitivity cases. Ground temperature and coating properties remain assumptions.

## Matching the simplified ANSYS Steady-State Thermal model

1. Suppress the 22 component block bodies for this board-only baseline, and
   suppress/delete their IHG and Heat Flow loads. The budget records remain in
   Python. Do not also apply component loads alongside full-board IHG.
2. Apply **full-board** IHG from the table to the three PCB bodies. The earlier
   small-parts-only residual IHG values are for a different, component-loaded
   model and are not the complete heat loads here.
3. Use the actual tube/coupler geometry and confirm thermal contacts. Scope
   external solar/convection/radiation only to genuinely exposed surfaces.
   Match the modeled axial domain or account explicitly for extra solid regions.
4. Hot-case external BCs: inward solar flux 239.910 W/m²; convection h=10,
   outside air 42°C; radiation To Ambient, emissivity 0.90, effective sky/ground
   surroundings 48.53°C. Convection and radiation use different temperatures.
5. Use internal convection on exposed board and inner-enclosure faces only.
   The Python hot-case air temperature 67.56°C is a starting bulk-temperature
   guess, not a reservoir temperature that remains valid for every geometry.
6. Iterate the shared internal bulk temperature using
   `T_air = sum(h_i*A_i*Tmean_i) / sum(h_i*A_i)` until the net internal convection
   heat flow is near zero. `Tmean_i` must be an area-weighted face temperature.
7. Use actual internal radiation view factors. Perfect enclosure requires the
   real closed bay boundaries, which are missing from the supplied STEP. Do not
   prescribe inner-wall temperature as a radiation bath.
8. Check electronics + absorbed solar = net exterior convection + radiation
   (including any additional modeled endcap/mount heat paths).

A bare-board Mechanical conduction model need not match an isothermal Python
board exactly. The four-layer KiCad stackups confirm copper heat spreading matters;
use layered/effective directional PCB properties if resolving board gradients.
Uniform volumetric board heat is a deliberate homogenization of component heat.

## Verification and remaining limits

The solver reports signed per-board, air, inner-wall, exterior-region and complete
bay energy residuals; it rejects a solution whose largest residual exceeds
1e-6 W. Automated checks cover the budgets, zero-electronics equilibrium,
separate convection coefficients, finite interface resistance, optional mounts,
exposed extension and invalid inputs. They check numerical consistency, not
experimental validation.

Remaining assumptions: steady state, isothermal PCBs and inner wall, well-mixed
sealed air, prescribed convection coefficients, all board radiation directed to
the wall, no board-to-board obstruction/view factors, no local solar gradients,
no endcap/mount heat losses by default, unverified operating heat budget. It does
not predict heating time, component junction temperatures, ventilation, fan curves
or flight aerothermal effects. Transient requires heat capacities and initial
conditions; Fluent can replace the internal lumped-air/convection assumptions.

## Original reference list

The following references were carried forward from the original project README;
their inclusion does not certify every current numerical assumption.


1. ESRA, 2026 IREC Design, Test & Evaluation Guide V1.1, Appendix G. https://www.soundingrocket.org/what-is-irec.html
2. NOAA NCEI, Engineering Weather Data, Midland TX. https://www.ncei.noaa.gov/pub/data/EngineeringWeatherData_CDROM/engwx/midland_tx.pdf
3. Bergman, Lavine, Incropera, DeWitt, Fundamentals of Heat and Mass Transfer, 7th ed., Wiley: ch. 1, 3, 7 (Churchill-Bernstein), 9 (Churchill-Chu), 13, Table A.4.
4. PVEducation, Standard Solar Spectra (ASTM G173). https://www.pveducation.org/node/430
5. Wikipedia, Albedo. https://en.wikipedia.org/wiki/Albedo
6. Swinbank, W.C. (1963), Long-wave radiation from clear skies, Q. J. R. Meteorol. Soc. 89, 339-348. https://doi.org/10.1002/qj.49708938105
7. NASA POWER climatology, launch site 31.05 N, 103.55 W. https://power.larc.nasa.gov/api/temporal/climatology/point?parameters=TS_MAX,T2M_MAX,ALLSKY_SRF_ALB&community=RE&longitude=-103.5473&latitude=31.0498&format=JSON
8. FLIR, What emissivity settings should I use when looking at PCBs and electronics? https://www.flir.com/en-ca/support/instruments2/what-emissivity-settings-should-i-use-when-looking-at-pcbs-and-electronics/
9. Wikipedia, FR-4. https://en.wikipedia.org/wiki/FR-4
10. Datasheets: TL780 https://www.ti.com/lit/ds/symlink/tl780.pdf, TPS563200 https://www.ti.com/lit/ds/symlink/tps563200.pdf, TLV1117LV https://www.ti.com/lit/pdf/SBVS160, LMX2541 https://www.ti.com/lit/ds/symlink/lmx2541.pdf, HMC361 https://www.analog.com/media/en/technical-documentation/data-sheets/hmc361s8g.pdf, BMT321 https://www.berex.com/Products/HighPowerAmplifier, SX1281 https://www.mouser.com/datasheet/2/761/DS_SX1280-1_V2.2-1511144.pdf, ADV7280A https://www.analog.com/media/en/technical-documentation/data-sheets/adv7280a.pdf, SiT5000 https://www.sitime.com/datasheet/SiT5000, NEO-M9N https://content.u-blox.com/sites/default/files/NEO-M9N-00B_DataSheet_UBX-19014285.pdf, MMC5983MA https://www.memsic.com/Public/Uploads/uploadfile/files/20220119/MMC5983MADatasheetRevA.pdf, ISL83488 https://www.renesas.com/us/en/document/dst/isl83483-isl83485-isl83488-isl83490-isl83491-datasheet
