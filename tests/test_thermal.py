"""Physical limiting cases and conservation checks for the lumped network."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np

path = Path(__file__).resolve().parents[1] / 'thermal' / 'Avbay_Steady.py'
spec = importlib.util.spec_from_file_location('avbay', path)
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)


class ThermalChecks(unittest.TestCase):
    def assert_balanced(self, result):
        self.assertLess(result['max_residual_W'], 1e-6)

    def test_loads_accounted_once(self):
        np.testing.assert_allclose(m.board_powers(), [4.727, .691, 6.643], atol=1e-12)
        self.assertAlmostEqual(m.board_powers().sum(), 12.061)

    def test_scenarios_conserve_energy(self):
        for title, *env in m.CASES:
            with self.subTest(title=title): self.assert_balanced(m.solve(*env))

    def test_zero_electronics_equalizes_internal_nodes(self):
        r = m.solve(42, 60, 900, 100, 10, P=[0, 0, 0])
        self.assert_balanced(r)
        np.testing.assert_allclose(r['T_p'], r['T_w'], atol=1e-7, rtol=0)
        self.assertAlmostEqual(r['T_a'], r['T_w'], places=7)

    def test_contact_resistance_raises_wall_not_single_region_skin(self):
        a = m.solve(42, 60, 900, 100, 10)
        b = m.solve(42, 60, 900, 100, 10, geometry=replace(m.Geometry(), contact_h=10))
        self.assert_balanced(b)
        self.assertAlmostEqual(a['T_s'], b['T_s'], places=7)
        self.assertGreater(b['T_w'], a['T_w'])

    def test_extension_and_finite_contact_conserve_energy(self):
        r = m.solve(42, 60, 900, 100, 10,
                    geometry=replace(m.Geometry(), modeled_length_m=.508, contact_h=100))
        self.assertEqual(len(r['regions']), 2)
        self.assert_balanced(r)

    def test_unequal_h_and_mounts_conserve_energy(self):
        cfg = replace(m.Settings(), h_board=2, h_wall=4, mount_conductance=(1, 2, 3))
        self.assert_balanced(m.solve(42, 60, 900, 100, 10, settings=cfg))

    def test_reject_invalid_inputs(self):
        for p in ([1, 2], [-1, 0, 0], [float('nan'), 0, 0]):
            with self.subTest(P=p), self.assertRaises(ValueError):
                m.solve(42, 60, 900, 100, 10, P=p)
        with self.assertRaises(ValueError):
            m.solve(42, 60, 900, 100, 10, geometry=replace(m.Geometry(), modeled_length_m=.6))


if __name__ == '__main__': unittest.main()
