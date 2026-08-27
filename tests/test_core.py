from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

from paper2_core import extract_continuous_drawdown_events_nd, fit_step_distributions, xy_from_latlon


class CoreTests(unittest.TestCase):
    def test_local_projection_origin_is_centered(self) -> None:
        x, y = xy_from_latlon(np.array([-46.0, -46.01, -45.99]), np.array([51.0, 51.02, 50.98]))
        self.assertAlmostEqual(float(x.mean()), 0.0, places=8)
        self.assertAlmostEqual(float(y.mean()), 0.0, places=8)

    def test_continuous_drawdown_finds_last_radial_maximum(self) -> None:
        points = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [1.4, 0.0], [1.0, 0.0]])
        events = extract_continuous_drawdown_events_nd(points, delta=0.8)
        self.assertGreaterEqual(len(events), 1)
        self.assertEqual(events[0].start_idx, 0)
        self.assertEqual(events[0].endpoint_idx, 2)
        self.assertEqual(events[0].trigger_idx, 4)

    def test_distribution_fits_are_finite(self) -> None:
        result = fit_step_distributions(np.array([1.0, 1.5, 2.0, 3.0, 5.0, 8.0]))
        self.assertEqual(result["n_steps"], 6)
        self.assertTrue(np.isfinite(result["aic_exp"]))
        self.assertTrue(np.isfinite(result["aic_lomax"]))

    def test_frozen_family_a_identity(self) -> None:
        frame = pd.read_csv(ROOT / "results/last_record_decomposition/family_a_event_metrics.csv.gz")
        error = np.abs(
            frame[["E_high", "E_low", "E_union"]].to_numpy(float)
            - frame[["R_high", "R_low", "R_union"]].to_numpy(float)
            - frame[["L_high", "L_low", "L_union"]].to_numpy(float)
        ).max()
        self.assertLessEqual(float(error), 1e-12)


if __name__ == "__main__":
    unittest.main()
