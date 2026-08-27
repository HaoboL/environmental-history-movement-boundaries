from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative: str):
    specification = importlib.util.spec_from_file_location(name, ROOT / relative)
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


class DownloadWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.chl = load_module("download_chl_for_test", "scripts/download_chl.py")
        cls.public = load_module("download_public_for_test", "scripts/download_public_data.py")

    def test_all_chl_manifests_are_executable(self) -> None:
        manifests = sorted((ROOT / "config/chl_requests").glob("*.csv"))
        self.assertEqual(len(manifests), 5)
        total = 0
        for manifest in manifests:
            rows = self.chl.read_manifest(manifest)
            self.assertGreater(len(rows), 0)
            total += len(rows)
            outputs: set[str] = set()
            for row in rows:
                command = self.chl.command_for(row, ROOT / "external_inputs/environment/test")
                self.assertEqual(command[:2], ["copernicusmarine", "subset"])
                self.assertNotIn(str(row.get("output_dir", "")), command)
                self.assertNotIn(row["output_filename"], outputs)
                outputs.add(row["output_filename"])
        self.assertEqual(total, 604)

    def test_public_source_inventory_is_frozen(self) -> None:
        self.assertEqual(self.public.UESAKA_DOI, "10.5061/dryad.tx95x6b2j")
        self.assertEqual(self.public.SHEARWATER_DOI, "10.5061/dryad.j9k60")
        self.assertEqual(len(self.public.USGS_REQUIRED), 5)
        self.assertTrue(self.public.GOTO_REPOSITORY.endswith("Data_Wandering_Albatross.git"))


if __name__ == "__main__":
    unittest.main()
