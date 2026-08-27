# Environmental history and movement boundaries

Reproducibility release for the study **“Environmental history shapes movement boundaries and search scaling.”** The repository contains the frozen analysis programs, independent audits, figure builders, source data, derived result tables, and machine-readable manifests supporting the manuscript.

## What is included

- `analysis/`: the seven formal Paper 2 analyses and the short-tailed-shearwater behaviour bridge. Public filenames describe the scientific task rather than internal run numbers, and portable copies use `PAPER2_PROJECT_ROOT`.
- `audits/`: independent checks of the last-record decomposition and behavioural analyses.
- `figures/`: builders for the four main and five supplementary figures, plus the Supplementary Data packager.
- `data/source_data/`: CSV source data underlying every main and supplementary figure.
- `data/derived/`: the submitted Source Data workbook and four Supplementary Data packages.
- `results/`: frozen formal result tables and audit records read by the public figure and verification workflows.
- `metadata/`: frozen analysis specifications and the release manifest.

The repository does **not** redistribute original tracking or satellite products. Those data remain available from their original public repositories and should be cited under their own terms.

## Quick verification

Use Python 3.10 or newer. Create an isolated environment and install the listed dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install "git+https://github.com/HaoboL/radial-drawdown.git"
```

Then run the release checks:

```bash
python scripts/validate_release.py
python -m unittest discover -s tests -v
```

The first command verifies every published file against `metadata/FILE_MANIFEST.csv`, checks the exact event-level identities \(E=R+L\) and \(C=\rho-\tau\), and checks the public Source Data workbook against its CSV exports.

To rebuild all figures from the frozen formal results:

```bash
python figures/build_main_figures.py
python figures/build_supplementary_figures.py
```

Outputs are written to `output/figures/`. These figure builders do not rerun movement segmentation, environmental extraction, bootstrapping, randomisation, or model fitting.

## Full analysis reruns

The formal analysis programs retain the original public-data input layout and frozen seeds. Set `PAPER2_PROJECT_ROOT` to a project tree containing the downloaded inputs and run the relevant file in `analysis/`. `metadata/ANALYSIS_MANIFEST.csv` maps manuscript components to programs, frozen results and required external inputs.

Movement segmentation uses [`radial-drawdown`](https://github.com/HaoboL/radial-drawdown). The companion methods repository is a dependency; it is not duplicated or modified here.

## Original public data

- Goto *et al.* wandering-albatross tracks: [PNAS article and accompanying repository](https://doi.org/10.1073/pnas.2312851121)
- Uesaka *et al.* Crozet high-frequency multisensor data: [Dryad](https://doi.org/10.5061/dryad.tx95x6b2j)
- Laysan albatross and red-footed booby biologging data: [USGS data release](https://doi.org/10.5066/P9NTEXM6)
- Short-tailed shearwater GPS and accelerometry: [Dryad](https://doi.org/10.5061/dryad.j9k60)
- Daily chlorophyll-*a*: Copernicus Marine product `OCEANCOLOUR_GLO_BGC_L4_MY_009_104`, dataset `cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D`, frozen version `202603` ([product DOI](https://doi.org/10.48670/moi-00281))

See `data/README.md` for redistribution boundaries and expected input roles.

## License and citation

Code is released under the MIT License. Data files in this repository are derived source data for the associated article; the original observations remain governed by their source repositories. Cite the article when available and use `CITATION.cff` for this release in the interim.
