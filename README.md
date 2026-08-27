# Environmental history and movement boundaries

Reproducibility repository for **“Environmental history shapes movement boundaries and search scaling.”** It contains the frozen formal analyses, independent audits, source data, derived result tables, and code for every main and supplementary figure.

There are two supported reproduction paths:

| Path | Downloads | Purpose | Typical time |
|---|---:|---|---:|
| **Frozen-result reproduction** | none | Validate the release, rerun audits and rebuild all figures from the article's frozen tables | minutes |
| **Public-source reconstruction** | about 3 GB animal data plus track-specific CHL subsets | Retrieve the original public observations and satellite inputs used to construct the analysis tables | hours, network-dependent |

The two paths are deliberately separate. A reader can verify every plotted value without downloading tens of gigabytes of sensor data, while a full reconstruction remains possible from the original public sources.

## Repository contents

- `analysis/`: formal statistical analyses used by the article.
- `audits/`: independent implementations of the last-record and behavioural checks.
- `figures/`: builders for all four main and five supplementary figures, including their source-data exporters.
- `data/source_data/`: CSV source data underlying every main and supplementary figure.
- `data/derived/`: the Source Data workbook and four Supplementary Data packages.
- `results/`: frozen formal result tables and machine-readable audit records.
- `config/chl_requests/`: the exact track-specific Copernicus Marine requests; these avoid global downloads.
- `metadata/`: analysis specifications, analysis-to-input mapping and release file manifest.
- `scripts/`: environment setup, public-data download, CHL download, input checking and release validation.

Original animal observations and complete satellite products are not redistributed in this Git repository. They remain governed by their source repositories and licenses.

## One-command installation

Requirements are Python 3.10 or newer, `git`, and enough free space for the reproduction path selected below. On Linux or macOS:

```bash
git clone https://github.com/HaoboL/environmental-history-movement-boundaries.git
cd environmental-history-movement-boundaries
bash scripts/setup_environment.sh
```

The setup script creates `.venv`, installs all analysis and download dependencies, installs [`radial-drawdown`](https://github.com/HaoboL/radial-drawdown), validates the release manifest and runs the unit tests. It does **not** download raw data.

The same operation is available as:

```bash
make setup
```

To activate the environment in a later shell:

```bash
source .venv/bin/activate
```

## Path 1: reproduce the submitted results without raw downloads

Run the release validation and unit tests:

```bash
make verify
```

This checks:

1. the SHA-256 and byte count of every published repository file;
2. the exact event-level identities \(E=R+L\) and \(C=\rho-\tau\);
3. one-to-one event keys and the public Source Data workbook;
4. the absence of private paths and credential-like strings;
5. the independent analysis tests.

Rebuild every main and supplementary figure:

```bash
make figures
```

The builders write PDF, SVG, PNG, colour-vision previews and build audits to `output/figures/`. They read only frozen formal results, so they do not rerun movement segmentation, CHL extraction, bootstrapping, randomisation or model fitting.

Equivalent explicit commands are:

```bash
.venv/bin/python scripts/validate_release.py
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python figures/build_main_figures.py
.venv/bin/python figures/build_supplementary_figures.py
```

## Path 2: download the original public observations

The default public-data command downloads the four animal-data sources used in the article:

```bash
make download-public
```

This is equivalent to:

```bash
.venv/bin/python scripts/download_public_data.py \
  --dataset all \
  --uesaka-profile gps
```

Downloads are written beneath `external_inputs/raw_sources/`, use `.part` files for interruption-safe resumption, and are checked against repository byte counts and Dryad SHA-256 digests when those are available. Existing valid files are skipped.

### Animal-data download profiles

| Source | Default content | Approximate transfer | Authentication |
|---|---|---:|---|
| Goto *et al.* wandering albatross | authors' complete public Git repository | small | none |
| Uesaka *et al.* Crozet albatross | all 55 GPS `_G.csv` records | 2.40 GB | none |
| USGS Hawaiian seabirds | deployments, Laysan e-obs, and red-footed-booby GPS/FastLog/WetDry files | 371 MB | none |
| Short-tailed shearwater | complete Dryad record | about 51 MB | none |

Preview the resolved file inventory without transferring data:

```bash
.venv/bin/python scripts/download_public_data.py --dataset all --dry-run
```

Download one source only:

```bash
.venv/bin/python scripts/download_public_data.py --dataset goto
.venv/bin/python scripts/download_public_data.py --dataset uesaka --uesaka-profile gps
.venv/bin/python scripts/download_public_data.py --dataset usgs
.venv/bin/python scripts/download_public_data.py --dataset shearwater
```

The complete Uesaka accelerometer, GPS, magnetometer and pressure archive is about 63.5 GB. It is not the default because the article's landing reconstruction uses the GPS subset. A deliberate full download requires an explicit size acknowledgement:

```bash
.venv/bin/python scripts/download_public_data.py \
  --dataset uesaka \
  --uesaka-profile full \
  --accept-large-download
```

The downloader never reads or stores API tokens: all animal repositories used here allow public anonymous download.

## Download the track-specific chlorophyll-*a* subsets

CHL comes from Copernicus Marine product `OCEANCOLOUR_GLO_BGC_L4_MY_009_104`, dataset `cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D`, frozen product version `202603`. The repository contains the exact spatial and temporal request rows for Goto albatrosses, Laysan albatrosses, red-footed boobies and short-tailed shearwaters.

Copernicus Marine requires a free account. Configure credentials locally with the official CLI; credentials are stored by the Copernicus client outside this repository:

```bash
source .venv/bin/activate
copernicusmarine login
```

Inspect the first two requests without downloading:

```bash
.venv/bin/python scripts/download_chl.py --all --dry-run --max-downloads 2
```

Run all frozen requests:

```bash
make download-chl
```

Each manifest is resumable at the file level: rerunning the command skips non-empty outputs. Progress lines report completed/total requests, percentage, elapsed time and ETA. Status CSVs and a JSON summary are written to `external_inputs/download_status/chl/`.

To run or resume a single collection:

```bash
.venv/bin/python scripts/download_chl.py \
  --manifest config/chl_requests/shearwater_chl_requests.csv
```

To resume at a particular manifest row or run a small batch:

```bash
.venv/bin/python scripts/download_chl.py \
  --manifest config/chl_requests/goto_chl_requests.csv \
  --start-row 200 \
  --max-downloads 25
```

No global CHL field is downloaded. `download_chl.py` ignores the historical machine-specific output paths retained in the frozen manifests and writes portable outputs to `external_inputs/environment/<collection>/`.

## Check downloaded inputs

After the animal and CHL downloads finish:

```bash
make check-inputs
```

The checker reports the number of expected public source files and the exact number of completed CHL requests for each collection. A missing source produces a non-zero exit status.

The resulting raw-source layout is:

```text
external_inputs/
├── raw_sources/
│   ├── goto_wandering_albatross/
│   ├── uesaka_crozet_multisensor/{metadata,raw}/
│   ├── usgs_hawaiian_seabirds/{metadata,raw}/
│   └── short_tailed_shearwater/{metadata,raw}/
├── environment/
│   ├── goto_chl/
│   ├── laysan_chl/
│   ├── laysan_dateline_chl/
│   ├── red_footed_booby_chl/
│   └── shearwater_chl/
└── download_status/chl/
```

## From public records to formal analyses

The computational chain is:

```text
public tracks and behavioural sensors
        +
track-specific daily CHL subsets
        ↓
quality control, temporal alignment and radial-drawdown segmentation
        ↓
frozen event-level analysis inputs
        ↓
formal analyses and independent audits
        ↓
source-data tables and article figures
```

`metadata/ANALYSIS_MANIFEST.csv` maps each manuscript component to its program, frozen result directory and external input requirement. Analysis programs use `PAPER2_PROJECT_ROOT` or paths relative to the repository; no private absolute path is required.

Several formal programs start from frozen event-level tables rather than repeating source-specific telemetry cleaning. This separation preserves the exact event population subjected to confirmatory statistics. The frozen results and all plotted source values are included in the repository, while the public-source downloader reconstructs the original observation layer. Before submission, the versioned release will permanently archive the larger analysis-input layer alongside the tagged code release so that both the raw-to-event and event-to-figure paths are fixed to one DOI.

For the current living manuscript repository, the immediately executable guarantees are therefore:

- all release checks, independent audits and figures run without external data;
- all original public animal and CHL records can be downloaded with the commands above;
- formal analysis programs, specifications, seeds and required input filenames are public;
- large analysis-level inputs will be frozen in the submission archive rather than committed as ordinary Git blobs.

## Public data provenance

- Goto *et al.* wandering-albatross tracks: [PNAS article](https://doi.org/10.1073/pnas.2312851121) and [authors' public repository](https://github.com/YusukeGoto510/Data_Wandering_Albatross)
- Uesaka *et al.* Crozet high-frequency multisensor data: [Dryad](https://doi.org/10.5061/dryad.tx95x6b2j)
- Laysan albatross and red-footed booby biologging data: [USGS ScienceBase](https://doi.org/10.5066/P9NTEXM6)
- Short-tailed shearwater GPS and accelerometry: [Dryad](https://doi.org/10.5061/dryad.j9k60)
- Daily chlorophyll-*a*: [Copernicus Marine product DOI](https://doi.org/10.48670/moi-00281)

The download implementations follow the [Dryad REST API documentation](https://github.com/datadryad/dryad-app/blob/main/documentation/apis/README.md), the USGS ScienceBase item API, and the official [Copernicus Marine `subset` interface](https://help.marine.copernicus.eu/en/articles/7972861-copernicus-marine-toolbox-cli-subset). Copernicus credential setup is documented by [Copernicus Marine](https://help.marine.copernicus.eu/en/articles/8185007-copernicus-marine-toolbox-credentials-configuration).

See `data/README.md` for redistribution boundaries and the roles of the included derived files.

## License and citation

Code is released under the MIT License. Data files in this repository are article-specific derived source data; original observations remain governed by their source repositories. Cite the article when available and use `CITATION.cff` for this living release in the interim. The submission-matched version will receive an immutable tag and permanent archive DOI after figures and text are frozen.
