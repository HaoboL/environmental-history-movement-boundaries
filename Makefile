PYTHON := .venv/bin/python

.PHONY: setup verify figures download-public download-chl check-inputs

setup:
	bash scripts/setup_environment.sh

verify:
	$(PYTHON) scripts/validate_release.py
	$(PYTHON) -m unittest discover -s tests -v

figures:
	$(PYTHON) figures/build_main_figures.py
	$(PYTHON) figures/build_supplementary_figures.py

download-public:
	$(PYTHON) scripts/download_public_data.py --dataset all --uesaka-profile gps

download-chl:
	$(PYTHON) scripts/download_chl.py --all

check-inputs:
	$(PYTHON) scripts/check_external_inputs.py
