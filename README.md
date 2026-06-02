# Spatiotemporal Tumour Evolution: CPM and the 3F Model

Core computational framework accompanying the manuscript
"Field–Flow–Front mapping of breast cancer evolutionary dynamics in situ".
This repository implements the Cancer Progression Metric (CPM) and the
Field–Flow–Front (3F) model — the central frameworks for mapping tumour
evolutionary dynamics in situ.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── run_all.R                  # master script: runs the full CPM + 3F pipeline
├── data/
│   └── example_input/
│       └── Integrated_Data_region5.csv   # small demo dataset (one region)
├── scripts/
│   └── 01_run_CP_3F.R         # CPM construction and 3F (Field/Flow/Front) analysis
├── results/                   # output directory (figures and tables)
└── env/                       # environment / dependency specification
    ├── sessionInfo.txt        # full R sessionInfo() output
    └── requirements.txt       # Python package versions
```

## 1. System requirements

### Software dependencies and versions
- R (v4.5.1), Python (v3.13.9)
- Key R packages: Seurat (v5.4.0), Monocle2 (v2.24.1), spacexr/RCTD (v2.2.0),
  SingleR (v2.12.0), infercnv (v1.26.0), clusterProfiler (v4.18.4),
  scTenifoldNet (v1.3), survival (v3.8.3), survminer (v0.5.1), rms (v8.1-0),
  timeROC (v0.4), randomForest (v4.7-1.2), lme4 (v2.0.1), lmerTest (v3.2.1),
  mgcv (v1.9.3), hdf5r (v1.3.12), shiny (v1.11.1), plotly (v4.11.0),
  dplyr (v1.1.4), tidyr (v1.3.1), data.table (v1.17.8), ggplot2 (v4.0.1),
  sf, dbscan, FNN
- Key Python libraries: opencv-python (v4.13.0), numpy (v2.3.5),
  tifffile (v2026.1.14), nibabel (v5.3.3), networkx (v3.1),
  matplotlib (v3.10.8), plotly (v6.5.2), scipy (v1.16.3),
  scikit-image (v0.26.0)

### Operating systems (tested on)
- Windows 10/11 (development and testing environment)

### Non-standard hardware
No non-standard hardware is required. A standard desktop/workstation is
sufficient to run the demo. For full-resolution spatial datasets
(Xenium / Visium HD / MIBI), ≥64 GB RAM is recommended. No GPU is required.

## 2. Installation guide

### Instructions
1. Install R (v4.5.1) and Python (v3.13.9).
2. Install the required R packages (versions listed above; see
   `env/sessionInfo.txt`).
3. Install the required Python packages:
```bash
   pip install -r env/requirements.txt
```

### Typical install time
Approximately 20–30 minutes on a normal desktop computer, including
dependency download and installation.

## 3. Demo

A small demo dataset for a single region is provided in
`data/example_input/Integrated_Data_region5.csv`.

### Instructions to run
```bash
Rscript run_all.R
```
This runs the CPM construction and 3F (Field/Flow/Front) analysis on the demo
region.

### Expected output
Figures and tables written to `results/`, including per-cell CPM values, the
CPM field map, and the Field/Flow/Front designation for the demo region.

### Expected run time
Approximately 5 minutes on a normal desktop computer for the provided demo
region.

## 4. Instructions for use

To run on your own data, replace the input file in `data/example_input/` with
your own cell-by-feature matrix and spatial coordinates in the same format,
then run:
```bash
Rscript run_all.R
```
Individual steps can be run separately via the scripts in `scripts/`.

### Reproduction instructions
The provided scripts reproduce the core quantitative results of the manuscript
for the CPM and 3F analyses, following the parameters detailed in the Methods
section of the manuscript.

## License
This project is released under the MIT License (see `LICENSE`).
