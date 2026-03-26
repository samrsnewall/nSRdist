# nSRdist

MATLAB code for estimating the probability distribution of normalised sedimentation rates (NSR) from radiocarbon-dated marine sediment cores.

## Scientific Background

The accumulation rate of sediment on the ocean floor varies over time and across locations, and understanding its distribution is important for interpreting paleoclimate records. This code uses radiocarbon dates from sediment cores to estimate **normalised sedimentation rates (NSR)** — the sedimentation rate at each point in a core divided by the mean sedimentation rate for that core — and fits probability distributions to them across a large ensemble of cores.

The approach builds on the method of [Lin et al. (2014)](https://doi.org/10.1002/2013PA002539), which used Bchron age-depth models to estimate NSR. The primary development in this codebase is a probabilistic extension that propagates the full calibrated age uncertainty of each radiocarbon date, rather than using only the mode or median of the age model. Several NSR estimation methods are implemented and compared.

The fitted NSR distributions are intended for use in paleoclimate age modelling applications such as [BIGMACS](https://doi.org/10.5194/gchron-5-235-2023).

## Pipeline Overview

The code operates in two stages:

**Stage 1 — `calcData.m`**

Reads core metadata and radiocarbon data, filters dates, handles doubly-dated depths and age reversals by constructing scenarios, runs Bchron age-depth modelling, and computes NSR histories using several methods. Results are saved to a `.mat` file in `Results/`.

**Stage 2 — `fitData.m`**

Loads the saved results from Stage 1 and fits four candidate probability distributions to the NSR data: a 2-component Mixed Log-Normal (MLN), a Log-Normal (LN), a Gamma, and an Inverse Gamma. Fitting is performed under various subset choices and weighting schemes, and results are saved to a new `.mat` file.

### NSR Estimation Methods

| Method | Description |
|---|---|
| **BMedian** | Uses the median of the Bchron age-depth model at each radiocarbon dated depth. Analogous to the Lin et al. (2014) approach. |
| **BMode** | Uses the mode of the Bchron age-depth model. |
| **BSamp (BChAR / BChIR)** | Samples directly from the Bchron posterior age-depth distribution across many runs. |
| **RSR0 / RSR500 / RSR1000 / RSR1500** | Random sampling from calibrated radiocarbon age probability distributions, with an optional minimum age-difference restriction (0, 500, 1000, or 1500 years) between consecutive dates. |

---

## Requirements

### MATLAB
- MATLAB (tested on recent versions; no specific minimum version identified)
- **Statistics and Machine Learning Toolbox** — required for `fitgmdist` (mixture model fitting)
- **Parallel Computing Toolbox** — required for `parfor` loops in `calcData.m` (optional; the code will run without it if `parfor` is replaced with `for`)

### Third-Party MATLAB Function: MatCal
Radiocarbon calibration is performed using **MatCal** (Lougheed & Obrochta, 2019). The functions `matcal` and `matcalq` must be on your MATLAB path. MatCal is available at:
[https://github.com/bryanlougheed/matcal](https://github.com/bryanlougheed/matcal)

### R and Bchron
Bchron age-depth modelling is run via an R script called from MATLAB. You will need:
- **R** installed and accessible from your system PATH (so that `Rscript` can be found automatically)
- The **Bchron** R package: `install.packages("Bchron")`
- The **IntCal** R package: `install.packages("IntCal")`

---

## Setup

### 1. Clone or download the repository

```bash
git clone https://github.com/samrsnewall/nSRdist.git
```

### 2. Set the sandbox path

Open `calcData.m` and set `S.sandboxPath` to the path of the repository on your machine:

```matlab
S.sandboxPath = "/path/to/your/nSRdist_code";
```

All other paths (Bchron folders, Lin2014 data, etc.) are constructed relative to this.

> **Note:** A future update will replace this with an automatic path detection using `fileparts(mfilename('fullpath'))`, removing the need for manual setup.

### 3. Confirm Rscript is on your PATH

The code will attempt to find `Rscript` automatically. If it is not on your PATH, set `S.RscriptPath` manually in `calcData.m`:

```matlab
S.RscriptPath = "/path/to/Rscript";  % e.g. "/usr/local/bin/Rscript" on Mac/Linux
```

### 4. Add MatCal to your MATLAB path

Either place the MatCal functions in `Functions/` or add them to your MATLAB path via `addpath`.

### 5. Provide input data

The code reads core metadata from an Excel spreadsheet (`DataSheets/COPYcore40MetadataAndLin2014_2.xlsx`) and radiocarbon data from the Mulitza et al. (2022) World Atlas dataset. Set the path to the World Atlas data in `calcData.m`:

```matlab
S.WApath = "path/to/WA_Foraminiferal_Isotopes_2022";
```

---

## Running the Code

### Stage 1: Calculate NSR data

Open and run `calcData.m`. Key settings at the top of the script control which cores are analysed, which methods are run, and where outputs are saved. The main output is a `.mat` file saved to `Results/`.

### Stage 2: Fit distributions

Open `fitData.m` and set `filepath` to point to the Stage 1 output file you want to analyse. Run the script to fit distributions and save results.

---

## Repository Structure

```
nSRdist_code/
├── calcData.m              # Stage 1 main script
├── fitData.m               # Stage 2 main script
├── Functions/              # All helper functions
├── ResultsScripts/         # Scripts for producing results figures
├── PlottingScripts/        # Standalone plotting scripts
├── SyntheticTests/         # Synthetic data validation scripts
├── OtherFiles/             # Miscellaneous and exploratory scripts
├── Lin2014Code/            # Reference code from Lin et al. (2014)
├── BchronFolders/          # Bchron inputs and outputs (generated, not tracked in git)
├── DataSheets/             # Excel metadata spreadsheet
└── Results/                # Saved .mat output files (not tracked in git)
```

---

## Glossary

| Term | Meaning |
|---|---|
| NSR / nSR | Normalised sedimentation rate — sedimentation rate divided by the mean rate for that core |
| SR | Sedimentation rate (cm/kyr) |
| MSPF | Monospecific planktonic foraminifera — the material type used for radiocarbon dating |
| DDD | Doubly-dated depth — a sediment depth with more than one radiocarbon date |
| Scenario | One possible age-depth sequence, constructed to avoid doubly-dated depths and age reversals |
| WA / WA2022 | The Mulitza et al. (2022) World Atlas of foraminiferal isotope data |
| BMode | NSR estimated using the mode of the Bchron age-depth model |
| BMedian | NSR estimated using the median of the Bchron age-depth model |
| BSamp | NSR estimated by sampling from the Bchron posterior |
| RSR | Random Sampling sedimentation Rate — NSR estimated by sampling from calibrated radiocarbon age PDFs |
| MLN | Mixed Log-Normal distribution (2-component) |
| LN | Log-Normal distribution |
| invGam | Inverse Gamma distribution |
| MSI | Mean Spacing Index — a measure of the resolution of radiocarbon dates in a core |
| BIGMACS | A paleoclimate age model that uses the NSR distribution as a prior |
| CFR | Confirmed (no) reversals — a flag used internally during scenario construction |
| S | The settings structure passed throughout the pipeline |
| fitS | The fitting settings structure used in Stage 2 |

---

## Key References

- Lin, L., Khider, D., Lisiecki, L. E., & Lawrence, C. E. (2014). Probabilistic sequence alignment of stratigraphic records. *Paleoceanography*, 29(6), 976–989. https://doi.org/10.1002/2014PA002713
- Mulitza, S., et al. (2022). World Atlas of late Quaternary foraminiferal oxygen and carbon isotope ratios. *Earth System Science Data*, 14, 2553–2611. https://doi.org/10.5194/essd-14-2553-2022
- Haslett, J., & Parnell, A. (2008). A simple monotone process with application to radiocarbon-dated depth chronologies. *Journal of the Royal Statistical Society: Series C*, 57(4), 399–418. https://doi.org/10.1111/j.1467-9876.2008.00623.x (Bchron)
- Lougheed, B. C., & Obrochta, S. P. (2019). A rapid, deterministic age-depth modelling routine for geological sequences with inherent depth uncertainty. *Paleoceanography and Paleoclimatology*, 34(1), 122–133. https://doi.org/10.1029/2018PA003457 (MatCal)
