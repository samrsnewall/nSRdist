# Lin2014Code copy

This folder contains of the MATLAB files and some associated data files from the Lin et al. (2014) study that were used to better understand their methods. We include it here for reference and to support understanding and reproduction of results.

MATLAB files that were directly from the authors of Lin et al. (2014) are: `calc_14c_sedrates5.m`, `calc_sedrates_trans4.m`, `func_combine_files.m`, `read_radiocarbon.m`, `read_14c_calibrated.m`, and `read_14c_ranges.m`.

The following files were modified for debugging and to better understand the process: `sedrates_lin_and_lee.m`, `func_combine_files_copy.m`, 

The Input folder contains the radiocarbon data used by Lin et al. (2014). 

The Output folder contains the Bchron output for two cores. Only two cores are provided due to the large file size of the outputs. These are 

---

## Scripts

### `calc_14c_sedrates5.m`
The main Lin et al. (2014) analysis script. Loops over a hard-coded list of 55 sediment cores, reads each core's radiocarbon depths and Bchron-derived calibrated ages from `Input/` and `Output/`, computes normalised sedimentation rates (SR divided by each core's mean SR), then pools all cores and produces a duration-weighted histogram of the SR ratio distribution.

### `calc_sedrates_trans4.m`
Computes sedimentation rate ratios in 1-kyr increments and calculates their autocorrelation. Uses the same core list as `calc_14c_sedrates5.m`.

### `func_combine_files.m`
Helper function used to combine input data from cores that have been split across multiple sub-segment files (e.g. `MD07-3076a`, `MD07-3076b`, ...).

### `read_radiocarbon.m`
Helper function. Reads a Bchron-formatted radiocarbon input file from `Input/<core>_radiocarbon.txt` and returns the sample depths and lab IDs.

### `read_14c_calibrated.m`
Helper function. Reads the calibrated age ranges file from `Output/<core>_radiocarbonCalibratedRanges.txt` and returns the 2.5%, 50%, and 97.5% calibrated ages for each date.

### `read_14c_ranges.m`
Helper function. Reads the Bchron age model file from `Output/<core>_radiocarbonranges.txt` and returns depths and the 50% (median) Bchron model age at each depth.

### `sedrates_lin_and_lee.m`
Edited (and debugged) version of `calc_14c_sedrates5.m` to better understand the process of Lin et al. (2014).

### `func_combine_files_copy.m`
Edited (debugged) version of `func_combine_files.m`.


---

## Data files

### `Lin2014_dts.mat`

NSR, duration, and length for the SR estimates after running `sedrates_lin_and_lee.m`. Used in `results_Newall`.

---

## Input/

Contains 252 Bchron-formatted radiocarbon input files, one per core or core variant. Files are tab-delimited with a header row and the following columns:

| Column | Description |
|---|---|
| `id` | Laboratory ID for the radiocarbon date |
| `Age` | Conventional radiocarbon age (yr BP) |
| `Error` | 1-sigma uncertainty on the radiocarbon age (yr) |
| `Depth` | Sample depth in core (cm) |
| `Thickness` | Sample thickness (cm); typically 0 |
| `Outlier1` | Prior outlier probability passed to Bchron |
| `Outlier2` | Secondary outlier parameter |
| `Type` | Sample type flag; 1 = included |

### File naming conventions

Most cores have a single input file (`<CoreName>_radiocarbon.txt`). Cores with doubly-dated depths or alternative date selections appear as variants:

| Suffix | Meaning |
|---|---|
| *(none)* | Standard input file for the core |
| `a`, `b`, `c`, … | Sub-segments of a core|
| `r` | Variant with a different reservoir age correction applied |
| `rescorr2` | Variant with a second reservoir correction applied |
| `all` | Variant including all available dates (e.g. `KNR31-GPC5all`) |
| `young`, `complete`, `comp` | Other date-selection variants specific to individual cores |

---

## Output/

Contains Bchron output files that were used in Lin et al. (2014). These are included to highlight what information was extracted from Bchron during the Lin et al. (2014) study, as it influenced the method used. Output files for only two cores (`AHF16832` and `V19-30`) are provided, as providing all files would make for a large file size.

Each core produces six output files. Only two are actually read by `sedrates_lin_and_lee.m` (via `read_14c_ranges.m` and `read_14c_calibrated.m` respectively); the remaining four are not used in the analysis.

| File | Used By Script | Description (Best Guess) |
|---|:---:|---|
| `<core>_radiocarbonranges.txt` | ✓ | Bchron age model summary: depth, 2.5, 50, and 97.5 percentile of model ages (kyr BP) at each sampled depth (Sampled depths are clearly specified such that they represent the depths at which the median of the age-depth model intersects integer kiloyears) |
| `<core>_radiocarbonCalibratedRanges.txt` | ✓ | Calibrated age percentiles for each radiocarbon date: lab ID, 2.5, 50, and 97.5 percentile of calibrated radiocarbon ages (kyr BP) (these have not been influenced by Bchron)|
| `<core>_radiocarbonTrueDates.txt` | — | Bchron age ensemble at the depths of the radiocarbon dates specifically (ages in kyr BP) |
| `<core>_radiocarbonoutliers.txt` | — | Per-date outlier probabilities from Bchron: lab ID, `Prob1`, `Prob2` |
| `<core>_radiocarbonchrons.txt` | — | Bchron age ensemble: each row is one MCMC iteration, columns are the same depths as used in `<core>_radiocarbonranges.txt` |
| `<core>_radiocarbonpars.txt` | — | ? |
