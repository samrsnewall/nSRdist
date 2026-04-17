# DataSheets

This folder contains the metadata spreadsheets and ancillary data files that control which sediment cores and radiocarbon dates are used in the analysis pipeline.

These spreadsheets are necessary because the Lin et al. (2014) and Mulitza et al. (2022) do not provide information about the material dated for every radiocarbon age, hence this sheet tells the code which data to use.

It also facilitates the excluding of some data manually if desired. For example, we manually removed outliers that we considered obvious to reduce the burden on the reversal choosing procedure used in the analysis pipeline. We also removed ages that were very far apart to help reduce the spread of the durations for which SR is measured across.

---

## datasheet.xlsx

This is the metadata spreadsheet read by `calcData.m` (via `extract3.m`). It contains one row per sediment core across the combined World Atlas 2022 (WA2022) and Lin et al. (2014) databases — 342 cores in total — regardless of whether a core is ultimately selected for analysis.

The spreadsheet serves two roles: (1) it determines *which cores* enter the analysis through a set of flag columns, and (2) it records *which individual radiocarbon dates* within each core should be excluded, and why.

### Core identification and geography

| Column | Used by pipeline | Description |
|---|:---:|---|
| `CoreName` | ✓ | Sediment core identifier (e.g. `RC13-228`). Must match the filename in the World Atlas data folder or the Lin2014 BchronInput folder. |
| `LongitudeDec` | ✓ | Longitude in decimal degrees (negative = West) |
| `LatitudeDec` | ✓ | Latitude in decimal degrees (negative = South) |
| `WaterDepthM` | ✓ | Water depth at the core site (metres) |
| `Basin` | ✓ | Ocean basin. Values: `Atlantic Ocean`, `Pacific Ocean`, `Indian Ocean`, `Mediterranean Sea`, `Red Sea`, `Bay of Bengal`, `Sea of Japan`. Used to apply different latitude limits to Atlantic vs. other basins. |
| `References` | — | Bibliographic reference(s) for the core's radiocarbon data |
| `Filename` | — | Filename of the source NetCDF file in the World Atlas database |
| `NumDates` | — | Total number of radiocarbon dates in the source file |

### Core selection flags

These columns determine whether a core is included in the analysis. Cells left blank (NaN) mean the flag does not apply.

| Column | Used by pipeline | Description |
|---|:---:|---|
| `UseChoice` | — | Working notes on whether a core was assessed for inclusion. Values include `PF` (planktonic foraminifera), `MSPF` (mixed-species planktonic foraminifera), `0` (excluded), `?` (uncertain). Informational only; does not control pipeline behaviour. |
| `DataLocations` | — | Notes on where radiocarbon data for this core were found. Values include `WA` (World Atlas), `Lin2014`, `WA, Lin2014`, `New`. Informational only. |
| `UseComment` | — | Free-text notes explaining why a core was or was not included |
| `WAuse` | ✓ | `1` = include this core using the World Atlas data source. Blank = do not include via this route. |
| `PForMSPF` | — | Notes whether the core's dates are pure planktonic foraminifera (`PF`) or mixed-species planktonic foraminifera (`MSPF`). Informational; does not control pipeline behaviour directly. |
| `Lin2014` | ✓ | `1` = this core is in the Lin et al. (2014) dataset and should be included when `S.useLin` is true. |
| `Lin2014Keep` | ✓ | `1` = this core is in both datasets but should be included using Lin2014 data even when running in WA mode. `0` = present in Lin2014 but excluded. |

### Date exclusion columns (World Atlas cores)

These columns list radiocarbon dates to exclude from WA cores before analysis. Lab IDs are stored as comma-separated strings (e.g. `KIA554, KIA4110`). Depth values are in cm.

| Column | Used by pipeline | Description |
|---|:---:|---|
| `Materials` | — | Description of the foraminifera species or material type used for each date. Informational. |
| `NonPlanktonicIDs` | ✓ | Lab IDs of dates on non-planktonic material (e.g. benthic foraminifera, bulk carbonate) to be excluded |
| `ReversalIDs` | ✓ | Lab IDs of dates that produce an age reversal and should be removed |
| `ReversalDepths` | ✓ | Depths of reversal dates, used as an alternative to lab IDs when IDs are unavailable |
| `AgeGapIDs` | ✓ | Lab IDs of dates that produce an unusually large age gap. Excluded only when `S.removeLargeGaps` is `true`. |
| `AgeGapDepths` | ✓ | Depths of age-gap dates. Used alongside `AgeGapIDs` when `S.removeLargeGaps` is `true`. |
| `MiscRemovalIDs` | ✓ | Lab IDs excluded for miscellaneous reasons not covered by the categories above |
| `Comment` | — | Free-text notes on the core or its exclusions |

### Date exclusion column (Lin2014 cores)

| Column | Used by pipeline | Description |
|---|:---:|---|
| `Lin2014ManualExclude` | ✓ | Lab IDs to manually exclude from Lin2014 cores. Applied when `S.modifyLin2014Data` is `true`. Comma-separated. |

### MSPF auxiliary columns

These columns support a separate species-selection workflow and are **not read by the main pipeline** (`calcData.m` / `extract3.m`). They are retained here for reference alongside the `COPYcorechoices_MSPF*.xlsx` files.

| Column | Description |
|---|---|
| `MSPF` | `1` = a usable mixed-species planktonic foraminifera record exists; `0` = does not; `?` = uncertain |
| `MSPF_species` | The foraminifera species used for the MSPF analysis |
| `MSPF_includeID` | Lab IDs of MSPF dates to include (often `all`) |
| `MSPF_includeDepth` | Depths of MSPF dates to include when lab IDs are unavailable |
| `MSPF_ReversalIDs` | Lab IDs of MSPF dates to exclude as reversals |
| `MSPF_ReversalDepths` | Depths of MSPF reversal dates |
| `MSPF_AgeGapIDs` | Lab IDs of MSPF dates to exclude as large age gaps |
| `MSPF_AgeGapDepths` | Depths of MSPF age-gap dates |

