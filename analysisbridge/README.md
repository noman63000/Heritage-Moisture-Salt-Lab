# AnalysisBridge 1.2.0

AnalysisBridge prepares very large research/model-output folders and ZIP archives for reliable analysis without manual extraction, hand-made manifests, or uploading giant hourly files.

## What is new in 1.2

Version 1.2 adds a streaming scientific-reduction engine for `.csv.gz` / `.tsv.gz` outputs. It reads the compressed table directly and creates compact analysis-ready derivatives while leaving the original file untouched.

For large hourly model outputs it can automatically generate:

- run-level descriptive statistics and approximate quantiles;
- monthly and seasonal summaries;
- exact top/bottom extreme events for key variables;
- threshold counts and longest temperature / rain / ice spells when those variables are detected;
- daily reductions from hourly data;
- wetting/drying-cycle diagnostics for moisture variables;
- rainfall-response windows when precipitation variables exist;
- duration-matched run summaries;
- duration-matched ensemble summaries by scenario/period;
- change-versus-historical-baseline tables when a historical/ERA baseline is detected;
- an inferred schema describing exactly which columns were treated as temperature, moisture, sulfate/salt, precipitation, ice, metadata, and time.

This is particularly useful for files such as `MASTER_HOURLY_KEY_METRICS.csv.gz`: the raw file can remain hundreds of MB locally while the scientifically useful derivatives are included in the much smaller ChatGPT upload packs.

## Core safety and provenance features

- Accepts folders, individual files, and ZIP archives.
- Recursively catalogs content and safely extracts ZIPs in an isolated workspace.
- Blocks unsafe ZIP paths and records corrupt/unreadable files.
- Calculates SHA-256 for every source file and flags exact duplicates.
- Summarizes CSV/TSV, JSON/JSONL, TXT, LOG, Markdown, YAML/XML-like text files.
- Recognizes compressed CSV/TSV files instead of treating them as opaque binaries.
- Creates `INVENTORY.csv`, `ISSUES.csv`, `MASTER_SUMMARY.json`, `FINAL_AUDIT.json`, and `analysis_catalog.sqlite`.
- Builds sequential upload packs in `CHATGPT_UPLOAD`.
- Never modifies original input files.
- Derived scientific files do not alter the source run fingerprint.

## Easiest Windows workflow

1. Extract the AnalysisBridge software ZIP once.
2. Double-click `run_analysis_bridge.bat`.
3. Click **Add Folder** or **Add Files / ZIPs**.
4. Select the original raw result folder or the original large result ZIP. Do **not** manually extract research/model ZIPs.
5. Leave **Scientific reduction of large .csv.gz outputs** checked.
6. Leave **Fast** mode selected unless you specifically need deep validation of huge uncompressed JSON/CSV files.
7. Click **Build Analysis Packs**.
8. Open the generated `CHATGPT_UPLOAD` folder.
9. Upload `ANALYSIS_INDEX.zip` first, then `ANALYSIS_PACK_001.zip`, `ANALYSIS_PACK_002.zip`, etc.

For the Mohenjo-daro run, select the original `mohenjo_daro_FINAL_RESULTS_....zip` (or its parent results folder). AnalysisBridge 1.2 will stream the internal `MASTER_HOURLY_KEY_METRICS.csv.gz` automatically.

## Important note about time axes

If an hourly table has a real timestamp/date column, AnalysisBridge uses it. If it has only `elapsed_hours` plus a period such as `2040-2060`, AnalysisBridge reconstructs an analysis time axis beginning January 1 of the first period year. That reconstruction is clearly documented in the generated schema and may not exactly reproduce non-Gregorian calendars.

## Scientific outputs

Each compressed scientific table gets its own folder inside `SCIENTIFIC_REDUCTION/`. Depending on detected columns, typical files include:

- `HOURLY_SCHEMA.json`
- `HOURLY_RUN_STATS.csv`
- `HOURLY_MONTHLY_STATS.csv`
- `HOURLY_SEASONAL_STATS.csv`
- `HOURLY_EXTREME_EVENTS.csv`
- `HOURLY_THRESHOLD_COUNTS.csv`
- `HOURLY_TOP_THRESHOLD_SPELLS.csv`
- `HOURLY_DAILY_REDUCTION.csv`
- `WETTING_DRYING_CYCLES.csv`
- `RAIN_RESPONSE_EVENTS.csv`
- `MATCHED_DURATION_RUN_SUMMARY.csv`
- `MATCHED_DURATION_ENSEMBLE_SUMMARY.csv`
- `MATCHED_DURATION_CHANGE_VS_BASELINE.csv`
- `SCIENTIFIC_REDUCTION_REPORT.json`

Names use `DAILY_`, `PHASE_`, or `TABLE_` instead of `HOURLY_` when appropriate.

## Build a one-file Windows EXE

Double-click `build_windows_exe.bat` on a Windows machine with Python installed. It installs PyInstaller and creates:

`dist\AnalysisBridge.exe`

After that, the EXE can be used like a normal standalone Windows application on that computer.

## Command-line mode

```text
python analysis_bridge.py "D:\my_results" -o "D:\AnalysisBridge_Output" --mode fast --pack-mb 180
```

Disable scientific reduction only if required:

```text
python analysis_bridge.py "D:\my_results" -o "D:\AnalysisBridge_Output" --no-science
```

## Privacy

All preparation and scientific reduction happens locally. AnalysisBridge does not upload data itself. Only the files you manually upload from `CHATGPT_UPLOAD` leave your computer.
