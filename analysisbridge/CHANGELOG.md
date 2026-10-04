# Changelog

## 1.2.0

- Added streamed scientific reduction for `.csv.gz` and `.tsv.gz` files.
- Added automatic inference of run ID, metadata, timestamps/elapsed-hours, and scientific variable roles.
- Added run, monthly, seasonal, daily, extreme-event, threshold-spell, wetting/drying, rainfall-response, matched-duration, ensemble, and baseline-change outputs.
- Added fallback time-axis reconstruction from `elapsed_hours` plus a `YYYY-YYYY` period.
- Added scientific-output counts to `MASTER_SUMMARY.json`.
- Added scientific outputs to upload packs while excluding giant source `.gz` files.
- Preserved source-only run fingerprints and SHA-256 provenance.
- Fixed duplicate run-ID metadata columns in scientific tables.
- Improved compressed-tabular file recognition.

## 1.1.0

- Initial Windows-first GUI/CLI release.
- ZIP extraction, hashing, duplicate detection, structured-file summaries, catalog database, audit, and upload packs.
