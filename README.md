# Heritage Moisture & Salt Lab

Reproducible climate-to-material exposure modelling for heat, moisture and dissolved sulfate transport in porous archaeological heritage, developed and demonstrated with a 21-case historical/future climate ensemble for Mohenjo-daro, Pakistan.

## Scope

This repository contains the publication-facing software and verification materials used in the associated Mohenjo-daro study. The scientific model is a reduced one-dimensional, thickness-averaged heat-moisture-Na2SO4-equivalent transport framework. It is designed for auditable multi-decadal scenario experiments rather than as a site-calibrated structural-damage predictor.

It does **not** claim to resolve mixed-salt thermodynamics, crystallization pressure, cracking, pore clogging, or deterministic material loss. The Mohenjo-daro application uses a literature-derived historic-brick material proxy and explicit lower-boundary scenario assumptions; these limitations are part of the study design and are reported in the manuscript.

## Publication release

Repository publication bundle: **v1.5.0**.

The bundle consolidates the exact component lineages used in the study:

- Heritage Moisture & Salt Lab scientific core: v1.2.0, with v1.2.1 preparation-interface update;
- parallel/resumable execution layer from the verified all-cases workflow;
- final pore-liquid/equilibrium freeze-thaw execution extension (v1.4 lineage);
- AnalysisBridge v1.2.0 for output audit and scientific reduction.

A Zenodo DOI will be added to this README and `CITATION.cff` when the release is archived.

## Repository map

- `src/hmsl/` - heat, moisture, sulfate, material, climate, project and workflow scientific core.
- `app/` - local application/server and browser interface.
- `execution/hlab_performance/` - long-run worker, checkpoint/resume, safety, cold-state and phase-extension modules.
- `analysisbridge/` - AnalysisBridge v1.2.0 scientific reduction/audit tool.
- `tests/` - core and phase-extension regression tests.
- `benchmarks/HAMSTAD_Benchmark_2/` - standardized HAMSTAD Benchmark 2 adapter, reference comparison outputs and summary.
- `validation/disaggregation/` - leave-one-year-out validation of the daily-to-hourly climate reconstruction workflow.
- `validation/boundary_sensitivity/` - lower-boundary sensitivity calculations used in the paper.
- `validation/event_consistency/` - independent historical-event consistency evidence.
- `configs/` - frozen Mohenjo-daro configuration, locked input workbook and exact historic-brick material file used by the study.
- `climate_processing/` - author-developed future daily-to-hourly forcing preparation scripts.
- `data_manifests/` - source identities, climate-forcing QA and period summaries. Large climate and simulation data are not stored in GitHub.
- `docs/` - methods, testing documentation, source-data links and verification reports.

## Installation

Python 3.11 was used for the production study. A reproducible dependency set is provided in `requirements.txt` and `environment.yml`.

With Conda:

```bash
conda env create -f environment.yml
conda activate heritage-moisture-salt-lab
```

Or with a Python virtual environment:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```

## Verification

The study used staged numerical and software verification, including mesh and maximum-step checks, common historical spin-up, wet-end stress tests, water/sulfate conservation ledgers, regression tests, extreme forcing windows and cold/phase checks.

A standardized HAMSTAD Benchmark 2 adapter is included under `benchmarks/HAMSTAD_Benchmark_2/`. It verifies the conservative finite-volume/BDF moisture numerical pattern against the benchmark analytical solution; it does not constitute site-specific validation of the Mohenjo-daro material or boundaries.

## Climate forcing and large data

Original ERA5-Land and NASA NEX-GDDP-CMIP6 products are external datasets and are credited to their respective producers. This repository includes processing code and exact manifests/QA records, not a duplicate of every large external source file.

The full R00-R20 processed research dataset and large simulation outputs are intended for a separate persistent data archive. The software DOI and data DOI will be cross-linked in the publication record.

## Reproducibility boundaries

The repository supports reproducibility of the implemented computational workflow. It should not be interpreted as evidence of full predictive validation at Mohenjo-daro. Site-specific archaeological brick functions, time-resolved internal moisture/temperature measurements and mixed-ion chemistry remain important future validation/development needs.

## Citation

GitHub will read citation metadata from `CITATION.cff`. After Zenodo archives release `v1.5.0`, the software DOI should be added there and here before journal submission.

## License

MIT License. See `LICENSE`.
