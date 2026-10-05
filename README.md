# Heritage Moisture & Salt Lab

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23151541.svg)](https://doi.org/10.5281/zenodo.23151541)

**Heritage Moisture & Salt Lab** is a reproducible computational framework for climate-to-material exposure modelling in porous archaeological heritage.

The software implements a reduced one-dimensional, thickness-averaged framework for coupled heat, moisture, and dissolved sulfate transport and is designed for auditable multi-decadal historical and future-climate simulations.

The framework was developed and demonstrated through a 21-case climate ensemble for **Mohenjo-daro, Pakistan**, comprising one historical ERA5-Land case and twenty future NEX-GDDP-CMIP6 simulations.

---

## Software DOI

The publication release **v1.5.0** is permanently archived on Zenodo:

**DOI:** [10.5281/zenodo.23151541](https://doi.org/10.5281/zenodo.23151541)

**GitHub repository:**  
https://github.com/noman63000/Heritage-Moisture-Salt-Lab

**Author:** Muhammad Nouman Akhtar  
**ORCID:** https://orcid.org/0000-0002-3367-1607

---

## Scientific scope

Heritage Moisture & Salt Lab provides an auditable climate-to-material modelling workflow for investigating how long-term climatic forcing can modify:

- material temperature;
- thermal-extreme exposure;
- moisture state;
- water transport;
- dissolved sulfate transport;
- near-surface sulfate concentration;
- sulfate redistribution with depth;
- water and sulfate mass balances;
- seasonal exposure behaviour;
- cold-state and pore-water/ice phase conditions;
- differences between historical and future climate scenarios.

The production implementation is a **reduced 1-D vertical, thickness-averaged transport framework**.

It is intended primarily for:

- heritage-science research;
- archaeological conservation studies;
- climate-change exposure assessment;
- long-duration scenario modelling;
- reproducible computational experiments;
- numerical verification and sensitivity analysis.

---

## What the software does not claim

The software should **not** be interpreted as a fully calibrated predictive digital twin of Mohenjo-daro.

The present publication version does not explicitly resolve:

- mixed-electrolyte thermodynamics;
- separate ionic transport for all individual salts;
- cation exchange;
- crystallization pressure;
- pore clogging;
- mechanically induced cracking;
- material loss;
- nucleation or supercooling hysteresis;
- detailed frost-damage mechanics;
- snow accumulation;
- groundwater-flow dynamics;
- deterministic conservation-damage prediction.

The Mohenjo-daro application therefore reports **material exposure and transport behaviour**, rather than direct prediction of physical damage.

---

# Publication release

This repository corresponds to:

**Heritage Moisture & Salt Lab v1.5.0**

The publication bundle consolidates the software lineages used in the Mohenjo-daro study, including:

- Heritage Moisture & Salt Lab scientific core;
- numerical updates developed during model verification;
- production multi-case execution workflow;
- checkpoint and restart functionality;
- cold-state stability checks;
- equilibrium water/ice phase extension;
- AnalysisBridge v1.2.0;
- numerical regression tests;
- standardized benchmark material;
- climate-disaggregation validation;
- lower-boundary sensitivity analysis;
- frozen publication configuration;
- climate-input manifests and QA records;
- provenance and reproducibility documentation.

The exact archived publication release is available at:

https://doi.org/10.5281/zenodo.23151541

---

# Repository structure

```text
Heritage-Moisture-Salt-Lab/
│
├── analysisbridge/
│   └── AnalysisBridge scientific reduction and audit workflow
│
├── app/
│   └── Local application/server and browser-interface components
│
├── benchmarks/
│   └── HAMSTAD_Benchmark_2/
│       ├── benchmark implementation
│       ├── analytical/reference solution
│       └── comparison outputs
│
├── climate_processing/
│   └── Author-developed climate preparation and
│       daily-to-hourly reconstruction workflows
│
├── configs/
│   ├── frozen Mohenjo-daro publication configuration
│   ├── locked model input workbook
│   └── exact material-definition files
│
├── data_manifests/
│   ├── forcing manifests
│   ├── forcing QA
│   ├── period summaries
│   └── source-data identities
│
├── docs/
│   ├── methodology
│   ├── verification documentation
│   ├── source-data information
│   └── testing records
│
├── execution/
│   └── hlab_performance/
│       ├── long-duration execution
│       ├── checkpoint/resume
│       ├── safety checks
│       ├── cold-state handling
│       └── phase-extension components
│
├── src/
│   └── hmsl/
│       ├── heat transport
│       ├── moisture transport
│       ├── sulfate transport
│       ├── material functions
│       ├── boundary conditions
│       ├── climate forcing
│       ├── project configuration
│       └── workflow components
│
├── tests/
│   └── numerical and software regression tests
│
├── validation/
│   ├── disaggregation/
│   ├── boundary_sensitivity/
│   └── event_consistency/
│
├── .gitignore
├── CHANGELOG.md
├── CITATION.cff
├── LICENSE
├── README.md
├── SHA256SUMS.txt
├── environment.yml
└── requirements.txt
