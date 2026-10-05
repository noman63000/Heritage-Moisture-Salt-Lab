# Heritage Moisture & Salt Lab

[![Software DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23151541.svg)](https://doi.org/10.5281/zenodo.23151541)
[![Dataset DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23151771.svg)](https://doi.org/10.5281/zenodo.23151771)

**Heritage Moisture & Salt Lab** is a reproducible computational framework for climate-to-material exposure modelling in porous archaeological heritage.

The software implements a reduced one-dimensional, thickness-averaged framework for coupled heat, moisture, and dissolved sulfate transport and is designed for auditable multi-decadal historical and future-climate simulations.

The framework was developed and demonstrated through a 21-case historical/future climate ensemble for **Mohenjo-daro, Pakistan**, comprising one historical ERA5-Land simulation and twenty future NEX-GDDP-CMIP6 simulations.

---

# Permanent research records

The software and research dataset associated with the Mohenjo-daro study are permanently archived as separate Zenodo research objects.

## Software

**Heritage Moisture & Salt Lab v1.5.0**

Zenodo DOI:

**https://doi.org/10.5281/zenodo.23151541**

GitHub repository:

**https://github.com/noman63000/Heritage-Moisture-Salt-Lab**

## Research dataset

**Mohenjo-daro R00–R20 Climate-Material Simulation Dataset v1.0.0**

Zenodo DOI:

**https://doi.org/10.5281/zenodo.23151771**

The dataset contains the complete R00–R20 simulation archive, processed scientific outputs, climate-forcing documentation, frozen inputs, numerical audit records, validation and sensitivity evidence, provenance information, supplementary tables, and complete hourly model outputs.

---

# Author

**Muhammad Nouman Akhtar**

University of Genoa, Italy

ORCID:

**https://orcid.org/0000-0002-3367-1607**

GitHub:

**https://github.com/noman63000**

---

# Scientific scope

Heritage Moisture & Salt Lab provides an auditable climate-to-material modelling workflow for investigating how historical and future climatic forcing can modify:

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
- numerical verification;
- sensitivity analysis;
- research-software development;
- climate-to-material exposure assessment.

---

# What the framework does not claim

The software should **not** be interpreted as a fully calibrated predictive digital twin of Mohenjo-daro.

The present publication version does not explicitly resolve:

- complete mixed-electrolyte thermodynamics;
- separate transport of all individual ionic species;
- cation exchange;
- crystallization pressure;
- pore clogging;
- mechanically induced cracking;
- deterministic material loss;
- nucleation and supercooling hysteresis;
- detailed frost-damage mechanics;
- snow accumulation;
- dynamically coupled groundwater flow;
- complete site-specific salt phase assemblages;
- deterministic conservation-damage prediction.

The Mohenjo-daro application therefore reports **material exposure, hygrothermal state, dissolved-sulfate behaviour, and transport response**, rather than direct prediction of physical deterioration.

---

# Publication software release

The publication software version is:

**Heritage Moisture & Salt Lab v1.5.0**

The frozen release is permanently archived at:

**https://doi.org/10.5281/zenodo.23151541**

The publication bundle consolidates the software lineages used in the Mohenjo-daro study, including:

- Heritage Moisture & Salt Lab scientific core;
- numerical updates developed during model verification;
- production multi-case execution workflow;
- checkpoint and restart functionality;
- long-duration execution utilities;
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

The GitHub repository can continue to evolve after publication, but the Zenodo v1.5.0 record preserves the exact software release associated with the study.

---

# Research dataset

The complete research dataset associated with the publication is archived separately as:

**Mohenjo-daro R00–R20 Climate-Material Simulation Dataset**

**Version:** 1.0.0

**DOI:**  
https://doi.org/10.5281/zenodo.23151771

The published Zenodo dataset contains approximately **862 MB** of research material, including:

- complete hourly R00–R20 simulation archives;
- core simulation results;
- annual summary tables;
- case-level summary tables;
- supplementary manuscript tables;
- detailed scientific data tables;
- prepared climate-forcing archive;
- climate-forcing manifest;
- climate-forcing QA;
- climate-period summaries;
- locked publication model input set;
- exact material-definition file used in the simulations;
- numerical validation evidence;
- lower-boundary sensitivity evidence;
- climate-disaggregation validation;
- provenance records;
- audit records;
- upload manifests;
- SHA-256 checksums;
- dataset documentation.

The research dataset DOI is:

**https://doi.org/10.5281/zenodo.23151771**

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
