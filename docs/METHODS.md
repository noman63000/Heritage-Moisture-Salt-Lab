# Methods and limitations - version 1.2.0

## Status and purpose

This is an independently implemented **reduced, one-dimensional, thickness-averaged heat/moisture/Na2SO4-equivalent transport model**. It is not DELPHIN, a two-dimensional wall/soil model, a complete salt thermodynamic model, or a calibrated damage predictor. Successful numerical integration and local mesh QA are not field validation.

Version 1.2.0 changes the internal liquid-face numerical discretisation, not the source material or the reviewed 1.1 exposed-face drainage assumption. It is not an acceptance-tolerance patch. Source tables, physical scalar inputs, initial QA state, selected QA windows, 0.10 m band support, normalisation floors and the 5% criterion are unchanged. Previous-scheme QA, spin-up and freeze records cannot certify this numerical revision; the installer archives them and preserves original data and prior results.

## 1. Source identity and input conventions

The reported repair uses the user's actual `project.json` and original `18_15_IBEET_3.m6`, with the recorded ERA5-Land climate checksum. The source material is a **Latvian historic-brick analogue**, not a measured Mohenjo-daro specimen. Material scalar rounding is checked but original tables are not overwritten, refitted or extrapolated. The user-provided source is now import-tested; this is not a fresh verification of the remote archive download.

The shared table coverage for the uploaded desorption material is **theta = 0.0000143836 to 0.349795 m3/m3**. Its effective-saturation scalar is **0.35**. The gap is **0.000205 m3/m3, or 0.0585714% of effective saturation**. Open porosity is a different quantity (0.38016). None of these quantities is silently substituted for another. A non-closed model with a wet-end coverage gap greater than 0.1% of effective saturation is rejected: that is an explicit implementation safeguard, not an experimental uncertainty estimate or universal physical criterion.

The supported forward/inverse retention pair is checked for consistency. The chosen desorption branch stays fixed; scanning hysteresis is not modelled. `pC` means log10(-capillary pressure / Pa). `lgKl` is converted to its pressure-driven liquid conductivity units (seconds), not read as diffusivity. Original vapor and thermal functions are used over their supported range.

Each climate record represents an hour [t,t+1). Fields are temperature in C, RH in percent, rain in mm/hour, solar and long-wave radiation in W/m2, wind in m/s. Each of the 21 original climate cases is separately checksummed and interpreted in its own calendar; no GCM averaging, repeated bias correction, missing-hour filling, or climate clipping is performed here. Calendar/schema success is separate from the model's physical applicability.

The sulfate input field uses **mg SO4/L**, not mg Na2SO4/L. Dividing by 96.06 g/mol gives sulfate-equivalent moles per m3 liquid. All inventories in reports use a square metre of horizontal wall footprint, not a square metre of facade.

## 2. Geometry and unchanged user assumptions

The verification uses H = 1.5 m, thickness b = 0.279 m, 200 uniform base cells, two exposed broad faces and an exposed crown. The QA reference has 400 cells. These are the user's reduced-domain assumptions, not independently measured wall geometry. Uniform200 gives dz = 0.0075 m; uniform400 gives dz = 0.00375 m. A developer-only transformed-grid option exists but is NOT used in the supplied actual-site verification or recommended base configuration.

The user's optical values, exposure multipliers, face rain fraction 0.05, crown rain, lower supply 0.2 mm/day, source sulfate 80.4 mg/L, reference D = 2.91e-10 m2/s, pore-concentration convention, zero added saturation exponent and zero initial sulfate are unchanged. The source sulfate is a regional proxy; groundwater depth is not converted into this imposed supply. Different sites require their own justified values.

## 2A. Version 1.2.0 liquid-face discretisation

The selected homogeneous brick supplies pressure retention p(theta) and log10 liquid conductivity log10 K(theta). The legacy face coefficient was the harmonic mean of its two endpoint conductivities. Version 1.2 uses a source-integral mean:

    K_face = [integral from p_left to p_right of K(p) dp] / (p_right - p_left)
    q_face = -(K_face/rho_water) * [(p_right - p_left)/distance + rho_water*g]

The same mean multiplies both the capillary-pressure gradient and the gravity term. This preserves zero flux for a discrete hydrostatic pressure gradient. The equal-pressure limit is continuous. Internal fluxes are still shared with opposite signs by adjoining finite volumes; water and advected salt are not added or removed by changing the average.

For this source, pC(theta)=log10(-p/Pa) and log10 K(theta) are piecewise linear. Their union of tabulation knots partitions the integral into segments. On each segment, with slopes a=dpC/dtheta and b=dlogK/dtheta, the integrand K dp/dtheta is -ln(10)*a*10**(pC+logK); its integral is evaluated analytically with stable expm1 limits. No source point is fitted, rescaled or replaced. Coefficients are evaluated at the existing source bounds; independently accepted state bounds are unchanged. This numerical choice is motivated by resolving steep constitutive changes across a cell, not by changing the measured liquid conductivity.

The implementation deliberately accepts only pressure-based retention with lgKl(Theta_l) for this scheme. It stops rather than silently adapting RH-based retention or a diffusivity-only source. Such projects retain the documented harmonic_legacy setting until a compatible implementation is independently tested. A different material needs its own QA; the Mohenjo checks do not establish universal accuracy.

The configuration setting is numerics.liquid_face_scheme = source_integral. The user setting, source file bytes and actual implementation source bytes are hashed. The legacy option is retained for explicit comparisons, not to allow cross-scheme reuse of a certificate. The installer selects source_integral only for projects requiring brick identity 18_15.

### Why the spin-up hypothesis was rejected

The supplied old spin-up NPZ matched its metadata checksum and exact saved configuration. Applying it to the old 200-cell method and conservatively splitting its cells for a common-state 400-cell sensitivity test did not resolve R01, R03, R07 or R19. Those diagnostics still exceeded the 5% gate. That test was not an independently converged fine-mesh spin-up and was not the actual transient state decades into the future climate. It was used only to test the initialization hypothesis. No QA PASS was created from it.

The final revision retains the original, deliberately abrupt local QA tests: generic RH70%, T25 C and specified initial salt at the beginning of every chosen 72h interval, independently on each mesh. The new source-integral method is tested against those SAME cases. Changing initialization to obtain a more favourable pass is not part of this release.

## 3. Conservative finite-volume water balance and retained v1.1 wet-end closure

Cell state is theta, T and total sulfate-equivalent moles n per m3 bulk. Internal liquid face flux is the pressure-gradient-plus-gravity law using the source-integral mean described in Section 2A for this release's selected brick. The legacy harmonic mean is an explicit alternate setting, not the selected method. Vapor flux remains driven by vapor pressure. Paired internal fluxes cancel in the domain inventory.

Let F_i be the water rate from internal divergence plus captured rain minus vapor exchange before the new drainage term. The revised equation is:

    dtheta_i/dt = F_i - s_i

At an exposed cell's supported wet endpoint theta_hi:

    s_i = max(F_i, 0) if that endpoint constraint is active;
    s_i = 0 otherwise.

`s_i` is liquid volume leaving per bulk volume per second. An event-driven active set locates first wet-end contact, holds the cell only while net influx is positive, and releases it when net influx turns negative. The corresponding runoff/seepage ledger receives rho_water*s_i*cell_width. Dissolved sulfate leaves at c_i*s_i*cell_width. This is **not post-step clipping or hidden water deletion**.

This closure is an **explicit reduced free-draining exposed-face assumption at the supported table endpoint**. It is not a measured Mohenjo drainage law, a saturated positive-pressure solve, a proof that the table endpoint equals pore saturation, or a resolved film/runoff route. Exposed salt can leave with seepage. Full pressure/ponding and spatial drainage physics would require a different model. Users must explicitly review this assumption. The tests establish its numerical implementation, not its adequacy for every wall.

The earlier illustrative rainfall and basal-supply capture factor `clip(1-(theta/theta_hi)^4,0,1)` remains; it is not calibrated infiltration data. Unabsorbed rain is accounted for separately. Condensation and internal influx can activate the new wet-end constraint even when rain capture is almost zero. Closed/internal cells are not given fictitious exposed drainage.

At the supported dry endpoint, outgoing donor water fluxes and evaporation are limited rather than allowing water removal from an unsupported state. A still-unsupported state is rejected, not mass-clipped.

## 4. Heat approximation

    (rho_brick*cp_brick + rho_water*cp_water*theta) dT/dt
        = conduction divergence + exposed sensible/radiative exchange - Lv*evaporation

The same source thermal function, hc = 5.7 + 3.8*wind, declared Lewis-number-1 air-side analogy, b/4 lumped material resistance, sky-view and solar-projection factors are retained. These are reduced model assumptions, not measured local transfer functions. Crown exchange is distributed into the top cell; broad-face exchange is distributed through the thickness-averaged column.

Water-carried sensible enthalpy (including drainage enthalpy), hydrate water/heat and ice are absent. The new exact water and sulfate ledgers do not claim full thermodynamic enthalpy conservation.

## 5. Salt and unsupported interpretations

    dn/dt = -div(q_liquid*c_upwind - mobility*grad(c)) - c*seepage

For the selected pore convention, mobility is theta*D, with no added empirical saturation power. A bulk convention is selectable only with explicit evidence/assumption review. The term 'effective diffusion coefficient' alone does not identify this convention. No activation energy is invented.

The user's project has no solubility file. Consequently c = n/theta is a **dissolved-equivalent transport indicator**, not a chemically valid concentrated-solution composition at all concentrations. This release does not provide validated crystallized mass, mirabilite/thenardite cycles, supersaturation pressure, mixed-ion chemistry, NaCl transport or physical damage. Optional user-sourced single-salt solubility partitioning still lacks full activities, hydration water and kinetics. The zero initial salt condition is a scenario starting inventory, not a statement that historical masonry was salt-free.

## 6. Integration, bounds and independently recorded exchanges

Each climate hour is a constant-forcing BDF initial-value interval; accepted physical states are continuous across hours. Event restarts handle active wet constraints. The maximum step remains the user's 1800 s; the adaptive solver usually takes smaller steps. Four smaller-step retries remain for numerical failures, not missing-physics errors.

The effective integration relative tolerance is **min(requested rtol, 3e-8)**. The uploaded project requests 1e-6; the repair tightens rather than relaxes it. Physical absolute tolerances are theta 1e-10, temperature 1e-6 C and bulk salt 1e-11. This tightening was needed by the independent diffusion benchmark after hourly solver restarts. QA's 5% criterion is a different tolerance and is unchanged.

A callable sparse finite-difference Jacobian exploits the physical nearest-neighbor pattern. Local external-exchange derivatives are assembled separately for global accounting rows; those rows no longer force a dense finite-difference grouping. Ordinary laptop launch defaults to one BLAS thread to avoid oversubscription; `HMSL_NUM_THREADS` explicitly overrides this.

Nine independently integrated ledgers record net water, net salt, seepage water, seepage salt, unabsorbed rain, accepted rain, net evaporation (negative for condensation), basal water and basal salt. Stored inventory changes are checked against external flux integrals in every committed chunk. These are not reconstructed from inventory changes and therefore can detect bookkeeping defects. Small closure error demonstrates accounting consistency, not site validity or continuum convergence.

Accepted hourly states are checked against supported theta bounds with 1e-9 allowance and salt nonnegativity with 1e-10 roundoff allowance. Tiny IEEE underflow (around -1e-322) is numerical roundoff, not a modeled negative salt pool. No larger unphysical state is silently accepted. A solved material cell below 0 C, checked on accepted internal BDF states, stops transport as outside model scope.

## 7. Local numerical QA - what a PASS means

For R00, selected 72-hour windows start at rows 0, 3461 and 110674. Each starts from the same stated RH70%, T25 C and zero salt, not an interpolated coarse solution or the state from a previous year. The wet/hot dates label climate windows; these QA windows are not historically conditioned predictions of those dates.

Metrics are time-mean theta, final water inventory, time-mean temperature, final sulfate inventory and the peak dissolved-equivalent concentration over **fixed 0.10 m vertical bands**. The last metric averages liquid concentration over a fixed vertical support; it is not a resolved outer facade layer. The single-cell maximum remains diagnostic only. Width-weighted cell intersections are used across band edges. No metric, band width, initial test state or normalisation floor was changed by version1.2.

Errors are normalized by max(abs(reference), floor), where floors in metric order are [0.01, 1, 1, 1e-5, 0.001] in their respective units. Thus small theta changes may be floor-normalized rather than pure percentage relative error. All gate metrics must be <=0.05. Mesh200 vs400 is tested first; then max-step1800/900/450s at200. A zero reported difference can occur because all three maximum step caps are above the adaptive solver's chosen steps; it does not prove an independently resolved temporal truncation error or independence from hourly forcing.

The release verification report records the actual version1.2 case-by-case QA results and its independently executed tests. A local pass does not validate all hours of a multi-decade trajectory, every spatial maximum, different boundary conditions, field damage or every climate extreme. Additional 100-cell, 800-cell and tighter-tolerance experiments are diagnostic evidence, not an automatic continuum-error certificate. Further independent benchmark/field validation and longer conditioned comparisons remain necessary for physical predictions.

Completed QA variants are source/config/code/window hashed. Changed code or input cannot reuse an incompatible prior pass. True download links expose the latest QA record and error log; a heartbeat is explicitly distinguished from confirmed hourly progress.

## 8. Spin-up, full runs and locking

Water/heat spin-up repeats the complete first model-calendar year, checking maximum cell changes against theta1e-4 and temperature0.05 C for the uploaded project. Salt is set separately to the specified initial inventory after spin-up; no salt steady state is claimed. Year-start and current-state hashes are saved with every24h committed chunk, so cancellation/server restart can resume without losing entire completed years.

Frozen runs reuse the resulting baseline state. This common-initial-state comparison is a study-design assumption; long future scenario windows are not automatically equilibrated to their period-specific climate. Freeze is an input/code/state reproducibility lock, not validation. Old controls are archived by the update, not rewritten to say PASS.

Full transport has chunk checkpoints and code/source identity checks. A test-window report is not a full-period R00 result. Time-normalize inventories/flux totals when comparing unequal historical/future durations and calendars. Do not label all 21 cases complete unless each has its own successful full-run output and applicable physics.

## 9. Important full-ensemble compatibility result

All21 input files pass schema/calendar/checksum checks, but **R01, R13, R14, R16, R17, R18, R19 and R20 contain subzero air temperatures**. Air below0 C is the existing conservative ice-free scope guard, not proof that the masonry actually froze. Version1.1 identifies this before a long batch. Material cooling below0 C also stops at runtime, even for positive-air cases.

No cold hours or GCMs were removed, no temperature was replaced with zero and no freezing property was invented. These eight complete cases require a documented change in model scope or a freeze-thaw-capable model before claiming the original full21 experiment. Running a warm72h window inside an otherwise unsupported case does not establish full-case compatibility.

## 10. Sources, implementation documentation and evidence

Primary documentation relevant to units and integrator implementation:

- SciPy solve_ivp: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
- DELPHIN material editor manual: https://bauklimatik-dresden.de/delphin/2nd/doc/DELPHIN6_1_MaterialEditor_en.pdf
- Source archive identity retained from input selection: https://zenodo.org/records/5656966
- Historic-brick data publication: https://doi.org/10.3390/buildings13102529

These references do not validate our new reduced drainage closure. Its equations, source endpoint gap and limitations are declared above. Numerical claims are supported by the actual logs, checksums and test scripts in `verification/repair_110`. Read `REPAIR_1_1.html` only for the earlier boundary-repair history. Current test scope and installation instructions are in `NUMERICAL_1_2.html` and `VERIFICATION_1_2.html`. Old methods are retained separately as superseded audit history.

## 11. Verification sources and release evidence

Read NUMERICAL_1_2.html and verification/numerical_120 for executed results. The old methods and older reports are retained as explicitly superseded history. This change is a numerical-method revision; new and legacy outputs must not be pooled as if generated by one identical frozen implementation.

Methodological reference: NASA/NPARC, Examining Spatial (Grid) Convergence, https://www.grc.nasa.gov/WWW/wind/valid/tutorial/spatconv.html . Systematic refinement is stronger than a bare two-grid comparison; neither is a substitute for physical validation. No formal Grid Convergence Index is claimed by this software.

Integrator reference: SciPy solve_ivp, https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html . max_step is an upper bound; adaptive steps may already be below all tested caps. Separate stricter-rtol checks in the release evidence address this limitation locally, not universally.

Kirchhoff-transform background: https://arxiv.org/html/2210.00260v2 . This is a mathematical background reference only. That paper's numerical method is not the finite-volume implementation here and does not validate this application.
