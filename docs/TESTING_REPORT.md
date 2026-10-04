# Heritage Moisture & Salt Lab 1.1.0 - repair verification

**Test basis:** the user's actual project, original Latvia 18_15 material and supplied source climate. **Not site validation or a completed paper-results ensemble.**

## Main outcome

The old v1.0.4 solver's exact wet-weather failure was reproduced at climate-window row 11,2022-08-17 21:00,400 cells. Its endpoint moisture overshoot was the same as the user's diagnostic. The revised formulation completed that hour, the complete 72-hour QA window and an extended 168-hour wet-weather test.

R00 local QA now **passes with 200 base cells versus 400 reference cells**, the same 5% gate, same 0.10 m metric and same physical numeric inputs. This success follows a **declared change to the wet-end boundary formulation**, not merely increasing cell count or weakening a check. There is no guarantee of error-free execution for every data set or unsupported physical regime.

## What was repaired

The prior external rain/supply limiter did not stop internal redistribution or condensation from pushing a cell past its supported moisture tables. Smaller steps could not fix that equation-level inconsistency. The source upper table limit 0.349795 is slightly below its effective-saturation scalar 0.35. Neither was raised or substituted.

The new event-driven wet-end closure sends positive excess net inflow out of an **exposed** cell as liquid seepage, including its dissolved salt. This reduced free-draining-face assumption is explicitly reviewable and must be reported in methods. It is not a measured drainage law, a saturated pressure solution or a resolved runoff film. Water/salt leaving are independently recorded; no post-step moisture clipping is used. Material curves, initial inventory and climate were not changed.

Unchanged completed QA tasks were reused in an additional 12-variant cache test with integration deliberately disabled. The actual calculations had already been executed; no incompatible cache was accepted.

A local sparse Jacobian removes avoidable global finite-difference work. Integration's effective relative tolerance is tightened to min(user setting,3e-8); QA remains 5%. New output ledgers distinguish seepage water/salt, unabsorbed/accepted rain, vapor exchange and basal inputs. Long spin-up now has resumable, input-hashed chunk checkpoints. Stale browser tokens recover once without disabling authorization. Real QA/error download links and climate compatibility checks are present.

## Exact local QA result

| Climate window start | Maximum normalized 200-vs 400 difference | Gate |
|---|---:|---|
| 2010-01-01 | 1.071961% | PASS |
| 2010-05-25 | 4.934294% | PASS |
| 2022-08-17 | 3.857848% | PASS |

May is close to the 5% threshold; it is not a large safety margin. Normalization floors and all metric vectors are retained in `qa_actual_R00_200_400.json`. The largest wet-window difference is mean moisture, not a concealed concentration maximum.

Maximum-step 1800/900/450 s comparisons also pass. In this run they return identical reported metrics because the adaptive solver's steps are smaller than those ceilings. This is a limitation of a max-step comparison, not proof of zero temporal discretization error. A local 72-hour PASS does not establish convergence of the entire multi-decade ensemble or single-cell maxima. QA starts from the same RH/T/salt, not historical warm-start states.

## Regression and interface checks

**57 automated tests passed**, including the previous test suite, real-format import, water/heat/salt analytic or discrete-diffusion benchmarks, bound events, independent seepage accounting, initial-state rejection, dry/wet controls, Jacobian comparison, spin-up checkpoint/resume and missing-freezing safeguards. The earlier v1.0.4 suite passed 40 while the separate wetting reproduction failed; passing unit tests alone was insufficient.

**27 interface/API checks passed**. The real application served the actual project and material, rejected tokenless/foreign-origin/path-traversal requests, ran the synthetic 24-hour interface smoke test, produced the compatibility report, displayed the 200-cell input and QA download link, and recovered a deliberately stale token. No uncaught JavaScript errors were observed.

Chromium's direct localhost navigation was blocked by container policy. The actual frontend was therefore exercised with browser fetch calls bridged to the running local HTTP API, not a mocked solver. This is not a native Windows launch test.

## Actual-data stress runs

- Historical 2022-08-17 window: **168 hours at 400 cells completed**, including the original failing hour. No rejected-hour retries. Maximum hourly theta 0.349795. Maximum chunk water-closure residual approximately 1.64e-13 kg/m2 footprint; salt residual approximately 1.88e-18 mol/m2 footprint.
- Highest total-rainfall 72-hour future window, R16,2091-08-03: **completed at 400 cells**.
- Hottest mean-temperature 72-hour future window, R04,2100-07-19: **completed at 400 cells**.
- Coldest selected future window, R14: **correctly refused** by the ice-free scope guard. The raw test result retains its exception as `FAILED`; it is an expected scope rejection, not an accepted transport calculation.

Exact logs, input windows and solver statistics are supplied. These checks are integration stress tests, not future-case mesh convergence or full-period results. The warm R16 window does not make the full R16 file supported.

## Historical water/heat spin-up

**PASS** after 2 full annual cycles.

| Cycle | Max theta change (m3/m3) | Max temperature change (C) |
|---|---:|---:|
| 1 | 0.0214651096 | 17.0752075 |
| 2 | 4.68375339e-17 | 6.21724894e-15 |

The test uses the full 2010 climate year, 200 cells, the user's theta 1e-4 and temperature 0.05 C convergence tolerances, and maximum 5 cycles. Salt is not equilibrated; it is set to the user's explicit zero inventory afterward. A passed baseline spin-up is not a salt steady state or field validation. The update does not install these test controls as a local user PASS; the installed application must perform its own matching checks.

The matching QA and spin-up also passed the actual freeze gate. A subsequent **72-hour R00 transport run using that spun-up initial state completed**, and a repeated request correctly reused its unchanged completed output. This is a deliberately truncated workflow test, not the full 2010-2023 historical run. See `post_spinup_workflow.json`.

Two further **168-hour historical cold-weather windows at 400 cells** also completed: the window surrounding the lowest air temperature and the window surrounding the coldest 72h mean. Minimum internal material temperatures were approximately 0.859 C and 1.572 C. These positive-temperature tests do not remove the guard for colder future scenarios. See `historical_cold_checks.json`.

## Full 21-case scope audit - important for the paper

All 21 hourly climate files pass schema/calendar/source-checksum checks: **3,804,024 records**. But these **eight complete future cases have subzero air temperatures**:

| Case | Minimum air temperature(C) | Subzero hours |
|---|---:|---:|
| R01 | -1.452094 | 4 |
| R13 | -0.216498 | 2 |
| R14 | -3.052374 | 18 |
| R16 | -0.393134 | 2 |
| R17 | -3.285132 | 9 |
| R18 | -3.382208 | 14 |
| R19 | -0.818884 | 2 |
| R20 | -3.511359 | 10 |

The existing reduced model has no ice/freezing physics. Negative air temperature is its conservative scope guard, not evidence that every masonry cell froze. Runtime material temperatures are now checked too. The new compatibility report identifies unsupported full cases **before** beginning a selected batch.

Do not clip temperatures to zero, discard cold hours or silently omit GCMs to get the planned comparison. The original full 21-run paper experiment still requires a documented scope resolution or a freeze-thaw-capable model. No complete 21-case transport ensemble was executed in this repair. Default sulfate is dissolved-equivalent only: no verified crystal mass, hydrate cycles or damage is generated.

## Installation and preservation

Extract the update ZIP. Close the old server/terminal. Double-click **INSTALL_UPDATE_WINDOWS.bat** inside the extracted update folder. Choose the ORIGINAL application folder containing **START_WINDOWS.bat**. The installer verifies payload hashes, backs up replaced code and old controls, and does not replace your project.json, source material, climate or outputs. Restart normally; confirm version 1.1.0. Read/accept the new wet-end assumption only for a documented reduced scenario. Keep 200 cells, save, and run R00 QA followed by water/heat spin-up and a new freeze.

Python dependency versions have not been changed. Native Windows execution was not available here. **Eleven installer preservation, checksum and rollback checks passed.** The update was also applied to a disposable copy of the original application: all 57 tests passed there; the material remained valid and the new review/QA/spin-up/freeze gates were not bypassed. Records are in `installer_checks.json` and `installed_application_smoke.json`. The update does not contain a fabricated QA pass for installation into your project.

## Reproducibility identity

- Uploaded material SHA-256: `b3280e6bfd9fd0baa2db335c18fac10cbb2d34a93e832ebc057efa2dabff0f48`
- Original uploaded project SHA-256: `07077ab6f9d65bbb1b95aa99be09b2d461bc7c31f92b60d7ffb50ae6549953bc`
- R00 climate SHA-256: `ba4a600d70e3912a124f80c52c057c5772d3243540ff060a34e0d23f9ca14c8f`
- Verified revised physics hash: `b5ef41447b2c76a13d438cb1cff71d11c44d38023571b6db98b9d10a1deca6cf`
- Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, pandas 2.2.3; Linux; one BLAS thread.

Code hashes, raw QA, stress-window hourly values, test logs and the original failure reproduction are included under `verification/repair_110`. Numeric source settings were compared field-for-field. Only the explicit new boundary-review flag was added to the test copy. Read `docs/METHODS.md` for equations, thresholds, output units and limitations.
