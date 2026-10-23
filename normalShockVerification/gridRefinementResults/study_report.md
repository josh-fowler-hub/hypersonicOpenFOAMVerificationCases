# Normal-shock grid refinement study

Analyzed 30 grids. Uniform 1D mesh size is h=L/N (not cbrt(V/N)).

Targets (fractions): field cell-average L1/jump ≤ 0.001; shock error/L ≤ 0.0001; scalar GCI/scale ≤ 0.001. These are configurable engineering choices, not universal acceptance standards.

| Cells | h (m) | Max normalized L1 | Shock error/L | Accuracy passes | GCI supported |
|---:|---:|---:|---:|:---:|:---:|
| 100 | 1.000000e-02 | 1.144653e-02 | 2.942695e-03 | false | false |
| 121 | 8.264463e-03 | 8.238487e-03 | 9.288060e-04 | false | false |
| 146 | 6.849315e-03 | 6.723630e-03 | 8.225136e-04 | false | false |
| 176 | 5.681818e-03 | 6.843368e-03 | 5.099326e-04 | false | false |
| 212 | 4.716981e-03 | 4.715753e-03 | 5.176872e-04 | false | false |
| 255 | 3.921569e-03 | 4.379504e-03 | 1.139749e-03 | false | false |
| 307 | 3.257329e-03 | 2.695178e-03 | 1.825514e-04 | false | false |
| 369 | 2.710027e-03 | 3.690830e-03 | 1.422923e-03 | false | false |
| 443 | 2.257336e-03 | 3.293500e-03 | 1.562750e-03 | false | false |
| 532 | 1.879699e-03 | 1.919445e-03 | 2.192511e-04 | false | false |
| 639 | 1.564945e-03 | 2.079790e-03 | 8.084991e-04 | false | false |
| 767 | 1.303781e-03 | 1.114502e-03 | 1.657269e-04 | false | false |
| 921 | 1.085776e-03 | 1.245567e-03 | 6.897974e-05 | false | false |
| 1106 | 9.041591e-04 | 8.600759e-04 | 8.445458e-05 | true | false |
| 1328 | 7.530120e-04 | 1.045171e-03 | 4.846321e-04 | false | false |
| 1594 | 6.273526e-04 | 7.376461e-04 | 2.808112e-04 | false | false |
| 1913 | 5.227392e-04 | 6.754228e-04 | 3.549928e-04 | false | false |
| 2000 | 5.000000e-04 | 4.820742e-04 | 1.353469e-04 | false | false |
| 2296 | 4.355401e-04 | 3.853823e-04 | 1.255483e-05 | true | false |
| 2401 | 4.164931e-04 | 3.914773e-04 | 5.920065e-07 | true | false |
| 2756 | 3.628447e-04 | 2.499218e-04 | 1.851065e-05 | true | false |
| 2882 | 3.469813e-04 | 2.530239e-04 | 5.288324e-05 | true | false |
| 3308 | 3.022975e-04 | 3.630509e-04 | 1.771351e-04 | false | false |
| 3459 | 2.891009e-04 | 2.909008e-04 | 1.137783e-04 | false | false |
| 4151 | 2.409058e-04 | 1.810144e-04 | 1.418008e-05 | true | false |
| 4981 | 2.007629e-04 | 1.161910e-04 | 2.497596e-06 | true | false |
| 5978 | 1.672800e-04 | 1.943105e-04 | 1.004919e-04 | false | false |
| 7174 | 1.393922e-04 | 1.343048e-04 | 5.595407e-05 | true | false |
| 8609 | 1.161575e-04 | 1.183345e-04 | 4.066540e-05 | true | false |
| 10331 | 9.679605e-05 | 5.955828e-05 | 2.837767e-06 | true | false |

## Mesh choice

No tested grid meets all accuracy and supported-GCI criteria across itself and all finer levels. No GCI-qualified mesh recommendation is justified.

Accuracy-only candidate: 7174 cells, h=0.00013939224979089704 m. This uses known analytical errors and is not a spatial-uncertainty certification.

## Finest-grid scalar uncertainty

All scalar values are signed analytical defects; exact=0. Shock uses domain length; fields use their analytical jump. Physical GCI units are meters for shock and the field's units for mean/plateau defects.

| Quantity | p | GCI/scale | Fourth-grid ratio | Supported | Exact covered |
|---|---:|---:|---:|:---:|:---:|
| T_mean_defect | NaN | NaN | NaN | false | false |
| T_plateau_defect | 1.0142 | 3.982971e-05 | 0.9961 | true | true |
| p_mean_defect | NaN | NaN | NaN | false | false |
| p_plateau_defect | 1.0402 | 1.805161e-05 | 0.9896 | true | true |
| rho_mean_defect | NaN | NaN | NaN | false | false |
| rho_plateau_defect | 1.0262 | 4.833380e-05 | 0.9931 | true | true |
| shock_offset | NaN | NaN | NaN | false | false |
| u_mean_defect | NaN | NaN | NaN | false | false |
| u_plateau_defect | 1.0580 | 1.028024e-05 | 0.9856 | true | true |

## Method and limits

Errors are recomputed from profiles, avoiding the old log-parser field-overwrite bug. Both point-sampled and exact cell-average step references are tabulated. Field norms are normalized by the jump, not mixed physical units. Cell averaging of primitive T and u is a diagnostic comparison, not a proof of conservative variable accuracy.

The exact step location and time come from each verification log. Signed shock and domain-mean defects subtract the moving exact solution at that time. This removes known exact motion but does not establish common-time numerical equivalence. Logs round shock locations to 1e-9 m; smaller errors are not resolved by these inputs.

Triplets use actual unequal spacing ratios. We solve (phi3-phi2)/(phi2-phi1)=(h3^p-h2^p)/(h2^p-h1^p), extrapolate phi0=phi1+(phi1-phi2)/(r21^p-1), and estimate U1=1.25*abs(phi1-phi2)/(r21^p-1). Intervals are phi1 ± U1. GCI is an uncertainty indicator, not a statistical confidence interval or total uncertainty.

Oscillating, stagnant, or non-positive-order sequences receive no extrapolation. Local fitted orders above 2 are retained for inspection but excluded from mesh selection. A fitted triplet alone is not evidence of asymptotic convergence: an independent fourth-grid difference must agree with the fitted model within a ratio of 0.8–1.25. This threshold and the order ceiling are conservative study heuristics. Exact-error coverage is reported separately and required for selection.

Adjacent and every-other-grid triplets are reported. The latter increase ratio separation to approximately 1.44 for a 20% sequence. Selection uses adjacent triplets consistently; every-other triplets are sensitivity diagnostics. No pointwise GCI is applied at the discontinuity. Linf and shock-position errors may oscillate with shock/cell alignment; a small value on one grid is insufficient.

Adaptive timesteps vary with mesh. These data contain combined space/time effects; repeat at reduced timesteps or Courant numbers on at least the candidate and finer grids. Compare changes in each quantity with its tolerance and spatial estimate before calling the result a spatial GCI. Model, input, and temporal uncertainty are not quantified or combined here. Mean defects are not conservation-balance residuals; those require boundary flux histories.

For a follow-up study use systematic grids with ratios around 1.4–2, keep physical inputs fixed, and compare a common saved physical time. Further cell-count predictions from oscillating sequences are deliberately omitted; no reliable required N follows from them.

Method reference: [NASA NPARC spatial convergence tutorial](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).

Outputs: [mesh summary](grid_refinement_summary.csv), [field errors](field_errors.csv), [GCI triplets](gci_triplets.csv), [mesh assessment](mesh_selection.csv).
