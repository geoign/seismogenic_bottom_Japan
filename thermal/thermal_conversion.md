# Thermal conversion of the base of brittle seismicity

Made by `scripts/10_thermal_conversion.py`. D_b10 is treated as the depth of an isotherm whose temperature T_c is calibrated against onshore heat-flow observations; a conductive geotherm through the surface and (D_b10, T_c) then gives heat flow, gradients and isotherm depths. Independent work based on the JMA unified hypocenter catalogue; not a product of JMA or of the providers of the heat-flow and strain-rate data.

## 1. Calibration

* Heat-flow observations: 612 onshore entries of the IHFC Global Heat Flow Database (Release 2021) in 122-148E, 24-46N; 249 of them fall on 161 nodes with a usable D_b10 (quality_flag 0, onshore strain-rate coverage, node median 20-300 mW/m2).
* Spearman correlation between observed heat flow and D_b10 at these nodes: -0.54.
* Geotherm parameters: k = 2.5 W/m/K, A0 = 2.0 microW/m3, hr = 10.0 km, T_s = 16 - 0.9 (lat - 35) deg C.
* T_c model: 1/T_c = 1/T_ref + (R/Q) [n ln(z/15 km) - ln(eps/eps_ref)], n = 3, Q = 356 kJ/mol, eps_ref = 89 nanostrain/yr (median second invariant over the domain).
* **Calibrated T_ref = 427 deg C (95 % bootstrap interval 395-476 deg C)**: the temperature at D_b10 when D_b10 = 15 km and the strain rate is at its median.
* Residual ln(q_obs / q_pred): median -0.000, robust s.d. 0.44 (a factor of 1.55); Spearman correlation of predicted and observed heat flow +0.57. The scatter is that of single heat-flow sites against a 10-40 km average of seismicity, so it is an upper bound on the error of the grid.

The same calibration applied to other depth statistics of the same product (constant strain-rate and depth terms as above):

| Depth statistic | Median depth over the domain (km) | Calibrated T_ref (deg C) | 95 % interval |
| --- | --- | --- | --- |
| D_b10 (10 % of the plateau rate) | 14.8 | 427 | 395-479 |
| classic D90 (family-weighted) | 13.5 | 397 | 359-435 |
| z_b (centre of the transition) | 12.4 | 360 | 335-392 |
| D_b90 (top of the transition) | 10.1 | 277 | 262-303 |

## 2. What the heat-flow data say about the depth and strain-rate terms

Fitting the two coefficients freely instead of fixing them from the flow law:

| Coefficient | Flow law (wet anorthite) | Flow law (wet quartz, n = 4, Q = 135 kJ/mol) | Free fit | 95 % interval |
| --- | --- | --- | --- | --- |
| depth term n R/Q (1/K) | 7.01e-05 | 2.46e-04 | 4.68e-04 | 2.41e-05 to 7.54e-04 |
| strain-rate term R/Q (1/K) | 2.34e-05 | 6.16e-05 | 2.63e-04 | 7.43e-05 to 4.08e-04 |

Mean absolute residual: 0.349 (flow-law coefficients, this product), 0.322 (free fit). The depth term is not constrained by the data. The strain-rate term has the sign the flow law predicts but is several times larger; since hot, weak crust also strains faster, and the geodetic field (2019-2023) contains elastic strain from plate coupling and post-seismic deformation of the 2011 Tohoku-oki earthquake, the free-fit value is read as a statistical association with heat flow, not as rheology, and is not used.

## 3. Sensitivity to the geotherm parameters

T_ref is re-calibrated for each parameter set, so the predicted heat flow is anchored to the observations; the temperature assigned to D_b10, and with it the depth of the 350 deg C isotherm, is what moves.

| k (W/m/K) | A0 (microW/m3) | T_ref (deg C) | median heat flow (mW/m2) | median depth of 350 deg C (km) | median depth of 580 deg C (km) |
| --- | --- | --- | --- | --- | --- |
| 2.0 | 2.0 | 530 | 79 | 9.3 | 16.4 |
| 2.5 | 1.0 | 455 | 79 | 11.2 | 19.2 |
| 2.5 | 2.0 | 427 | 79 | 11.8 | 20.9 |
| 2.5 | 3.0 | 403 | 80 | 12.5 | 22.9 |
| 3.0 | 2.0 | 358 | 79 | 14.4 | 25.5 |

## 4. The grids

4070 nodes (0.1 degree; quality_flag 0 and onshore strain-rate coverage). Quantiles over these nodes:

| Quantity | 5 % | 25 % | 50 % | 75 % | 95 % |
| --- | --- | --- | --- | --- | --- |
| D_b10, smoothed input (km) | 10.4 | 12.7 | 14.8 | 17.2 | 24.1 |
| T at D_b10 (deg C) | 409 | 420 | 426 | 433 | 444 |
| surface heat flow (mW/m2) | 53.6 | 69.6 | 79.2 | 91.0 | 110.4 |
| surface gradient (deg C/km) | 21.4 | 27.8 | 31.7 | 36.4 | 44.2 |
| mean gradient to D_b10 (deg C/km) | 16.4 | 23.6 | 27.9 | 33.0 | 41.1 |
| depth of 350 deg C (km) | 8.1 | 10.0 | 11.8 | 13.9 | 19.9 |
| depth of 580 deg C (km) | 14.0 | 17.6 | 20.9 | 24.9 | 36.4 |
| depth of 350 deg C / D_b10 | 0.77 | 0.79 | 0.80 | 0.81 | 0.83 |

With T_ref at the ends of its 95 % interval the median heat flow is 74 and 88 mW/m2 (bands 7-8). The depth of 350 deg C depends more on the conductivity than on the calibration error: its median is 10.5 km with T_ref at the upper end of the interval (band 9) and 14.4 km with k = 3.0 W/m/K, for which T_ref recalibrates to 358 deg C (band 10; the surface gradient is then heat flow / 3.0).

## 5. Independent check against a Curie-point depth model

Li et al. (2017) global Curie-point depth (magnetite, about 580 deg C) at the same 4070 nodes: median 25.3 km (25-75 %: 22.5-28.1 km). The 580 deg C depth predicted here: median 20.9 km (17.6-24.9 km). Spearman correlation between the two maps: +0.11; between the Curie depth and D_b10: +0.11. The Curie model carries no node-scale information here, so only the regional medians can be compared. With k = 3.0 W/m/K (T_ref = 358 deg C) the predicted median is 25.5 km. The heat-flow calibration with k = 2.5 and the Curie median therefore bracket the temperature at D_b10 between about 358 and 427 deg C, and the 350 deg C isotherm between 0.80 and 0.97 times D_b10. The heat-flow grid is the same in both cases.

## 6. Limits

* Steady-state, one-dimensional conduction. Not valid next to hydrothermal systems and young intrusions (where the shallow D_b10 reflects advected heat), nor in the fore-arc cooled by the slab. The 0.025-degree product should not be converted this way.
* T_c is tied to observations only through one regional constant. Lithology (mafic lower crust, D_b10 deeper than about 25 km) and pore pressure are not modelled.
* Onshore heat-flow sites over-represent geothermal areas; node medians and the 20-300 mW/m2 window limit, but do not remove, the resulting upward bias of T_ref.
* Depths are JMA routine depths (1-D velocity model); a systematic depth error of 1-3 km maps into T_ref.
* k is constant; its decrease with temperature would steepen the geotherm at depth and make the 580 deg C depth shallower than given here.
