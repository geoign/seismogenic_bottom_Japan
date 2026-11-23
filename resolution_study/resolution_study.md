# Resolution study of the seismogenic-bottom map

## 1. Data-limited resolution (smallest radius with ESS >= 30 families and n >= 60 events)

Cumulative fraction of 0.025-degree nodes whose smallest admissible radius is <= r:

| r (km) | land nodes of Honshu, Kyushu, Shikoku (130-142E, 31-41.5N) | all nodes that reach 40 km |
| --- | --- | --- |
| 2 | 14.8 % | 5.0 % |
| 3 | 26.6 % | 9.7 % |
| 4 | 39.0 % | 15.1 % |
| 5 | 50.3 % | 20.5 % |
| 7 | 68.3 % | 30.8 % |
| 10 | 85.0 % | 44.0 % |
| 15 | 96.1 % | 59.9 % |
| 20 | 98.9 % | 71.2 % |
| 30 | 99.8 % | 87.6 % |
| 40 | 99.9 % | 100.0 % |

## 2. Reproducibility between two independent halves of the catalogue (fixed radius, 0.1-degree nodes)

Families (declustered sequences) were split by parity of their id into two halves, and the fit was repeated on each half with ESS >= 15. The noise of the full-data estimate is RMS(half A - half B) / 2. The fine-scale texture is the deviation of D_b10 from its 0.5-degree running mean; its reproducibility is the correlation of that texture between the two halves, and its signal is its RMS in the full-data map after removing the noise.

| radius (km) | nodes | RMS(A-B) km | noise of D_b10 (RMS) km | noise (robust, MAD) km | median s.e.(z_b) from Hessian km | corr(A, B) | RMS texture of D_b10 km | corr of texture (A, B) | texture signal / noise |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 547 | 1.01 | 0.50 | 0.23 | 0.31 | 0.991 | 2.09 | 0.853 | 4.0 |
| 3 | 1058 | 0.85 | 0.43 | 0.26 | 0.31 | 0.994 | 2.51 | 0.904 | 5.8 |
| 5 | 2303 | 0.95 | 0.48 | 0.25 | 0.31 | 0.994 | 2.59 | 0.907 | 5.4 |
| 10 | 5116 | 1.08 | 0.54 | 0.25 | 0.30 | 0.992 | 2.51 | 0.900 | 4.5 |
| 20 | 8443 | 1.09 | 0.55 | 0.20 | 0.27 | 0.993 | 2.30 | 0.905 | 4.1 |

## 3. The 0.025-degree product

Valid nodes: 195040 (0.025 deg) vs 12284 (0.1 deg).

Search radius actually used by the valid 0.025-degree nodes:

| radius (km) | nodes | share | median s.e.(z_b) km | median ESS |
| --- | --- | --- | --- | --- |
| 3 | 19686 | 10.1 % | 0.34 | 122 |
| 5 | 21436 | 11.0 % | 0.44 | 88 |
| 7 | 20586 | 10.6 % | 0.46 | 79 |
| 10 | 26019 | 13.3 % | 0.51 | 84 |
| 15 | 30870 | 15.8 % | 0.56 | 93 |
| 20 | 21366 | 11.0 % | 0.66 | 82 |
| 30 | 31330 | 16.1 % | 0.58 | 101 |
| 40 | 23747 | 12.2 % | 0.53 | 89 |

Of the 12181 coincident nodes valid in both products, 8331 use a radius >= 10 km in the fine product and are therefore identical to the coarse product. Where the fine product uses a smaller radius, D_b10(0.025 deg) - D_b10(0.1 deg) is:

| fine radius (km) | nodes | median difference km | robust s.d. km | abs. difference > 3 km |
| --- | --- | --- | --- | --- |
| 3 | 1224 | -0.61 | 1.16 | 12.4 % |
| 5 | 1354 | -0.35 | 1.01 | 7.8 % |
| 7 | 1272 | -0.15 | 0.68 | 5.1 % |

## 4. Split-half test of the 0.025-degree product: kyushu (129.5-132.5E, 31-34N)

13417 nodes valid in both halves. RMS(A-B) = 1.10 km, robust = 0.62 km, so the noise of the full-data D_b10 is about 0.31 km (robust). corr(A, B) = 0.939; correlation of the anomalies from a 0.275-degree (11-node) running mean = 0.751 (how much of the fine texture is reproduced by independent data).

## 4. Split-half test of the 0.025-degree product: tohoku_kanto (138-142E, 35-40N)

27452 nodes valid in both halves. RMS(A-B) = 1.60 km, robust = 0.85 km, so the noise of the full-data D_b10 is about 0.43 km (robust). corr(A, B) = 0.987; correlation of the anomalies from a 0.275-degree (11-node) running mean = 0.820 (how much of the fine texture is reproduced by independent data).

