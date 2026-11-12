# Paper audit (53/54 pass)

| Where | Claim | Result | Value |
|---|---|---|---|
| IV setup | alpha_s = 4*pi*sigma/lambda^2 = 8.2e3 (sigma = -11.235 dBsm) | PASS | 8237, code 8237 |
| IV setup | p1 = 0.032 (S = 20 m) | PASS | 0.0320 |
| IV setup | p1 = 0.278 (S = 8 m) | PASS | 0.2779 |
| Contrib. 2 | p1 decays as S^-3 for wide S (asymptote within 5% at S = 40 m) | PASS | asym/exact at S=40: 1.000 |
| Fig. 2 valid. | all interferers change CR by at most 1.1% | PASS | 1.05% |
| Fig. 2 valid. | all interferers change SR by at most 0.4% | PASS | 0.37% |
| Fig. 2 valid. | Theorem 1/2 bounds below MC for all P | PASS |  |
| Fig. 2a | M=2: bound 29.24 vs MC 29.45 at 30 dBm | PASS | 29.239 / 29.453 |
| Fig. 2a | M=5: bound 43.6 vs MC 51.4, within about 15% | PASS | 43.61 / 51.43 (15.2%) |
| Fig. 2a | CR EDMA/TDMA = 1.85-1.9 x (M=2) | PASS | 1.845-1.902 |
| Fig. 2a | CR EDMA/TDMA = 2.74-2.85 x (M=5) | PASS | 2.741-2.853 |
| Fig. 2a | sum_m Pr{A_m^c} = 1.94 (M=2) | PASS | 1.936 |
| Fig. 2a | sum_m Pr{A_m^c} = 2.93 (M=5) | PASS | 2.934 |
| Fig. 2b | M=2: SR bound 2.22 vs MC 2.28 | PASS | 2.222 / 2.285 |
| Fig. 2b | M=5: SR bound 3.83 vs MC 5.61 | PASS | 3.828 / 5.614 |
| Fig. 2b | SR EDMA/TDMA = 2.7-7.1 x (M=2) | PASS | 2.720-7.072 |
| Fig. 2b | SR EDMA/TDMA = 4.3-5.7 x (M=5) | PASS | 4.261-5.724 |
| Fig. 2b | SR gain 7.1x at 0 dBm (M=2) | PASS | 7.072 |
| Fig. 2b | SR gain stays above 2.7x at 30 dBm | PASS | M=2 2.720, M=5 4.261 |
| Fig. 3a | CR peaks at 51.4 for M=5 (D_L=40) | PASS | M*=5, 51.43 |
| Fig. 3a | CR peaks at 29.3 for M=3 (D_L=20) | PASS | M*=3, 29.26 |
| Fig. 3a | CR peaks at 18.0 for M=2 (D_L=10) | PASS | M*=2, 17.95 |
| Fig. 3a | CR non-monotone in M (D_L=10) | PASS | signs of increments: +---+++ |
| Fig. 3a | CR non-monotone in M (D_L=20) | PASS | signs of increments: ++----- |
| Fig. 3a | CR non-monotone in M (D_L=40) | PASS | signs of increments: ++++--- |
| Fig. 3a | EDMA leads TDMA (CR) for all M (D_L=20) | PASS | min ratio 1.046 |
| Fig. 3a | EDMA leads TDMA (CR) for all M (D_L=40) | PASS | min ratio 1.073 |
| Fig. 3a | even where p1 > 1/2 (D_L=20, M>=4) | PASS | p1(S=5) = 0.505 |
| Fig. 3a | D_L=10, M>=3: TDMA is better | PASS | ratios 0.91 0.83 0.81 0.81 0.83 0.85 |
| Fig. 3a | M=3 at D_L=10 means S narrower than 3.3 m | PASS | 3.333 |
| Fig. 3a | isolate only 0.79-0.61 segments on average (D_L=10, M=3..8) | PASS | 0.793 0.685 0.633 0.611 0.606 0.610 |
| Fig. 3a | D_L=10, M=2: p1 = 0.51 | PASS | 0.5052 |
| Fig. 3a | D_L=10, M=2: EDMA 18.0 vs TDMA 17.0 | PASS | 17.95 / 17.02 |
| Sec. IV | p1 < 1/2 is necessary for M=2 (asymptotically: 2(1-p1) < 1 at p1=0.505) | PASS | 2(1-p1) = 0.990 |
| Fig. 3b | SR increases with M (D_L=40) | PASS | increments: +1.496 +1.455 +1.129 +0.745 +0.426 +0.208 +0.076 |
| Fig. 3b | SR saturates instead of decreasing (D_L=10; dip after peak < 5%) | PASS | M=2..8: 2.02-2.24, max dip from peak 0.0% |
| Fig. 3b | SR saturates instead of decreasing (D_L=20; dip after peak < 5%) | PASS | M=2..8: 2.57-3.49, max dip from peak 2.7% |
| Fig. 3b | M=1: EDMA 0.79 vs TDMA 0.52 (D_L=40) | PASS | 0.789 / 0.520 |
| Fig. 3b | D_L=40: SR 2.28 -> 6.32 for M=2..8 | PASS | 2.285 -> 6.324 |
| Fig. 3b | D_L=40: TDMA 0.84-1.51 | PASS | 0.840 -> 1.513 |
| Fig. 3b | D_L=40: p1 <= 0.51 for M<=8 | PASS | 0.505 |
| Fig. 3b | SR saturates near 3.5 for M>=4 (D_L=20) | PASS | 3.48 3.49 3.45 3.41 3.39 |
| Fig. 3b | SR near 2.0-2.2 (D_L=10, M>=2) | PASS | 2.022-2.239 |
| Fig. 3b | D_L=10, M=8: 2.24 vs 2.05 | PASS | 2.239 / 2.047 |
| Fig. 3b | ULA below EDMA and TDMA in all cases (vs M, CR and SR) | PASS | max ULA/TDMA CR 0.94, SR 0.71 |
| Info (no claim) | ULA below TDMA vs P (CR and SR) | **FAIL** | max ULA/TDMA SR M=2: 1.09 at 0 dBm |
| Conclusion | best M for communication grows with D_L (2, 3, 5) | PASS |  |
| Limit M=1 | bound exact at M=1 (D_L=10), CR and SR within 0.1% | PASS | CR 15.8351/15.8348, SR 1.4300/1.4299 |
| Limit M=1 | bound exact at M=1 (D_L=20), CR and SR within 0.1% | PASS | CR 15.1407/15.1405, SR 1.1625/1.1624 |
| Limit M=1 | bound exact at M=1 (D_L=40), CR and SR within 0.1% | PASS | CR 14.0164/14.0168, SR 0.7888/0.7892 |
| Limit phi large | phi=0.1: MC and bound approach M*Bc (M=2) | PASS | bound 30.23 <= MC 30.24 <= M*Bc 30.28 |
| Limit phi=0 | no blockage: EDMA CR below TDMA (M=2) | PASS | EDMA 8.57 vs TDMA 15.48; bound 0.00 |
| Limit phi large | phi=0.1: MC and bound approach M*Bc (M=5) | PASS | bound 77.08 <= MC 77.40 <= M*Bc 79.86 |
| Limit phi=0 | no blockage: EDMA CR below TDMA (M=5) | PASS | EDMA 11.27 vs TDMA 18.03; bound 0.00 |
