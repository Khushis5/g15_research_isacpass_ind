# Code

Python 3.11 with `numpy`, `scipy`, and `matplotlib`.

```bash
cd src
python plot_v3.py
python audit_paper.py
```

`plot_v3.py` writes the four evaluation figures (CR/SR vs. $P$, CR/SR vs. $M$) to `../../Results/figures/`. `audit_paper.py` is the independent numerical check behind every claim in the report; its output is saved in `../../Results/audit_report.md` (53/54 checks pass — see that file for the one flagged item, a low-power edge case the report does not rely on).

## Files

| File | Purpose |
|---|---|
| `src/edma_sim.py` | System model and Monte Carlo simulator: EDMA (pinching antenna at the midpoint), TDMA, and the conventional ULA benchmark. Holds every closed form ($I(a)$, $p_1(\phi)$, and the Theorem 1/2 lower bounds) and the AGV radar-cross-section parameterization ($\alpha_s = 4\pi\bar\sigma/\lambda^2$). |
| `src/plot_v3.py` | Regenerates the four report figures. |
| `src/audit_paper.py` | Independently re-derives or re-checks every numerical claim in the report (bound values, ratios, limiting cases). |
| `src/sim_cache_final.npz` (generated) | Cached Monte Carlo results, keyed by a hash of every model parameter, so a parameter change can never reuse stale data. Delete it to force a full recompute. |

## Parameters

This project uses no external dataset; see `../Data/README.md`. All data is generated from:

| Parameter | Value |
|---|---|
| Carrier frequency $f_c$ | 28 GHz |
| Waveguide height $d$ | 3 m |
| Waveguide length $D_L$ | 10, 20, 40 m |
| Number of segments $M$ | 1–8 |
| Noise power $\sigma^2 = \sigma_s^2$ | −90 dBm |
| Probe length $L$ | 5 |
| Target RCS $\bar\sigma$ (AGV) | −11.235 dBsm ($\alpha_s = 8237$) |
| Blockage parameter $\phi$ | 0.02 m$^{-2}$ |
| Transmit power $P$ | 0–30 dBm |

## Notes

- Monte Carlo: $10^5$ user/target drops (seed 2718). For each drop, the blockage expectation is evaluated exactly over every pattern of the cross links (and, for the ULA, every pattern of the $M$ shared element blockages).
- The closed-form bounds use the adjacent-only interference/clutter/leakage model; `edma_mc(..., nearest=False)` includes every segment instead, used for the ablation check in the report (Section 7/8).
