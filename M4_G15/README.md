# M4 (Final Submission) — Performance Analysis of Multi-Segment Pinching-Antenna ISAC Under LoS Blockage

**Category:** Research Project
**Student:** Khushi Shah (AU2440045), single-member group
**Mentor:** Prof. Dhaval Patel
**Milestones:** M1 Proposal & Scoping → M2 Technical Artifact → M3 Results & Evaluation → **M4 Final Submission**

## Project summary

This project analyzes how line-of-sight (LoS) blockage affects a multi-segment pinching-antenna system (PASS) that performs integrated sensing and communication (ISAC) under environment-division multiple access (EDMA): each waveguide segment activates one pinching antenna at the midpoint of its user and target, all segments operate concurrently, and blockage of the cross-segment links jointly governs communication interference, sensing clutter, and sensing leakage.

Across the four milestones the project:
- derived a closed-form mean cross-segment LoS probability $p_1(\phi)$ and closed-form lower bounds on the ergodic sum communication rate (CR) and sensing rate (SR), for any number of segments $M$ (M2);
- verified both against a Monte Carlo simulator (M2);
- ran the full numerical evaluation — CR/SR vs. transmit power $P$ and vs. $M$, for three waveguide lengths — under a physically grounded automated-guided-vehicle (AGV) target radar cross-section (RCS) and a conventional uniform linear array (ULA) benchmark, and validated the high-SNR design rule $p_1(\phi)<1/2$ for when EDMA outperforms TDMA (M3);
- consolidated the full project, including the final SOTA positioning, limitations, and future work (**M4, this submission**).

**Headline results:** EDMA nearly triples the sum CR of TDMA and raises the sum SR up to sevenfold; a conventional ULA trails both throughout (at most $94\%$ of TDMA's CR, $71\%$ of its SR). The design rule correctly flags EDMA's one failure regime (short waveguides, many segments) but is shown to be necessary, not sufficient, away from that regime.

## Repository structure

```
M4_G15/
├── README.md              <- this file
├── Report/
│   ├── M4_G15_Report.pdf  <- final consolidated report
│   └── src/                <- LaTeX source
├── Video/
│   ├── M4_G15_video_link.md  <- link to the 20-minute walkthrough (fill in after recording)
│   └── video_script.md       <- section-by-section speaking script
├── Code/
│   ├── src/                <- simulator, plotting, and audit scripts
│   └── README.md
├── Data/
│   └── README.md           <- no external dataset; see this file for how all numbers are generated
└── Results/
    ├── figures/             <- the four report figures (CR/SR vs. P and vs. M)
    └── audit_report.md      <- independent numerical verification of every claim (53/54 pass)
```

## Installation and execution

```bash
pip install numpy scipy matplotlib
cd Code/src
python plot_v3.py       # regenerates Results/figures/
python audit_paper.py   # independently re-checks every numerical claim in the report
```

No external dataset is used (see `Data/README.md`); every number and figure is generated from the fixed physical parameters listed in `Code/README.md` and the model in `Code/src/edma_sim.py`.

## Results at a glance

| Scheme | CR, $M=2$ | SR, $M=2$ | CR, $M=5$ | SR, $M=5$ |
|---|---|---|---|---|
| EDMA, PA at midpoint (MC) | 29.45 | 2.28 | 51.43 | 5.61 |
| Lower bound (Theorems 1–2) | 29.24 | 2.22 | 43.61 | 3.83 |
| TDMA | 15.48 | 0.84 | 18.03 | 1.32 |
| Conventional ULA | $\le$ TDMA throughout (at most $94\%$ of TDMA's CR, $71\%$ of its SR) | | | |

($P=30$~dBm, $D_L=40$~m; bits/s/Hz. See `Results/audit_report.md` for every verified number and `Report/M4_G15_Report.pdf` Sections 7–9 for the full discussion, including where the bounds loosen and where EDMA's design rule is only a necessary condition.)

## Individual contribution

All work in this single-member project — literature review, modeling, derivations, simulator, evaluation, critical analysis, and report, across all four milestones — was carried out by Khushi Shah under the mentorship of Prof. Dhaval Patel.
