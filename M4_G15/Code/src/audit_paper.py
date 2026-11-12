"""Audit: recompute every numeric claim of main_ver2_4pg_v3.tex from the current model and
report PASS / FAIL.  Rerun after any change to the model or the text.
Also runs the limiting-case (backward-compatibility) checks:
  M = 1          -> bounds are exact (single-segment ISAC-PASS)
  phi -> large   -> bounds and MC approach the isolated limit M*B
  phi = 0        -> blockage-free model (all cross links LoS)
Output: printed table + audit_report.md next to this file."""
import math, os
import numpy as np
import edma_sim as e
import plot_phi as pp

HERE = os.path.dirname(os.path.abspath(__file__))
rows = []


def check(where, claim, ok, value):
    rows.append((where, claim, "PASS" if ok else "FAIL", value))


def r(x, n):                       # paper rounding
    return round(float(x), n)


def iso_segments(D_L, M, n=100_000, seed=1):
    """sum_m Pr{all adjacent cross links of segment m blocked} (communication)."""
    xU, xT, t = e._geometry(D_L, M, n, seed)
    tot = 0.0
    for m in range(M):
        nb = e._neighbours(m, M, True)
        pc = np.exp(-e.PHI * ((t[:, nb] - xU[:, [m]]) ** 2 + e.D2))
        tot += np.prod(1 - pc, axis=1).mean()
    return tot


p1 = lambda S: e.p_k(e.PHI, S, 1)
pw = {M: e.power_sweep(M) for M in (2, 5)}
ms = {DL: e.m_sweep(DL) for DL in (10.0, 20.0, 40.0)}
ula_p = {M: e.ula_power_sweep(M) for M in (2, 5)}
ula_m = {DL: e.ula_m_sweep(DL) for DL in (10.0, 20.0, 40.0)}
P = e.P_DBM

# ------------------------------------------------------------ Section IV setup
lam = e.C0 / e.FC
a_s = 4 * math.pi * 10 ** (-11.235 / 10) / lam ** 2
check("IV setup", "alpha_s = 4*pi*sigma/lambda^2 = 8.2e3 (sigma = -11.235 dBsm)",
      abs(a_s - 8.2e3) / 8.2e3 < 0.01 and e.ALPHA_S == 8237.0, f"{a_s:.0f}, code {e.ALPHA_S:g}")
check("IV setup", "p1 = 0.032 (S = 20 m)", r(p1(20), 3) == 0.032, f"{p1(20):.4f}")
check("IV setup", "p1 = 0.278 (S = 8 m)", r(p1(8), 3) == 0.278, f"{p1(8):.4f}")
check("Contrib. 2", "p1 decays as S^-3 for wide S (asymptote within 5% at S = 40 m)",
      abs(math.sqrt(math.pi) * math.exp(-e.PHI * e.D2) / (2 * e.PHI ** 1.5 * 40 ** 3) / p1(40) - 1) < 0.05,
      f"asym/exact at S=40: {math.sqrt(math.pi)*math.exp(-e.PHI*e.D2)/(2*e.PHI**1.5*40**3)/p1(40):.3f}")

# ------------------------------------------------------------ validation paragraph
dc = max(np.max(np.abs(pw[M][6] - pw[M][0]) / pw[M][0]) for M in (2, 5)) * 100
ds = max(np.max(np.abs(pw[M][7] - pw[M][1]) / pw[M][1]) for M in (2, 5)) * 100
check("Fig. 2 valid.", "all interferers change CR by at most 1.1%", r(dc, 1) <= 1.1, f"{dc:.2f}%")
check("Fig. 2 valid.", "all interferers change SR by at most 0.4%", r(ds, 1) <= 0.4, f"{ds:.2f}%")
lbok = all((pw[M][4] <= pw[M][0] + 1e-9).all() and (pw[M][5] <= pw[M][1] + 1e-9).all() for M in (2, 5))
check("Fig. 2 valid.", "Theorem 1/2 bounds below MC for all P", lbok, "")

# ------------------------------------------------------------ Fig. 2a
check("Fig. 2a", "M=2: bound 29.24 vs MC 29.45 at 30 dBm",
      r(pw[2][4][-1], 2) == 29.24 and r(pw[2][0][-1], 2) == 29.45, f"{pw[2][4][-1]:.3f} / {pw[2][0][-1]:.3f}")
gap5 = 1 - pw[5][4][-1] / pw[5][0][-1]
check("Fig. 2a", "M=5: bound 43.6 vs MC 51.4, within about 15%",
      r(pw[5][4][-1], 1) == 43.6 and r(pw[5][0][-1], 1) == 51.4 and gap5 <= 0.155,
      f"{pw[5][4][-1]:.2f} / {pw[5][0][-1]:.2f} ({gap5*100:.1f}%)")
for M, lo, hi in ((2, 1.85, 1.90), (5, 2.74, 2.85)):
    rat = pw[M][0] / pw[M][2]
    check("Fig. 2a", f"CR EDMA/TDMA = {lo}-{hi} x (M={M})",
          r(rat.min(), 2) == lo and r(rat.max(), 2) == hi, f"{rat.min():.3f}-{rat.max():.3f}")
for M, D_L, val in ((2, 40.0, 1.94), (5, 40.0, 2.93)):
    g = iso_segments(D_L, M)
    check("Fig. 2a", f"sum_m Pr{{A_m^c}} = {val} (M={M})", r(g, 2) == val, f"{g:.3f}")

# ------------------------------------------------------------ Fig. 2b
check("Fig. 2b", "M=2: SR bound 2.22 vs MC 2.28",
      r(pw[2][5][-1], 2) == 2.22 and r(pw[2][1][-1], 2) == 2.28, f"{pw[2][5][-1]:.3f} / {pw[2][1][-1]:.3f}")
check("Fig. 2b", "M=5: SR bound 3.83 vs MC 5.61",
      r(pw[5][5][-1], 2) == 3.83 and r(pw[5][1][-1], 2) == 5.61, f"{pw[5][5][-1]:.3f} / {pw[5][1][-1]:.3f}")
for M, lo, hi in ((2, 2.7, 7.1), (5, 4.3, 5.7)):
    rat = pw[M][1] / pw[M][3]
    check("Fig. 2b", f"SR EDMA/TDMA = {lo}-{hi} x (M={M})",
          r(rat.min(), 1) == lo and r(rat.max(), 1) == hi, f"{rat.min():.3f}-{rat.max():.3f}")
rat2 = pw[2][1] / pw[2][3]
check("Fig. 2b", "SR gain 7.1x at 0 dBm (M=2)", r(rat2[0], 1) == 7.1, f"{rat2[0]:.3f}")
check("Fig. 2b", "SR gain stays above 2.7x at 30 dBm",
      min(pw[M][1][-1] / pw[M][3][-1] for M in (2, 5)) >= 2.7 - 0.05,
      f"M=2 {rat2[-1]:.3f}, M=5 {pw[5][1][-1]/pw[5][3][-1]:.3f}")

# ------------------------------------------------------------ Fig. 3a
Ms = np.arange(1, 9)
for DL, Mpk, val in ((40.0, 5, 51.4), (20.0, 3, 29.3), (10.0, 2, 18.0)):
    cr = ms[DL][0]
    check("Fig. 3a", f"CR peaks at {val} for M={Mpk} (D_L={DL:g})",
          int(np.argmax(cr)) + 1 == Mpk and r(cr.max(), 1) == val, f"M*={np.argmax(cr)+1}, {cr.max():.2f}")
for DL in (10.0, 20.0, 40.0):
    d = np.diff(ms[DL][0])
    check("Fig. 3a", f"CR non-monotone in M (D_L={DL:g})", (d > 0).any() and (d < 0).any(),
          "signs of increments: " + "".join("+" if x > 0 else "-" for x in d))
for DL in (20.0, 40.0):
    check("Fig. 3a", f"EDMA leads TDMA (CR) for all M (D_L={DL:g})",
          (ms[DL][0] > ms[DL][2]).all(), f"min ratio {np.min(ms[DL][0]/ms[DL][2]):.3f}")
check("Fig. 3a", "even where p1 > 1/2 (D_L=20, M>=4)", p1(5) > 0.5, f"p1(S=5) = {p1(5):.3f}")
lose = ms[10.0][0][2:] < ms[10.0][2][2:]
check("Fig. 3a", "D_L=10, M>=3: TDMA is better", lose.all(), "ratios " + " ".join(f"{x:.2f}" for x in ms[10.0][0][2:] / ms[10.0][2][2:]))
check("Fig. 3a", "M=3 at D_L=10 means S narrower than 3.3 m", r(10 / 3, 1) == 3.3, f"{10/3:.3f}")
gs = [iso_segments(10.0, M) for M in range(3, 9)]
check("Fig. 3a", "isolate only 0.79-0.61 segments on average (D_L=10, M=3..8)",
      r(gs[0], 2) == 0.79 and r(gs[-1], 2) == 0.61 and max(gs) < 1, " ".join(f"{g:.3f}" for g in gs))
check("Fig. 3a", "D_L=10, M=2: p1 = 0.51", r(p1(5), 2) == 0.51, f"{p1(5):.4f}")
check("Fig. 3a", "D_L=10, M=2: EDMA 18.0 vs TDMA 17.0",
      r(ms[10.0][0][1], 1) == 18.0 and r(ms[10.0][2][1], 1) == 17.0, f"{ms[10.0][0][1]:.2f} / {ms[10.0][2][1]:.2f}")
check("Sec. IV", "p1 < 1/2 is necessary for M=2 (asymptotically: 2(1-p1) < 1 at p1=0.505)",
      2 * (1 - p1(5)) < 1, f"2(1-p1) = {2*(1-p1(5)):.3f}")

# ------------------------------------------------------------ Fig. 3b
d = np.diff(ms[40.0][1])
check("Fig. 3b", "SR increases with M (D_L=40)", (d > 0).all(),
      "increments: " + " ".join(f"{x:+.3f}" for x in d))
for DL in (10.0, 20.0):
    sr = ms[DL][1]; drop = 1 - sr[1:].min() / sr[1:].max()
    dip = 1 - sr[np.argmax(sr):].min() / sr.max()
    check("Fig. 3b", f"SR saturates instead of decreasing (D_L={DL:g}; dip after peak < 5%)",
          dip < 0.05,
          f"M=2..8: {sr[1:].min():.2f}-{sr[1:].max():.2f}, max dip from peak {100*(1-sr[np.argmax(sr):].min()/sr.max()):.1f}%")
check("Fig. 3b", "M=1: EDMA 0.79 vs TDMA 0.52 (D_L=40)",
      r(ms[40.0][1][0], 2) == 0.79 and r(ms[40.0][3][0], 2) == 0.52, f"{ms[40.0][1][0]:.3f} / {ms[40.0][3][0]:.3f}")
check("Fig. 3b", "D_L=40: SR 2.28 -> 6.32 for M=2..8",
      r(ms[40.0][1][1], 2) == 2.28 and r(ms[40.0][1][-1], 2) == 6.32, f"{ms[40.0][1][1]:.3f} -> {ms[40.0][1][-1]:.3f}")
check("Fig. 3b", "D_L=40: TDMA 0.84-1.51", r(ms[40.0][3][1], 2) == 0.84 and r(ms[40.0][3][-1], 2) == 1.51,
      f"{ms[40.0][3][1]:.3f} -> {ms[40.0][3][-1]:.3f}")
check("Fig. 3b", "D_L=40: p1 <= 0.51 for M<=8", r(p1(5), 2) <= 0.51, f"{p1(5):.3f}")
check("Fig. 3b", "SR saturates near 3.5 for M>=4 (D_L=20)",
      np.all(np.abs(ms[20.0][1][3:] - 3.5) < 0.15), " ".join(f"{x:.2f}" for x in ms[20.0][1][3:]))
check("Fig. 3b", "SR near 2.0-2.2 (D_L=10, M>=2)",
      ms[10.0][1][1:].min() >= 1.95 and ms[10.0][1][1:].max() <= 2.25, f"{ms[10.0][1][1:].min():.3f}-{ms[10.0][1][1:].max():.3f}")
check("Fig. 3b", "D_L=10, M=8: 2.24 vs 2.05",
      r(ms[10.0][1][-1], 2) == 2.24 and r(ms[10.0][3][-1], 2) == 2.05, f"{ms[10.0][1][-1]:.3f} / {ms[10.0][3][-1]:.3f}")
check("Fig. 3b", "ULA below EDMA and TDMA in all cases (vs M, CR and SR)",
      all((ula_m[DL][k] < ms[DL][k]).all() and (ula_m[DL][k] < ms[DL][k + 2]).all() for DL in ula_m for k in (0, 1)),
      "max ULA/TDMA CR {:.2f}, SR {:.2f}".format(max(np.max(ula_m[DL][0] / ms[DL][2]) for DL in ula_m),
                                                  max(np.max(ula_m[DL][1] / ms[DL][3]) for DL in ula_m)))
check("Info (no claim)", "ULA below TDMA vs P (CR and SR)",
      all((ula_p[M][k] < pw[M][k + 2]).all() for M in (2, 5) for k in (0, 1)),
      "max ULA/TDMA SR M=2: {:.2f} at {} dBm".format(np.max(ula_p[2][1] / pw[2][3]), P[np.argmax(ula_p[2][1] / pw[2][3])]))

# ------------------------------------------------------------ Conclusion
check("Conclusion", "best M for communication grows with D_L (2, 3, 5)",
      [int(np.argmax(ms[DL][0])) + 1 for DL in (10.0, 20.0, 40.0)] == [2, 3, 5], "")

# ------------------------------------------------------------ limiting cases
for DL in (10.0, 20.0, 40.0):
    m = ms[DL]
    check("Limit M=1", f"bound exact at M=1 (D_L={DL:g}), CR and SR within 0.1%",
          abs(m[4][0] / m[0][0] - 1) < 1e-3 and abs(m[5][0] / m[1][0] - 1) < 1e-3,
          f"CR {m[0][0]:.4f}/{m[4][0]:.4f}, SR {m[1][0]:.4f}/{m[5][0]:.4f}")
for M in (2, 5):
    d = pp.phi_data(M); S = 40 / M; rho = 1 / e.SIG2
    Bc = (e.I(e.D2 + rho * e.ETA, S) - e.I(e.D2, S)).real
    check("Limit phi large", f"phi=0.1: MC and bound approach M*Bc (M={M})",
          d[4][-1] <= d[0][-1] <= M * Bc + 1e-9, f"bound {d[4][-1]:.2f} <= MC {d[0][-1]:.2f} <= M*Bc {M*Bc:.2f}")
    check("Limit phi=0", f"no blockage: EDMA CR below TDMA (M={M})", d[0][0] < d[6][0],
          f"EDMA {d[0][0]:.2f} vs TDMA {d[6][0]:.2f}; bound {d[4][0]:.2f}")

# ------------------------------------------------------------ report
w = max(len(c[1]) for c in rows)
lines = [f"| Where | Claim | Result | Value |", "|---|---|---|---|"]
for where, claim, res, val in rows:
    print(f"{res:4}  {where:16} {claim:<{w}}  {val}")
    lines.append(f"| {where} | {claim} | {'**FAIL**' if res == 'FAIL' else 'PASS'} | {val} |")
nf = sum(1 for x in rows if x[2] == "FAIL")
print(f"\n{len(rows) - nf}/{len(rows)} PASS, {nf} FAIL")
with open(os.path.join(HERE, "audit_report.md"), "w", encoding="utf-8") as f:
    f.write(f"# Paper audit ({len(rows) - nf}/{len(rows)} pass)\n\n" + "\n".join(lines) + "\n")
