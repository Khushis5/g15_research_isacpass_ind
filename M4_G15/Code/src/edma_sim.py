"""
Shared simulation + closed-form module for the four figures of the letter
(ISAC-PASS under EDMA geometry).  Model = Section II of main_ver2.tex:
  * M isolated waveguide segments of width S = D_L/M; PA_m at the midpoint of U_m and T_m
    (no in-waveguide propagation loss -- the antenna radiates at full free-space strength
    regardless of where it sits along the segment)
  * own links LoS; cross links blocked with P(LoS) = exp(-phi r^2)
  * PERFORMANCE-ANALYSIS MODEL (this file): users/targets on y = 0 and interference / clutter /
    leakage only from the two NEIGHBOURING segments m-1, m+1 (Section II of the paper keeps
    the general model with all segments and (x, y, 0) positions)
  * sensing: clutter (mu) + leakage (nu/L) over the same PA_m--T_i link
  * TDMA baseline: TWO isolated waveguides of length D_L (Section III-C).
      - Transmit waveguide: M pinching antennas FIXED at the centre of each segment,
        tau_m = (m + 1/2) S (they do not track U_m or T_m).
      - Receive waveguide: ONE pinching antenna FIXED at the centre of the whole
        waveguide, x_r = D_L/2, used only to collect the sensing echo.
      - Slot m: only PA_m (Tx) is active, with power M*P (equal per-chain average power).
        Communication is the ordinary point-to-point link Tx_m -> U_m. Sensing is now
        BISTATIC: Tx_m illuminates T_m, and the fixed Rx antenna collects the echo, so the
        round trip uses two different legs, Tx_m->T_m and T_m->Rx, instead of one leg twice.
      - No interference, clutter, or leakage (only one Tx/one target active per slot).
Run any figN script; Monte Carlo results are cached in sim_cache_final.npz under a key that
contains a hash of ALL model parameters, so a parameter change can never reuse stale data.
"""
import itertools, math, cmath, os, hashlib
import numpy as np

FC, C0 = 28e9, 3e8
ETA = (C0 / (4 * np.pi * FC)) ** 2
D, D2 = 3.0, 9.0                      # waveguide height d, d^2 (m)
SIG2 = 10 ** ((-90 - 30) / 10)        # noise power, W (-90 dBm) for both CR and SR
ALPHA_S, L = 8237.0, 5                # alpha_s = 4*pi*sigma/lambda^2, AGV RCS sigma = -11.235 dBsm at 28 GHz; frame length
PHI = 0.02                            # blockage severity
N_MC = 100_000                        # Monte Carlo realisations
MODEL_VERSION = "final-v2-baselines"            # bump when the model equations change
_PARAMS = repr((FC, C0, D, SIG2, ALPHA_S, L, PHI, N_MC, MODEL_VERSION))
TAG = hashlib.md5(_PARAMS.encode()).hexdigest()[:8]
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sim_cache_final.npz")

# ------------------------------------------------------------ Monte Carlo
def _geometry(D_L, M, n, seed):
    """x-positions measured from the left end of the waveguide (segment m spans
    [m*S, (m+1)*S)). Returns xU, xT (EDMA user/target draws) and t (EDMA antenna,
    the midpoint of xU, xT)."""
    rng = np.random.default_rng(seed)
    S = D_L / M
    st = np.arange(M) * S                          # left edge of each segment
    xU = st + S * rng.random((n, M)); xT = st + S * rng.random((n, M))
    t = 0.5 * (xU + xT)
    return xU, xT, t

def _neighbours(m, M, nearest):
    """Interferers of segment m: the two neighbours m-1, m+1 (nearest=True) or all other segments."""
    idx = (m - 1, m + 1) if nearest else range(M)
    return np.array([i for i in idx if 0 <= i < M and i != m], dtype=int)

def edma_mc(D_L, M, P, phi=PHI, n=N_MC, seed=2718, nearest=True, place="mid"):
    """Ergodic sum CR and SR of EDMA; blockage is averaged exactly over the 2^(#interferers)
    patterns. nearest=False -> interference/clutter/leakage from ALL other segments.
    place="mid": PA_m at the user-target midpoint (proposed); place="user": PA_m at x_U,m
    (EDMA placement of Ding et al.), which then also senses T_m monostatically."""
    xU, xT, t = _geometry(D_L, M, n, seed)
    if place == "user":
        t = xU.copy()
    elif place == "target":
        t = xT.copy()
    rho, rs = P / SIG2, P * L / SIG2
    dU = (t[:, :, None] - xU[:, None, :]) ** 2 + D2      # PA_i -> U_m  [i, m]
    dT = (t[:, :, None] - xT[:, None, :]) ** 2 + D2      # PA_m -> T_i  [m, i]
    Rc = np.zeros((n, M)); Rs = np.zeros((n, M))
    for m in range(M):
        nb = _neighbours(m, M, nearest)
        own_c = rho * ETA / dU[:, m, m]
        own_s = rs * ALPHA_S * ETA ** 2 / dT[:, m, m] ** 2
        ic = rho * ETA / dU[:, nb, m]; pc = np.exp(-phi * dU[:, nb, m])
        dmi, dii = dT[:, m, nb], dT[:, nb, nb]
        mu = ALPHA_S * ETA ** 2 / dmi ** 2                # clutter gain
        nu = ALPHA_S * ETA ** 2 / (dii * dmi)             # leakage gain
        isn = rs * (mu + nu / L); ps = np.exp(-phi * dmi)
        ac = np.zeros(n); as_ = np.zeros(n)
        for pat in itertools.product([0, 1], repeat=len(nb)):
            pat = np.array(pat)
            wc = np.prod(np.where(pat == 1, pc, 1 - pc), axis=1)
            ws = np.prod(np.where(pat == 1, ps, 1 - ps), axis=1)
            ac += wc * np.log2(1 + own_c / (1 + (ic * pat).sum(1)))
            as_ += ws * np.log2(1 + own_s / (1 + (isn * pat).sum(1)))
        Rc[:, m], Rs[:, m] = ac, as_ / L
    return Rc.sum(1).mean(), Rs.sum(1).mean()

def tdma_mc(D_L, M, P, n=N_MC, seed=2718):
    """TDMA with two fixed-antenna waveguides (Section III-C):
    Tx_m fixed at the centre of segment m; Rx fixed at the centre of the whole
    waveguide. Slot power M*P (equal per-chain average power). No interference/clutter;
    sensing is bistatic (Tx_m -> T_m -> Rx, two different legs)."""
    rng = np.random.default_rng(seed)
    S = D_L / M
    st = np.arange(M) * S
    xU = st + S * rng.random((n, M)); xT = st + S * rng.random((n, M))
    tau = st + S / 2                                    # fixed Tx antenna, one per segment
    xr = D_L / 2                                        # fixed Rx antenna (single, shared)
    rho, rs = M * P / SIG2, M * P * L / SIG2
    gc = rho * ETA / ((tau - xU) ** 2 + D2)
    r_tx2 = (tau - xT) ** 2 + D2                        # Tx_m -> T_m
    r_rx2 = (xr - xT) ** 2 + D2                         # T_m -> Rx (bistatic leg)
    gs = rs * ALPHA_S * ETA ** 2 / (r_tx2 * r_rx2)
    return np.log2(1 + gc).mean(), (np.log2(1 + gs) / L).mean()

# ------------------------------------------------------- closed-form bound
def I(a, S):
    """I(a) = E_w[log2(w^2 + a)], Eq. (I_def); a may be complex."""
    a = complex(a); c = S / 2
    return ((1 - 4 * a / S ** 2) * cmath.log(c ** 2 + a) + 4 * a / S ** 2 * cmath.log(a)
            - 3 + 8 * cmath.sqrt(a) / S * cmath.atan(S / (2 * cmath.sqrt(a)))) / math.log(2)

def p_k(phi, S, k):
    """Mean cross-link LoS probability, Eq. (ebar_closed).  theta = phi*S^2 > 0."""
    th, q = phi * S * S, k - 1
    if th < 1e-12:                                   # phi -> 0: every link is LoS
        return math.exp(-phi * D2)
    L0 = lambda x: 0.5 * math.sqrt(math.pi / th) * math.erf(math.sqrt(th) * x)
    L1 = lambda x: (1 - math.exp(-th * x * x)) / (2 * th)
    L2 = lambda x: (L0(x) - x * math.exp(-th * x * x)) / (2 * th)
    Ls = [L0, L1, L2]; tot = 0.0
    for (g2, g1, g0), (lo, hi) in zip([(2, 0, 0), (-2, 4, -1), (2, -8, 8)],
                                       [(0, .5), (.5, 1.5), (1.5, 2)]):
        a = [g2 * q * q - g1 * q + g0, g1 - 2 * g2 * q, g2]
        tot += sum(a[j] * (Ls[j](hi + q) - Ls[j](lo + q)) for j in range(3))
    return math.exp(-phi * D2) * tot

def lb(kind, P, D_L, M, phi=PHI):
    """Lower bounds of Theorems 1/2 for general M under nearest-neighbour interference
    (equal to the explicit M=2, M=5 forms): segment m pays |{m-1,m+1} within 1..M| * p_1."""
    S = D_L / M
    if kind == "c":
        rho = P / SIG2
        cl = (I(D2 + rho * ETA, S) - I(D2, S)).real
        mx = math.log2(1 + rho * ETA / D2)
    else:
        rs = P * L / SIG2; K = ETA * math.sqrt(rs * ALPHA_S)
        cl = (2 * I(D2 + 1j * K, S).real - 2 * I(D2, S).real) / L
        mx = math.log2(1 + rs * ALPHA_S * ETA ** 2 / D2 ** 2) / L
    p1 = p_k(phi, S, 1) if M > 1 else 0.0
    return sum(max(cl - mx * len(_neighbours(m, M, True)) * p1, 0.0) for m in range(M))

# ------------------------------------------------------------------ cache
def _load():
    return dict(np.load(CACHE, allow_pickle=True)) if os.path.exists(CACHE) else {}

def cached(key, fn):
    d = _load()
    key = f"{key}_{TAG}"
    if key not in d:
        d[key] = np.array(fn()); np.savez(CACHE, **d)
    return d[key]

P_DBM = np.arange(0, 31, 2)
def power_sweep(M, D_L=40.0):
    """rows: EDMA CR, EDMA SR, TDMA CR, TDMA SR, LB CR, LB SR  vs P (dBm)."""
    def run():
        out = []
        for pd in P_DBM:
            P = 10 ** ((pd - 30) / 10)
            ec, es = edma_mc(D_L, M, P); tc, ts = tdma_mc(D_L, M, P)
            fc, fs = edma_mc(D_L, M, P, nearest=False)
            uc, us = edma_mc(D_L, M, P, place="user")
            out.append([ec, es, tc, ts, lb("c", P, D_L, M), lb("s", P, D_L, M), fc, fs, uc, us])
        return np.array(out).T
    return cached(f"pow_M{M}_DL{D_L:g}", run)

def target_sweep(M, D_L=40.0):
    """rows: EDMA CR, EDMA SR with PA_m at the target T_m, vs P (dBm)."""
    def run():
        out = []
        for pd in P_DBM:
            P = 10 ** ((pd - 30) / 10)
            out.append(list(edma_mc(D_L, M, P, place="target")))
        return np.array(out).T
    return cached(f"tgt_M{M}_DL{D_L:g}", run)

def m_sweep_phi(D_L, phi, Ms=range(1, 9), P=1.0):
    """rows: EDMA CR, EDMA SR, TDMA CR, TDMA SR vs M for a given phi."""
    def run():
        out = []
        for M in Ms:
            ec, es = edma_mc(D_L, M, P, phi=phi); tc, ts = tdma_mc(D_L, M, P)
            out.append([ec, es, tc, ts])
        return np.array(out).T
    return cached(f"Mphi_DL{D_L:g}_phi{phi:g}", run)

def m_sweep(D_L, Ms=range(1, 9), P=1.0):
    """rows: EDMA CR, EDMA SR, TDMA CR, TDMA SR, LB CR, LB SR vs M at P = 30 dBm (1 W)."""
    def run():
        out = []
        for M in Ms:
            ec, es = edma_mc(D_L, M, P); tc, ts = tdma_mc(D_L, M, P)
            uc, us = edma_mc(D_L, M, P, place="user")
            out.append([ec, es, tc, ts, lb("c", P, D_L, M), lb("s", P, D_L, M), uc, us])
        return np.array(out).T
    return cached(f"M_DL{D_L:g}", run)

PHIS = np.logspace(-3, -1, 13)
def phi_sweep(M, D_L=40.0, P=1.0):
    """rows: EDMA CR, EDMA SR, TDMA CR, TDMA SR, LB CR, LB SR, EDMA-all CR, EDMA-all SR,
    PA-at-user CR, PA-at-user SR, p1  vs phi, at P = 30 dBm."""
    def run():
        out = []
        tc, ts = tdma_mc(D_L, M, P)
        for ph in PHIS:
            ec, es = edma_mc(D_L, M, P, phi=ph)
            fc, fs = edma_mc(D_L, M, P, phi=ph, nearest=False)
            uc, us = edma_mc(D_L, M, P, phi=ph, place="user")
            out.append([ec, es, tc, ts, lb("c", P, D_L, M, ph), lb("s", P, D_L, M, ph),
                        fc, fs, uc, us, p_k(ph, D_L / M, 1)])
        return np.array(out).T
    return cached(f"phi_M{M}_DL{D_L:g}", run)

# ------------------------------------------------- conventional-antenna (ULA) baseline
LAM = C0 / FC
ULA_LOS = True        # BS links of the conventional ULA are LoS (own links are LoS in EDMA too)

def ula_mc(D_L, M, P_list, phi=PHI, n=N_MC, seed=2718, x0=0.0):
    """Conventional-antenna benchmark: M antennas on top of the BS, a half-wavelength ULA
    at the BS, i.e., centred at (x0, 0, d) = (0, 0, d), the feed end of the waveguide
    (as in the network-model figure), along the waveguide direction. Antenna m transmits U_m's
    signal / the probe for T_m with power P (no precoding) and receives the echo of T_m
    (monostatic), i.e., the same SINR forms as Eqs. (sinr_comm_rate), (sensing_rate_closed)
    with the PA positions replaced by the antenna positions. The antennas are co-located,
    so every link to U_m (T_m) shares one blockage state, Pr = exp(-phi r^2); a blocked user
    or target gets zero rate. The blockage expectation is evaluated exactly over all
    2^M patterns. With ULA_LOS (paper setting) all BS links are LoS, like the EDMA own links,
    i.e., phi = 0 here. Returns arrays (CR, SR) over P_list."""
    if ULA_LOS:
        phi = 0.0
    rng = np.random.default_rng(seed)
    S = D_L / M
    st = np.arange(M) * S
    xU = st + S * rng.random((n, M)); xT = st + S * rng.random((n, M))
    xa = x0 + (np.arange(M) - (M - 1) / 2) * LAM / 2              # antenna x-positions
    P_arr = np.atleast_1d(np.asarray(P_list, float))
    rho, rs = P_arr / SIG2, P_arr * L / SIG2
    dU = (xa[None, :, None] - xU[:, None, :]) ** 2 + D2            # antenna i -> U_m [i, m]
    dT = (xa[None, :, None] - xT[:, None, :]) ** 2 + D2            # antenna m -> T_i [m, i]
    pU = np.exp(-phi * dU.mean(1)); pT = np.exp(-phi * dT.mean(1)) # shared blockage (n, M)
    cr = np.zeros(len(P_arr)); sr = np.zeros(len(P_arr))
    for pat in itertools.product([0, 1], repeat=M):
        pat = np.array(pat, bool)
        if not pat.any():
            continue
        wU = np.prod(np.where(pat, pU, 1 - pU), axis=1)            # (n,)
        wT = np.prod(np.where(pat, pT, 1 - pT), axis=1)
        for m in np.flatnonzero(pat):
            nb = np.array([i for i in np.flatnonzero(pat) if i != m], dtype=int)
            gc = ETA / dU[:, m, m]; ic = (ETA / dU[:, nb, m]).sum(1)
            dmi, dii = dT[:, m, nb], dT[:, nb, nb]
            gs = ALPHA_S * ETA ** 2 / dT[:, m, m] ** 2
            cs = (ALPHA_S * ETA ** 2 / dmi ** 2 + ALPHA_S * ETA ** 2 / (dii * dmi) / L).sum(1)
            rc = np.log2(1 + rho[:, None] * gc / (1 + rho[:, None] * ic))
            rsn = np.log2(1 + rs[:, None] * gs / (1 + rs[:, None] * cs)) / L
            cr += (rc * wU).mean(1); sr += (rsn * wT).mean(1)
    return cr, sr

def ula_power_sweep(M, D_L=40.0):
    """rows: ULA CR, ULA SR vs P (dBm)."""
    def run():
        cr, sr = ula_mc(D_L, M, 10 ** ((P_DBM - 30) / 10))
        return np.vstack([cr, sr])
    return cached(f"ulaLoS_pow_M{M}_DL{D_L:g}", run)

def ula_m_sweep(D_L, Ms=range(1, 9), P=1.0):
    """rows: ULA CR, ULA SR vs M at P = 30 dBm."""
    def run():
        return np.array([[v[0] for v in ula_mc(D_L, M, [P])] for M in Ms]).T
    return cached(f"ulaLoS_M_DL{D_L:g}", run)

# ------------------------------------------------------------- plot style
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
BLUE, ORNG, GREY = "#0B5FA5", "#D9480F", "#3B3B3B"
def style():
    mpl.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8.5,
        "axes.labelsize": 9.5, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "legend.fontsize": 7.6, "axes.linewidth": 0.8, "lines.linewidth": 1.5,
        "lines.markersize": 4.6, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "grid.linewidth": 0.4,
        "grid.alpha": 0.35, "savefig.dpi": 400})
def new_fig(): return plt.subplots(figsize=(3.55, 2.9))
def save(fig, name):
    os.makedirs("figs", exist_ok=True)
    fig.savefig(f"figs/{name}.png", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"figs/{name}.pdf", bbox_inches="tight", pad_inches=0.02)
