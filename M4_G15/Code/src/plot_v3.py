"""Revised result figures (reviewer round 1): four single-column figures.
  fig_cr_vs_P  - ergodic sum CR vs P, M = 2 (hollow markers) and M = 5 (filled markers)
  fig_sr_vs_P  - ergodic sum SR vs P, same convention
  fig_cr_vs_M  - ergodic sum CR vs M, D_L = 10, 20, 40 m (marker shape)
  fig_sr_vs_M  - ergodic sum SR vs M
Saved to results/v3/<name>.pdf/.png.  Axis labels are drawn here (not by LaTeX); legend entries
follow the top-to-bottom order of the curves; a zoom inset (ellipse + arrow) is added only where
an EDMA curve overlaps or nearly overlaps another curve."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
from edma_sim import (power_sweep, target_sweep, m_sweep, ula_power_sweep, ula_m_sweep, P_DBM)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "v3")
os.makedirs(OUT, exist_ok=True)

# "--half": panels for two sub-figures side by side in one IEEE column (0.49\columnwidth each)
HALF = "--half" in sys.argv
SFX = "_h" if HALF else ""
if HALF:
    FS = dict(tick=7.2, label=8.0, leg=6.4, ins=5.6, lw=0.85, ms=3.0, cross=3.6, alw=0.6,
              mew=0.7, axlw=0.6)
    PANEL = (1.80, 1.85)
else:
    FS = dict(tick=8.5, label=9.5, leg=7.2, ins=6.8, lw=1.15, ms=4.6, cross=5.2, alw=0.8,
              mew=0.95, axlw=0.8)
    PANEL = (3.45, 2.65)                 # inches, IEEE single column
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": FS["tick"],
    "axes.labelsize": FS["label"], "xtick.labelsize": FS["tick"], "ytick.labelsize": FS["tick"],
    "legend.fontsize": FS["leg"], "axes.linewidth": FS["axlw"], "lines.linewidth": FS["lw"],
    "lines.markersize": FS["ms"], "lines.markeredgewidth": FS["mew"],
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 3.2 * FS["axlw"] / 0.8, "ytick.major.size": 3.2 * FS["axlw"] / 0.8,
    "xtick.minor.size": 1.8 * FS["axlw"] / 0.8, "ytick.minor.size": 1.8 * FS["axlw"] / 0.8,
    "axes.grid": True, "grid.color": "0.86", "grid.linewidth": 0.5 * FS["axlw"] / 0.8,
    "grid.linestyle": "-", "legend.edgecolor": "0.15", "legend.fancybox": False,
    "legend.framealpha": 1.0, "savefig.dpi": 600})
# One style per scheme, identical in every figure.
SCHEME = {
    "user": dict(label=r"EDMA, PA at $U_m$",         color="#1F4E9E", marker="^", ls="-"),
    "tgt":  dict(label=r"EDMA, PA at $T_m$",         color="#1F4E9E", marker="v", ls="-"),
    "mid":  dict(label="EDMA, PA at midpoint",  color="#D7191C", marker="o", ls="-"),
    "full": dict(label="EDMA, all interferers",      color="k",       marker="x", ls="none"),
    "lb":   dict(label="Lower bound (Theorem 1)",    color="#D7191C", marker="d", ls="-."),
    "lbs":  dict(label="Lower bound (Theorem 2)",    color="#D7191C", marker="d", ls="-."),
    "tdma": dict(label="TDMA",                       color="#1A8A3A", marker=">", ls="-"),
    "ula":  dict(label="Conventional ULA",           color="k",       marker="s", ls="--"),
}


if HALF:                                 # shorter labels for the small sub-figures
    SCHEME["lb"]["label"] = "Lower bound (Th. 1)"
    SCHEME["lbs"]["label"] = "Lower bound (Th. 2)"


def draw(ax, x, y, key, filled=False, markevery=(1, 2), **kw):
    """Plot one curve in the scheme's style; filled markers mark M = 5, hollow M = 2."""
    s = SCHEME[key]
    cross = s["marker"] in ("x", "+")
    st = dict(color=s["color"], ls=s["ls"], marker=s["marker"], markevery=markevery,
              markerfacecolor=s["color"] if (filled or cross) else "white",
              markeredgecolor=s["color"])
    if cross:
        st.update(markersize=FS["cross"])
    st.update(kw)
    return ax.plot(x, y, **st)[0]


def proxy(key):
    s = SCHEME[key]
    cross = s["marker"] in ("x", "+")
    return Line2D([], [], color=s["color"], ls=s["ls"], marker=s["marker"],
                  markerfacecolor=s["color"] if cross else "white", markeredgecolor=s["color"],
                  markersize=FS["cross"] if cross else FS["ms"], label=s["label"])


def legend_top(ax, handles, ncol=2):
    """Framed legend above the axes; entries fill column by column (top-to-bottom order)."""
    leg = ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.045), ncol=ncol,
                    handlelength=2.3, borderpad=0.4, labelspacing=0.25, handletextpad=0.5,
                    columnspacing=1.0, borderaxespad=0.0)
    leg.get_frame().set_linewidth(0.6)
    return leg


def zoom(ax, curves, xwin, ywin, rect, arrow_to, log=False, efac=(1.15, 1.7),
         ylabels_right=False):
    """Inset magnifying xwin x ywin (curves: callable drawing on the inset); an ellipse marks
    the magnified region on the main axes and an arrow points from it to the inset."""
    ins = ax.inset_axes(rect)
    curves(ins)
    ins.set_xlim(*xwin); ins.set_ylim(*ywin)
    if not log:
        from matplotlib.ticker import MaxNLocator
        ins.yaxis.set_major_locator(MaxNLocator(nbins=3))
    if log:
        ins.set_yscale("log")
        if ywin[1] / ywin[0] < 2:            # narrow window: 2-3 evenly spaced labels
            from matplotlib.ticker import FixedLocator, MaxNLocator
            tv = [v for v in MaxNLocator(nbins=3).tick_values(*ywin) if ywin[0] <= v <= ywin[1]]
            ins.yaxis.set_major_locator(FixedLocator(tv))
        else:
            ins.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 1.5, 2, 3, 4, 5, 6, 8)))
        ins.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ins.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.3g}"))
        ins.yaxis.set_minor_formatter(NullFormatter())
    ins.tick_params(labelsize=FS["ins"], length=2 * FS["axlw"] / 0.8, width=FS["axlw"] * 0.8,
                    which="both", labelbottom=False, pad=1.5)
    if ylabels_right:
        ins.tick_params(axis="y", labelleft=False, labelright=True)
    ins.grid(True, which="major", color="0.88", linewidth=0.4)
    for sp in ins.spines.values():
        sp.set_linewidth(0.7)
    to_ax = lambda xx, yy: ax.transAxes.inverted().transform(ax.transData.transform((xx, yy)))
    (x0, y0), (x1, y1) = to_ax(xwin[0], ywin[0]), to_ax(xwin[1], ywin[1])
    ex, ey, ew, eh = (x0 + x1) / 2, (y0 + y1) / 2, (x1 - x0) * efac[0], (y1 - y0) * efac[1]
    ax.add_patch(Ellipse((ex, ey), ew, eh, transform=ax.transAxes, fill=False, ec="k", lw=FS["alw"],
                         zorder=6))
    vx, vy = arrow_to[0] - ex, arrow_to[1] - ey          # arrow starts on the ellipse edge
    t = 1.0 / np.sqrt((vx / (ew / 2)) ** 2 + (vy / (eh / 2)) ** 2)
    ax.annotate("", xy=arrow_to, xytext=(ex + t * vx, ey + t * vy), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="k", lw=FS["alw"],
                                mutation_scale=8 * FS["alw"] / 0.8), zorder=6)
    return ins


# fixed axes box (inches: left, bottom, right, top margins) so every panel has identical
# outer size and identical plotting area, independent of tick-label widths
MARGIN = (0.44, 0.33, 0.06, 0.04) if HALF else (0.62, 0.45, 0.08, 0.06)


def bold_legend(leg):
    import re
    for t in leg.get_texts():
        t.set_fontweight("bold")
        t.set_text(re.sub(r"\$([^$]+)\$", lambda m: r"$\boldsymbol{" + m.group(1) + "}$",
                          t.get_text()))


def place(ax):
    W, H = ax.figure.get_size_inches()
    l, b, r, t = MARGIN
    ax.set_position([l / W, b / H, 1 - (l + r) / W, 1 - (b + t) / H])


def finish(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, labelpad=1.5)
    ax.set_ylabel(ylabel, labelpad=2)
    ax.set_axisbelow(True)
    place(ax)
    if ax.get_legend() is not None:
        bold_legend(ax.get_legend())


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + SFX + ".pdf"))
    fig.savefig(os.path.join(OUT, name + SFX + ".png"))
    plt.close(fig)


# ------------------------------------------------------------------ CR / SR vs P
def power_data(kind):
    data = {}
    for M in (2, 5):
        d = power_sweep(M); u = ula_power_sweep(M)
        if kind == "cr":
            data[M] = dict(mid=d[0], full=d[6], lb=d[4], tdma=d[2], ula=u[0])
        else:
            data[M] = dict(mid=d[1], full=d[7], lbs=d[5], tdma=d[3], ula=u[1])
    return data


def m_key():
    """One legend row explaining the marker fill: filled = M = 5, hollow = M = 2."""
    from matplotlib.legend_handler import HandlerTuple
    f = Line2D([], [], color="0.3", ls="none", marker="o", markerfacecolor="0.3")
    h = Line2D([], [], color="0.3", ls="none", marker="o", markerfacecolor="white")
    lab = "$M=5$ / $M=2$" if HALF else "$M=5$ / $M=2$ (filled / hollow)"
    return (f, h), lab, {tuple: HandlerTuple(ndivide=None, pad=0.6)}


def vs_power(kind):
    P = np.asarray(P_DBM, float)
    data = power_data(kind)
    keys = list(data[2])

    def plot_all(a, markevery=(1, 2)):
        for M in (5, 2):
            for k in keys:
                draw(a, P, data[M][k], k, filled=(M == 5), markevery=markevery)

    fig, ax = plt.subplots(figsize=PANEL)
    place(ax)
    plot_all(ax)
    ax.set_xlim(0, 30); ax.set_xticks(range(0, 31, 5))
    # legend inside, top-left; entries in the top-to-bottom order of the curves
    order = sorted(keys, key=lambda k: -max(data[M][k][-1] for M in (2, 5)))
    hk, lk, hmap = m_key()
    leg = ax.legend(handles=[proxy(k) for k in order] + [hk],
                    labels=[SCHEME[k]["label"] for k in order] + [lk], handler_map=hmap,
                    loc="upper left", fontsize=FS["leg"], handlelength=2.2, borderpad=0.35,
                    labelspacing=0.22, handletextpad=0.45, borderaxespad=0.45)
    leg.get_frame().set_linewidth(0.6)
    every = lambda a: plot_all(a, markevery=1)
    d2, d5 = data[2], data[5]
    if kind == "cr":
        ax.set_ylim(0, 120 if HALF else 80)
        # M = 2: midpoint, all-interferer and bound curves within ~0.2 bits/s/Hz
        zoom(ax, every, (25, 29), (d2["lb"][-3] - 0.35, d2["mid"][-2] + 0.35),
             rect=[0.16, 0.33, 0.3, 0.17] if HALF else [0.61, 0.66, 0.3, 0.31],
             arrow_to=(0.46, 0.38) if HALF else (0.78, 0.655), ylabels_right=not HALF)
    else:
        # log axis: the SR spans ~3 decades (two-way path loss makes it small at low P);
        # the empty decades on top hold the legend
        ax.set_yscale("log"); ax.set_ylim(3e-3, 1e4 if HALF else 3e2)
        # M = 2: midpoint, all-interferer and bound curves within ~3%
        zoom(ax, every, (25, 29), (d2["lbs"][-3] * 0.97, d2["mid"][-2] * 1.03),
             rect=[0.64, 0.045, 0.33, 0.2], arrow_to=(0.8, 0.25), log=True,
             efac=(1.15, 2.0))
    finish(ax, "Transmit power $P$ (dBm)",
           "Ergodic sum CR (bits/s/Hz)" if kind == "cr" else "Ergodic sum SR (bits/s/Hz)")
    save(fig, f"fig_{kind}_vs_P")


# ------------------------------------------------------------------ CR / SR vs M
DL_MARK = {40.0: "o", 20.0: "s", 10.0: "^"}


def vs_M(kind):
    Ms = np.arange(1, 9)
    r = 0 if kind == "cr" else 1
    curves = []                                   # (y, scheme, D_L)
    for DL in DL_MARK:
        d = m_sweep(DL); u = ula_m_sweep(DL)
        curves += [(d[r], "mid", DL), (d[r + 2], "tdma", DL), (u[r], "ula", DL)]

    def plot_all(a, keys=("mid", "tdma", "ula")):
        for yy, key, DL in curves:
            if key not in keys:
                continue
            s = SCHEME[key]
            a.plot(Ms, yy, color=s["color"], ls=s["ls"], marker=DL_MARK[DL],
                   markerfacecolor="white", markeredgecolor=s["color"])

    fig, ax = plt.subplots(figsize=PANEL)
    place(ax)
    plot_all(ax)
    hs = [Line2D([], [], color=SCHEME[k]["color"], ls=SCHEME[k]["ls"], label=lab)
          for k, lab in (("mid", "EDMA, PA at midpoint"), ("tdma", "TDMA"),
                         ("ula", "Conventional ULA"))]
    hd = [Line2D([], [], color="0.25", ls="none", marker=mk, markerfacecolor="white",
                 label=f"$D_L={D:g}$ m") for D, mk in DL_MARK.items()]
    # two stacked legends at the top left: schemes, then the D_L marker key right below it
    kw = dict(labelspacing=0.28, borderpad=0.4, handletextpad=0.4, fontsize=FS["leg"])
    leg1 = ax.legend(handles=hs, loc="upper left", handlelength=2.0, borderaxespad=0.4, **kw)
    bold_legend(leg1)
    ax.add_artist(leg1)
    fig.canvas.draw()
    bb = leg1.get_window_extent().transformed(ax.transAxes.inverted())
    leg2 = ax.legend(handles=hd, loc="upper left", bbox_to_anchor=(bb.x0, bb.y0 - 0.02),
                     handlelength=1.0, borderaxespad=0.0, **kw)
    bold_legend(leg2)
    for lg in (leg1, leg2):
        lg.get_frame().set_linewidth(0.6)
    ax.set_xlim(0.7, 8.3); ax.set_xticks(Ms)
    if kind == "cr":
        ax.set_ylim(0, 110 if HALF else 80)
        # small M: EDMA (all D_L) and TDMA within ~3 bits/s/Hz; EDMA (D_L = 10 m) crosses TDMA
        zoom(ax, lambda q: plot_all(q, ("mid", "tdma")), (0.85, 3.15), (12.6, 18.6),
             rect=[0.5, 0.27, 0.47, 0.15] if HALF else [0.5, 0.385, 0.47, 0.19],
             arrow_to=(0.5, 0.30) if HALF else (0.5, 0.43), efac=(1.12, 1.45))
    else:
        ax.set_ylim(0, 12 if HALF else 9.5)
    finish(ax, "Number of segments $M$",
           "Ergodic sum CR (bits/s/Hz)" if kind == "cr" else "Ergodic sum SR (bits/s/Hz)")
    save(fig, f"fig_{kind}_vs_M")


def preview():
    """2 x 2 overview of the four figures."""
    import matplotlib.image as mpimg
    names = ["fig_cr_vs_P", "fig_sr_vs_P", "fig_cr_vs_M", "fig_sr_vs_M"]
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 6.6))
    for a, n, cap in zip(axs.flat, names, ["Fig. 2", "Fig. 3", "Fig. 4", "Fig. 5"]):
        a.imshow(mpimg.imread(os.path.join(OUT, n + SFX + ".png"))); a.axis("off")
        a.set_title(cap, y=-0.08, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "preview_all" + SFX + ".png"), dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    for kind in ("cr", "sr"):
        vs_power(kind)
        vs_M(kind)
    preview()
    print("saved to", os.path.abspath(OUT))
