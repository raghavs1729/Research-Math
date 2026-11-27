#!/usr/bin/env python3
"""Figure: bin maxima of T_{q,a}(x) against the universal threshold W(x), q = 3 and 4.

From merged scan output (merge_scan.py): for each bin (0.01 in log10 x) the prime p with the
largest rho(p) = n^2 phi log p/(e p m), with its n, m.  By Proposition 2,
  T(p) - W(p) = G(p) phi^2 L^2 / p^{3/2} = (rho - 1) e m phi L / sqrt(p),
so the bin maximum of T - W is obtained exactly from (rho, p, m) (up to the variation of the
weight within a bin).  Failure <=> T >= W.
Usage: plot_T.py merged.txt out.png
"""
import sys, math, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import W_of
merged, outpng = sys.argv[1], sys.argv[2]
CPHI = [1, 2,2, 2,2, 4,4,4,4, 4,4,4,4, 4,4,4,4]
data = collections.defaultdict(list)
for line in open(merged):
    w = line.split()
    if not w or w[0] != "BIN": continue
    c, b, rho, p, n, m = int(w[1]), int(w[2]), float(w[6]), int(w[7]), int(w[8]), int(w[9])
    if b < 900 or m == 0 or p == 0: continue
    L = math.log(p)
    data[c].append((math.log10(p), (rho - 1) * math.e * m * CPHI[c] * L / math.sqrt(p) + W_of(p)))
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6), sharey=True, facecolor=SURF)
panels = [(4, 3, 4, "q = 4"), (3, 1, 2, "q = 3")]
xmax = max(x for c in (1, 2, 3, 4) for x, _ in data[c])
xmax = math.ceil(xmax * 2 - 1e-9) / 2
us = [9 + 0.01 * i for i in range(int(round((xmax - 9) / 0.01)) + 1)]
Wc = [W_of(10 ** u) for u in us]
ytop = min(max(Wc) * 1.05, 16)
for ax, (q, cres, cnon, title) in zip(axes, panels):
    ax.set_facecolor(SURF)
    for c, col in ((cnon, ORANGE), (cres, BLUE)):
        xs, ys = zip(*sorted(data[c]))
        ax.plot(xs, ys, color=col, lw=1.1, zorder=2, solid_joinstyle="round")
    ax.plot(us, Wc, color=INK, lw=2, zorder=3)
    ax.set_xlim(9, xmax + 0.02); ax.set_ylim(-4.2, ytop)
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.tick_params(length=0)
    ax.set_xlabel("log\u2081\u2080 x")
    a_res, a_non = (1, 3) if q == 4 else (1, 2)
    ax.set_title(title, loc="left", fontsize=12, color=INK, pad=24)
    # legend row under the panel title (text in ink, colour carried by the swatch)
    y0 = 1.035
    ax.plot([0.00, 0.035], [y0, y0], color=ORANGE, lw=3, transform=ax.transAxes, clip_on=False, solid_capstyle="butt")
    ax.text(0.045, y0, f"a = {a_non}, non-residue", transform=ax.transAxes, va="center", fontsize=9.5, color=INK2)
    ax.plot([0.30, 0.335], [y0, y0], color=BLUE, lw=3, transform=ax.transAxes, clip_on=False, solid_capstyle="butt")
    ax.text(0.345, y0, f"a = {a_res}, quadratic residue", transform=ax.transAxes, va="center", fontsize=9.5, color=INK2)
    ax.plot([0.66, 0.695], [y0, y0], color=INK, lw=3, transform=ax.transAxes, clip_on=False, solid_capstyle="butt")
    ax.text(0.705, y0, "threshold W(x)", transform=ax.transAxes, va="center", fontsize=9.5, color=INK2)
axes[0].set_ylabel("largest T(x) in each bin of width 0.01 in log\u2081\u2080 x")
fig.suptitle("Ramanujan's inequality mod q fails exactly where T(x) reaches W(x)", x=0.012, ha="left", fontsize=13, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(outpng, dpi=160, facecolor=SURF)
print("wrote", outpng)
