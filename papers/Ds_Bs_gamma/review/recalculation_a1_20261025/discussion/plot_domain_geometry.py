"""Illustrate the remaining domain slope at fixed a=1 and s_eff.

Geometry only; no sum-rule functions or production modules are imported.
"""
from pathlib import Path
import json
import os
import sys

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(HERE / ".cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np


def main():
    mass_a, mass_p, mass_q, s_eff = 2.53512, 1.96835, 1.27, 9.5
    mass_ratio = (mass_a / mass_p) ** 2
    floor = mass_q**2 / s_eff
    cases = [(1.0, "Equal intercepts", "#0047ab"),
             (mass_ratio, "Mass-proportional intercepts", "#b52025")]
    plt.rcParams.update({
        "font.family": "sans-serif", "mathtext.fontset": "dejavusans",
        "font.size": 11, "axes.linewidth": 1., "svg.fonttype": "path",
        "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 5.6), sharex=True, sharey=True)
    fig.subplots_adjust(left=.085, right=.98, bottom=.25, top=.79, wspace=.20)
    records = []
    for ax, (eta, title, color) in zip(axes, cases):
        xp, ya = (1+eta)/eta, 1+eta
        ax.axvspan(0, floor, color="0.94", zorder=0)
        ax.axhspan(0, floor, color="0.94", zorder=0)
        polygon = [(floor, floor), ((ya-floor)/eta, floor),
                   (floor, ya-eta*floor)]
        ax.add_patch(Polygon(polygon, closed=True, color=color, alpha=.14, lw=0))
        grid = np.linspace(0, 2.8, 200)
        ax.plot(grid, grid, color="0.5", ls="--", lw=1.)
        ax.plot([0, xp], [ya, 0], color=color, lw=2.6)
        ax.plot([floor, floor], [floor, ya-eta*floor], color="0.55", lw=.8)
        ax.plot([floor, (ya-floor)/eta], [floor, floor], color="0.55", lw=.8)
        ax.scatter([1], [1], s=38, c="black", zorder=5)
        ax.annotate("Fixed intersection\n$(1,1)$", xy=(1,1), xytext=(1.68,1.37),
                    fontsize=10, ha="center", arrowprops={"arrowstyle":"->","lw":.8})
        ax.text(.50, .47, "Retained\nregion", color=color, fontsize=11,
                ha="center", va="center")
        ax.text(2.30, 2.47, "$s_A=s_P$", color="0.45", fontsize=10,
                ha="center", rotation=45)
        ax.set_title(title+"\n"+rf"$\eta={eta:.3f}$", fontsize=12, pad=10)
        ax.set_xlim(0,2.8)
        ax.set_ylim(0,2.8)
        ax.set_aspect("equal")
        ax.set_xticks([0,1,2])
        ax.set_yticks([0,1,2])
        ax.tick_params(direction="in", top=True, right=True)
        ax.set_xlabel(r"Pseudoscalar variable $s_P/s_{\rm eff}$")
        records.append({"assignment":title,"eta":eta,"s_eff_GeV2":s_eff,
                        "intercept_P_GeV2":xp*s_eff,"intercept_A_GeV2":ya*s_eff,
                        "line_at_diagonal":ya-eta})
    axes[0].set_ylabel(r"Axial variable $s_A/s_{\rm eff}$")
    fig.suptitle(r"Same $a=1$ and $s_{\rm eff}$; different boundary slopes", y=.97, fontsize=15)
    fig.text(.5,.10,r"Example: $D_{s1}(2536)\to D_s\gamma$, $s_{\rm eff}=9.5\,\mathrm{GeV}^2$.",
             ha="center",fontsize=11)
    fig.text(.5,.045,r"Grey strips lie below the soft-kernel lower limit $s_P,s_A\geq m_c^2$.",
             ha="center",fontsize=10,color="0.35")
    for ext in ("png","svg"):
        fig.savefig(HERE/f"domain_geometry.{ext}",dpi=180)
    plt.close(fig)
    (HERE/"domain_geometry.json").write_text(json.dumps({
        "purpose":"Geometric illustration only; no OPE or width evaluated",
        "mass_A_GeV":mass_a,"mass_P_GeV":mass_p,"cases":records},indent=2)+"\n")
    print(json.dumps(records,indent=2))


if __name__ == "__main__":
    main()
