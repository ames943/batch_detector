"""Shared style for all arXiv-version figures. Import and call `apply()` first."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Colorblind-safe palette (Okabe-Ito derived)
RESPONSE = "#0072B2"    # dark blue
COHORT = "#D55E00"      # orange/vermillion
CHANCE = "#7F7F7F"      # gray
GRID = "#dddddd"

# Marker/line-style redundancy for grayscale readability
RESPONSE_MARKER = "o"     # filled circle
COHORT_MARKER = "s"       # open square
RESPONSE_LINE = "-"
COHORT_LINE = (0, (3.5, 2))   # dashed
CHANCE_LINE = ":"

SINGLE_COL = 3.5
DOUBLE_COL = 7.1

# Qualitative, colorblind-safe (Okabe-Ito) palette for categorical groupings
# (individual cohorts / sites) — distinct from the RESPONSE/COHORT/CHANCE
# semantic-role colors above, which encode "which factor" not "which category".
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9", "#E69F00", "#000000"]


def apply():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 7.5,
        "axes.labelsize": 7.5,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "lines.linewidth": 1.1,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def panel_letter(ax, letter, x=-0.02, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="right")


def save(fig, path_noext):
    fig.savefig(f"{path_noext}.pdf")
    fig.savefig(f"{path_noext}.png", dpi=300)
    plt.close(fig)
    print(f"  Saved: {path_noext}.pdf / .png")
