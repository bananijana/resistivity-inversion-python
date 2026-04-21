"""
ERT / VES Resistivity Inversion Using pyGIMLi
-----------------------------------------------
Module 1 — Vertical Electrical Sounding (VES)
  - 1D Schlumberger sounding data (synthetic)
  - pyGIMLi VESManager inversion
  - Apparent resistivity curve
  - 1D resistivity-depth model
  - Layer interpretation

Module 2 — 2D ERT Profile
  - Synthetic Wenner array ERT data
  - pyGIMLi ERTManager inversion
  - Pseudosection (apparent resistivity)
  - Inverted 2D resistivity section
  - Sensitivity distribution

Geological context:
  Simulates a shallow aquifer investigation scenario —
  resistivity contrast between saturated alluvial sand,
  clay aquitard, and weathered/fractured basement.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import LogLocator, LogFormatter
import pygimli as pg
import pygimli.physics.ert as ert
from pygimli.physics import VESManager
import pygimli.meshtools as mt
import warnings
import os

warnings.filterwarnings("ignore")
pg.setVerbose(False)

OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 1 — VERTICAL ELECTRICAL SOUNDING (VES)
# ═══════════════════════════════════════════════════════════════════════════════

print("── Module 1: VES Inversion ──")

# ── Synthetic VES data (Schlumberger) ─────────────────────────────────────────
# AB/2 spacings (m)
ab2 = np.array([1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0,
                15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0,
                150.0, 200.0, 300.0])

# Synthetic apparent resistivity (Ω·m) — KH-type curve
# Simulates: topsoil (high) → clay (low) → sand aquifer (high) → basement (high)
rho_app = np.array([85, 72, 65, 52, 44, 40, 38, 42,
                    55, 68, 90, 110, 125, 145, 160,
                    180, 195, 210], dtype=float)

# Add ±5% noise
np.random.seed(14)
rho_app *= (1 + np.random.uniform(-0.05, 0.05, len(rho_app)))
rho_app  = np.round(rho_app, 2)

# Error vector (5%)
err = rho_app * 0.05

# Save VES data
ves_df = pd.DataFrame({"AB2_m": ab2, "Rho_app_ohmm": rho_app,
                        "Error_ohmm": err})
ves_df.to_csv(f"{OUTPUT_DIR}/ves_data.csv", index=False)
print(f"  AB/2 range: {ab2.min()} – {ab2.max()} m")
print(f"  Rho_app range: {rho_app.min():.1f} – {rho_app.max():.1f} Ω·m")

# ── VES Inversion ─────────────────────────────────────────────────────────────
mgr = VESManager()

# Starting model: 4 layers
# thicknesses (m), resistivities (Ω·m)
startModel = np.array([2.0, 5.0, 15.0,        # thicknesses
                        80.0, 15.0, 120.0, 200.0])  # resistivities

mn2 = ab2 / 3.0
inv_result = mgr.invert(data=rho_app, err=err,
                         ab2=ab2, mn2=mn2,
                         startModel=startModel,
                         lam=20, verbose=False)

print(f"  Inversion chi² = {mgr.inv.chi2():.3f}")

# Extract model — startModel has 3 thicknesses + 4 resistivities
n_layers = 4
model_arr = np.array(mgr.model)
thicknesses = model_arr[:n_layers - 1]
resistivities = model_arr[n_layers - 1:]

# Build depth array for step plot
max_inv_depth = float(np.sum(thicknesses) * 1.5)
depths     = np.concatenate([[0], np.cumsum(thicknesses), [max_inv_depth]])
depth_plot = np.repeat(depths, 2)[1:-1]
rho_plot   = np.repeat(resistivities, 2)
max_depth  = max_inv_depth

# Compute fit
rho_fit = mgr.inv.response

# ── Plot 1: VES Apparent Resistivity + Model Fit ──────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 6))

# Left — apparent resistivity curve
ax = axes[0]
ax.loglog(ab2, rho_app, "o", color="#1565C0", markersize=7,
          label="Observed", zorder=3)
ax.loglog(ab2, rho_fit, "-", color="#D32F2F", linewidth=2,
          label=f"Calculated  χ²={mgr.inv.chi2():.2f}", zorder=2)
ax.fill_between(ab2, rho_app - err, rho_app + err,
                alpha=0.2, color="#1565C0", label="±5% error")
ax.set_xlabel("AB/2 (m)", fontsize=11)
ax.set_ylabel("Apparent Resistivity (Ω·m)", fontsize=11)
ax.set_title("VES — Apparent Resistivity Curve\n(Schlumberger Array)",
             fontsize=11, fontweight="bold")
ax.legend(fontsize=9)
ax.grid(True, which="both", alpha=0.3)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Right — 1D resistivity-depth model
ax2 = axes[1]
ax2.semilogx(rho_plot, depth_plot, color="#2E7D32", linewidth=2.5)
ax2.fill_betweenx(depth_plot, 0, rho_plot, alpha=0.2, color="#2E7D32")
ax2.invert_yaxis()
ax2.set_xlabel("Resistivity (Ω·m)", fontsize=11)
ax2.set_ylabel("Depth (m)", fontsize=11)
ax2.set_title("Inverted 1D Resistivity–Depth Model",
              fontsize=11, fontweight="bold")
ax2.set_ylim(max_depth, 0)
ax2.grid(True, which="both", alpha=0.3)

# Annotate layers
layer_labels = ["Topsoil / Dry alluvium",
                "Clay / Silt aquitard",
                "Saturated sand aquifer",
                "Weathered basement"]
layer_colors = ["#8D6E63", "#A0785A", "#42A5F5", "#78909C"]
mid_depths   = [(depths[i] + depths[i+1])/2
                for i in range(len(depths)-1)]
mid_depths.append(depths[-1] * 1.15)

for lbl, col, mid in zip(layer_labels, layer_colors, mid_depths):
    ax2.axhline(mid, color=col, linewidth=0, alpha=0)
    ax2.text(resistivities.max() * 1.05, mid, lbl,
             fontsize=8, color=col, va="center", fontstyle="italic")

for d in depths[1:-1]:
    ax2.axhline(d, color="grey", linewidth=0.8,
                linestyle="--", alpha=0.6)

ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

fig.suptitle("VES Inversion — 1D Subsurface Resistivity Model",
             fontsize=13, fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/ves_inversion.png", dpi=180,
            bbox_inches="tight")
plt.close()
print("Saved: ves_inversion.png")

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 2 — 2D ERT PROFILE
# ═══════════════════════════════════════════════════════════════════════════════

print("\n── Module 2: 2D ERT Inversion ──")

# ── Synthetic 2D ERT geometry (Wenner array) ──────────────────────────────────
# 24 electrodes, 5 m spacing → 115 m profile
n_elec   = 24
spacing  = 5.0
elec_x   = np.arange(n_elec) * spacing
elec_pos = np.column_stack([elec_x, np.zeros(n_elec)])

# Create scheme (Wenner)
scheme = ert.createData(elec_pos, schemeName="wa")
print(f"  Electrodes: {n_elec}  |  Spacing: {spacing} m")
print(f"  Profile length: {elec_x.max():.0f} m")
print(f"  Data points: {scheme.size()}")

# ── Forward model — synthetic subsurface ─────────────────────────────────────
# Build geometry: clay layer over sand aquifer over basement
world  = mt.createWorld(start=[-5, 0], end=[elec_x.max()+5, -40],
                         worldMarker=True)
layer1 = mt.createRectangle(start=[-5, -0], end=[elec_x.max()+5, -5],
                             marker=2)   # clay
layer2 = mt.createRectangle(start=[-5, -5], end=[elec_x.max()+5, -20],
                             marker=3)   # sand aquifer
# Add a resistive anomaly (dry zone / buried structure)
anomaly = mt.createCircle(pos=[60, -12], radius=6, marker=4)

geom = world + layer1 + layer2 + anomaly
mesh_fwd = mt.createMesh(geom, quality=34, area=2)

# Assign resistivities
rhomap = [[1,  200],   # background basement
          [2,  15],    # clay
          [3,  80],    # saturated sand
          [4,  350]]   # resistive anomaly

# Simulate data
mgr2  = ert.ERTManager(sr=False, verbose=False)
data  = mgr2.simulate(mesh_fwd, scheme=scheme,
                       res=rhomap, noiseLevel=0.03,
                       noiseAbs=1e-6, seed=42)

print(f"  Simulated data range: "
      f"{data['rhoa'].array().min():.1f} – "
      f"{data['rhoa'].array().max():.1f} Ω·m")

# ── Inversion ─────────────────────────────────────────────────────────────────
mgr2.invert(data, lam=30, verbose=False)
print(f"  Inversion chi² = {mgr2.inv.chi2():.3f}")

# ── Plot 2: Pseudosection + Inverted Section ──────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# Top — apparent resistivity pseudosection
pg.show(data, data["rhoa"], ax=axes[0],
        label="Apparent Resistivity (Ω·m)",
        cMap="jet", logScale=True,
        xlabel="Distance (m)", ylabel="Pseudo-depth (m)",
        orientation="vertical")
axes[0].set_title("Apparent Resistivity Pseudosection (Wenner Array)",
                  fontsize=12, fontweight="bold")

# Bottom — inverted resistivity section
pg.show(mgr2.paraDomain, mgr2.model, ax=axes[1],
        label="Resistivity (Ω·m)",
        cMap="jet", logScale=True,
        xlabel="Distance (m)", ylabel="Depth (m)",
        orientation="vertical")
axes[1].set_title("Inverted 2D Resistivity Section",
                  fontsize=12, fontweight="bold")

# Annotate geology
axes[1].text(5, -2.5,  "Clay / Silt",
             fontsize=9, color="white", fontweight="bold")
axes[1].text(5, -12,   "Saturated Sand Aquifer",
             fontsize=9, color="white", fontweight="bold")
axes[1].text(5, -30,   "Weathered Basement",
             fontsize=9, color="white", fontweight="bold")
axes[1].text(57, -12,  "Resistive\nAnomaly",
             fontsize=8, color="white", fontweight="bold", ha="center")

fig.suptitle("2D ERT Survey — Pseudosection and Inverted Resistivity Section",
             fontsize=13, fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/ert_2d_section.png", dpi=180,
            bbox_inches="tight")
plt.close()
print("Saved: ert_2d_section.png")

# ── Plot 3: Coverage / Sensitivity ────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 5))
pg.show(mgr2.paraDomain, np.log10(mgr2.coverage()),
        ax=ax, cMap="plasma",
        label="log₁₀(Coverage)",
        xlabel="Distance (m)", ylabel="Depth (m)",
        orientation="vertical")
ax.set_title("ERT Inversion — Model Coverage / Sensitivity",
             fontsize=12, fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/ert_coverage.png", dpi=180,
            bbox_inches="tight")
plt.close()
print("Saved: ert_coverage.png")

# ── Plot 4: Forward model (true resistivity) ──────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 5))
pg.show(mesh_fwd, data=pg.solver.parseMapToCellArray(rhomap, mesh_fwd),
        ax=ax, cMap="jet", logScale=True,
        label="True Resistivity (Ω·m)",
        xlabel="Distance (m)", ylabel="Depth (m)",
        orientation="vertical")
ax.set_title("True Subsurface Resistivity Model (Forward Model)",
             fontsize=12, fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/ert_true_model.png", dpi=180,
            bbox_inches="tight")
plt.close()
print("Saved: ert_true_model.png")

# ── Save ERT data ─────────────────────────────────────────────────────────────
ert_df = pd.DataFrame({
    "a": np.array(data["a"]),
    "b": np.array(data["b"]),
    "m": np.array(data["m"]),
    "n": np.array(data["n"]),
    "rhoa": np.array(data["rhoa"]),
    "err":  np.array(data["err"]),
})
ert_df.to_csv(f"{OUTPUT_DIR}/ert_data.csv", index=False)
print("Saved: ert_data.csv")

print(f"\n✓ All outputs saved to: {OUTPUT_DIR}")
print(f"  Module 1 — VES: 1 plot + 1 CSV")
print(f"  Module 2 — 2D ERT: 3 plots + 1 CSV")
