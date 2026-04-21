# resistivity-inversion-python

![Python](https://img.shields.io/badge/Python-3.8+-blue)
![pyGIMLi](https://img.shields.io/badge/pyGIMLi-1.5+-orange)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Complete-brightgreen)

A Python workflow for subsurface resistivity inversion and imaging using pyGIMLi. Covers both 1D sounding interpretation (VES) and 2D profile imaging (ERT), applied to a groundwater exploration scenario in a shallow alluvial setting.

---

## Overview

Electrical resistivity surveys are routinely used in hydrogeological investigations to delineate aquifer geometry, map clay layers, and identify buried structures — all without any drilling. This project implements a full Python-based processing and inversion pipeline using pyGIMLi, the open-source geophysical inversion library.

The workflow is structured around a practical groundwater exploration problem: identifying a saturated sand aquifer sandwiched between a surface clay layer and weathered basement, with a resistive anomaly present at depth — a setting directly relevant to shallow aquifer investigations in fluvial and deltaic environments.

---

## Why pyGIMLi?

pyGIMLi (Python Library for Inversion and Modelling) is an open-source framework specifically designed for geophysical forward modelling and inversion. It provides:
- Robust regularised inversion (Tikhonov regularisation)
- Flexible mesh generation for complex geometries
- Built-in ERT and VES managers
- Full control over inversion parameters

---

## Modules

### Module 1 — Vertical Electrical Sounding (VES)
Schlumberger array 1D sounding inversion using `VESManager`. Resolves a 4-layer resistivity model from apparent resistivity measurements across AB/2 spacings of 1–300 m.

### Module 2 — 2D ERT Profile
Wenner array 2D tomography using `ERTManager`. Forward simulation on a mesh with four resistivity units, followed by smooth-model inversion and sensitivity analysis.

---

## Geological Context

The simulated subsurface represents conditions typical of shallow aquifer systems in alluvial plains:

| Layer | Depth | Resistivity | Hydrogeological Role |
|---|---|---|---|
| Topsoil / dry alluvium | 0–2 m | ~80 Ω·m | Vadose zone |
| Clay / silt aquitard | 2–7 m | ~15 Ω·m | Confining layer |
| Saturated sand aquifer | 7–22 m | ~80–120 Ω·m | Target aquifer |
| Weathered basement | >22 m | ~200 Ω·m | Basement |
| Resistive anomaly | ~12 m depth | ~350 Ω·m | Buried feature |

---

## Inversion Setup

**VES:**
- Array: Schlumberger
- AB/2 range: 1–300 m (18 measurements)
- Layers: 4
- Regularisation λ = 20
- Noise: 5%

**ERT:**
- Array: Wenner
- Electrodes: 24 × 5 m spacing → 115 m profile
- Data points: 84
- Regularisation λ = 30
- Noise: 3%

---

## Project Structure

```
resistivity-inversion-python/
│
├── ert_ves_inversion.py     # VES + ERT inversion script
└── README.md
```

---

## Dependencies

| Library | Version | Purpose |
|---|---|---|
| `pygimli` | ≥1.5 | Geophysical inversion (VES + ERT managers, mesh tools) |
| `numpy` | ≥1.21 | Numerical computation |
| `pandas` | ≥1.3 | Data handling and CSV export |
| `matplotlib` | ≥3.4 | Visualisation |

```bash
pip install pygimli numpy pandas matplotlib
```

---

## Usage

```bash
python ert_ves_inversion.py
```

---

## Outputs

| File | Description |
|---|---|
| `ves_inversion.png` | Left: apparent resistivity curve (observed vs calculated). Right: inverted 1D resistivity–depth model with layer labels |
| `ert_2d_section.png` | Top: apparent resistivity pseudosection. Bottom: inverted 2D resistivity section with geological annotations |
| `ert_coverage.png` | log₁₀ model coverage showing inversion sensitivity across the profile |
| `ert_true_model.png` | True forward model mesh for comparison with inverted result |
| `ves_data.csv` | AB/2 spacings, apparent resistivity values, error estimates |
| `ert_data.csv` | Electrode positions (a,b,m,n), apparent resistivity, error |

---

## Research Context

Resistivity surveys are a standard component of hydrogeological baseline studies — including aquifer characterisation, groundwater potential zone delineation, and EIA baseline investigations for mining and infrastructure projects. This workflow provides a reproducible Python alternative to commercial inversion software (RES2DINV, EarthImager), making resistivity analysis fully accessible in an open-source environment.

Field exposure to Earth Resistivity Meter (VES) measurements forms the practical basis for this computational implementation.

---

## Author

**Banani Jana**
M.Sc. Applied Geology (2nd Year), Presidency University, Kolkata
Email: bananijana2002@gmail.com

---

## License

MIT License
