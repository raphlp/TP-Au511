# Control of Aircraft - Practical Work (Au511)

## 📋 Project Overview

This practical work implements a **longitudinal autopilot** for a MIRAGE III class fighter aircraft. The project covers the complete design of a cascaded control system with three feedback loops (q, γ, z) and includes flight management simulation.

**Subject Number:** 82  
**Operating Point:** Altitude = 21,230 ft | Mach = 0.95

## 👥 Authors

- LAUPIES Raphaël
- CHARDON DU RANQUET Quentin

---

## 🚀 Quick Start

### Prerequisites

- Python 3.8+
- Required packages:

```bash
pip install numpy matplotlib scipy control
```

> **Note:** `slycot` is optional (used for interactive sisotool). The script works without it since all gains are pre-calculated.

### Running the Script

```bash
cd "/path/to/TP Au511"
python3 "TP 1.py"
```

### Expected Behavior

1. **Console output** with numerical results (equilibrium point, poles, damping ratios, settling times)
2. **Multiple plot windows** will open sequentially - **close each window to continue execution**:
   - Step response α/δm and q/δm
   - Washout filter comparison
   - Gamma loop step response
   - Z loop step response
   - Flight management simulation (4 phases)
   - Alternative flare comparison (ramp vs exponential)

---

## 📁 Project Structure

```
TP Au511/
├── TP 1.py                              # Main script (all sections)
├── atm_std.py                           # US Standard Atmosphere 1976 model
├── sisopy31.py                          # SISO tool for interactive gain tuning
├── ussa76ut86.txt                       # Atmospheric data table
├── PracticalWork_Control_of_aircraft.pdf # Assignment document
├── stepalphaq.pdf                       # Generated plot (step response)
└── README.md                            # This file
```

---

## 📐 Operating Point (Subject 82)

| Parameter | Value |
|-----------|-------|
| Altitude | 21,230 ft (6,471 m) |
| Mach number | 0.95 |
| Air density (ρ) | 0.628 kg/m³ |
| Speed of sound (a) | 314.5 m/s |
| True airspeed (V_eq) | 298.8 m/s |

### Aerodynamic Coefficients (M = 0.95)

| Coefficient | Value | Description |
|-------------|-------|-------------|
| Cx0 | 0.018 | Zero-lift drag |
| k | 0.28 | Polar coefficient |
| Cz_α | 2.85 rad⁻¹ | Lift gradient w.r.t. incidence |
| Cz_δm | 0.75 rad⁻¹ | Lift gradient w.r.t. elevator |
| δm0 | 0.002 rad | Equilibrium elevator for null lift |
| α0 | 0.013 rad | Incidence for null lift |
| f | 0.545 | Aerodynamic center (body+wings) |
| f_δ | 0.82 | Aerodynamic center (fins) |
| Cm_q | -0.55 s/rad | Damping coefficient |

---

## 🔧 Implementation Summary

### 1. Equilibrium Point Calculation
Iterative algorithm to find α_eq, δm_eq, and Fpx_eq at the operating point.

**Results:**
- α_eq = 2.59°
- δm_eq = -0.48°
- Fpx_eq = 19,141 N

### 2. State Space Model
6-state longitudinal model: X = [V, γ, α, q, θ, z]ᵀ

### 3. Open-Loop Mode Analysis
- **Short period:** ωn = 3.16 rad/s, ξ = 0.30
- **Phugoid:** ωn = 0.044 rad/s, ξ = 0.17

### 4. Controller Synthesis

| Loop | Gain | Damping Ratio | Specifications |
|------|------|---------------|----------------|
| q (pitch rate) | Kr = -0.10 | ξ ≈ 0.65 | Per design requirement |
| γ (flight path) | Kγ = 5.0 | ξ ≥ 0.5 | OS ≤ 15%, Ts minimized |
| z (altitude) | Kz = 0.001 | ξ ≥ 0.5 | OS ≈ 2.5%, Ts ≈ 5s |

### 5. Saturation Analysis
Maximum flight path angle γ_max calculated using bisection method to ensure load factor Δnz ≤ 3g.

### 6. Flight Management
Four-phase simulation:
1. **Ascent** - constant γ (gamma hold)
2. **Cruise** - constant z (altitude hold, ~100s)
3. **Descent** - constant γ (gamma hold)
4. **Flare** - altitude hold for landing

### 7. Alternative Flare Commands
- **Ramp command** - linear altitude decrease
- **Exponential command** - smooth exponential decay (τ = 10s)

---

## 📊 Generated Outputs

The script generates several PDF plots:
- `stepalphaq.pdf` - Step response of α and q

Console outputs include:
- Equilibrium point values
- Transfer functions
- Pole locations with damping ratios
- Settling times and overshoots
- γ_max saturation value

---

## ⚠️ Troubleshooting

### "No module named 'slycot'"
This is optional. The sisotool calls are commented out, and all gains are pre-calculated.

### Warnings about "invalid value" or "divide by zero"
These warnings are normal for some edge cases in the control library and don't affect results.

### Plots don't appear
Make sure to close each plot window to continue script execution.

---

## 📚 References

- Course: Au511 - Control of Aircraft
- Instructor: Jean-Pierre NOUAILLE
- Date: November 2025
- Aircraft model: MIRAGE III class fighter

---

## 📝 License

Academic project - ISAE-SUPAERO / ENSICA
