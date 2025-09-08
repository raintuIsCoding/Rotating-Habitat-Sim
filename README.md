# Rotating Space Habitat Simulator (MATLAB)

A real-time 3D simulator for artificial-gravity space stations. Adjust ring radius and RPM, visualize a torus habitat, and get human-factors readouts (g at head/feet for a 2 m person, Δg across the body, Coriolis cues). Includes one-click solvers to hit 1 g by solving for either radius or RPM.

![Demo](docs/demo.gif)

---

## Features
- **Interactive 3D station**: torus geometry with smooth camera + HUD.
- **Human-factors readouts**: head/feet g for a 2 m person, Δg across body, Coriolis at walking speed.
- **One-click 1-g solvers**: compute RPM for a target g at a given radius, or compute radius for a target g at a given RPM.
- **Design trade helpers**: quick sanity checks for artificial-gravity constraints before deeper analysis.

## Quick Start
1. Open the repo in MATLAB (R2023b+ recommended).
2. Run the main script/function:
   ```matlab
   >> spinStation_minimal
