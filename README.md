# N-Body Gravity Simulator

A two-dimensional Newtonian gravity simulation written in C++. The project compares numerical integration methods and direct force calculation with the Barnes-Hut algorithm.

## Current features

- Direct O(n²) gravitational force calculation
- Barnes-Hut quadtree force approximation
- Direct-versus-Barnes-Hut acceleration error measurement
- Euler and Velocity Verlet integration
- Energy and energy-drift measurements
- Random systems containing 1–100 bodies
- Three-body figure-eight example
- Browser-based trajectory visualisation

## Build and run

```bash
make
./nbody
```

The program asks for the number of bodies and writes their positions to `trajectory.csv`.

To run the figure-eight example:

```bash
./nbody figure-eight
```

## Visualisation

After running the simulator, start a local server:

```bash
python3 -m http.server 8000
```

Then open [http://localhost:8000/viewer.html](http://localhost:8000/viewer.html) in a browser.

## Roadmap

- Benchmark direct and Barnes-Hut force calculations
- Compare numerical accuracy and energy conservation
- Add automated tests
