# simDER: rod in confined space

This repository contains the legacy discrete-elastic-rod simulator that will be
used as the mechanics foundation for the mechanics-informed RRT described in
[`doc/README_mechanics_informed_RRT.md`](doc/README_mechanics_informed_RRT.md).

## Build

The project requires a C++17 compiler, Eigen3, LAPACK, OpenGL, GLU, and GLUT.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Run

Run the OpenGL demonstration from the repository root so the legacy relative
mesh paths resolve correctly:

```bash
./build/simDER option.txt
```

Command-line overrides follow `--`:

```bash
./build/simDER option.txt -- render 0 saveData 0 totalTime 0.05
```

## Test

```bash
ctest --test-dir build --output-on-failure
```

The current tests include a short headless regression check for the legacy
dynamic solver and a zero-load static-equilibrium residual check. Additional
mechanics validation tests for derivatives, magnetic loading, and contact will
be added as those interfaces are extracted.
