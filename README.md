# FEM 2D Heat Conduction Solver

MATLAB implementation of a 2D steady-state heat conduction problem using the Finite Element Method (FEM). The solver uses linear triangular elements to compute temperature distribution across a square steel plate.

## Features
- Mesh generation and element connectivity
- Global stiffness matrix assembly
- Application of Dirichlet boundary conditions
- Visualization of temperature distribution
- Clean and modular MATLAB code

## Files
- `FEM_1a.m`: When three sides of the plate are insulated and the 4th side is maintained at 25°C. Solved using 2 elements and 𝑞=1000*𝛿(𝑥−𝑥𝑐,𝑦−𝑦𝑐 ) W/m^2, where xc and yc are the coordinates of the center of the plate, and 𝛿 is the Dirac-delta function.
- `FEM_1b.m`: When three sides of the plate are insulated and the 4th side is maintained at 25°C. Solved using 4 elements intersecting at the center.
- `FEM_1c.m`: When three sides of the plate are insulated and the 4th side is maintained at 25°C. Solved using ~10 elements with heat source at the intersection of elements.
- `FEM_1d.m`: When three sides of the plate are insulated and the 4th side is maintained at 25°C. Solved using ~100 elements with heat source at the intersection of elements.
- `FEM_2.m`: When three sides of the plate are insulated and the 4th side is subjected to convective heat transfer such that 𝑞𝑛=ℎ(𝑇−𝑇0) where T0 = 25°C and h = 10 W/(m^2-K). Solved using ~100 elements with heat source at the intersection of elements.

## Author
Dipan Laha
