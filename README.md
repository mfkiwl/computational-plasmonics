This repository is dedicated to the numerical computation of surface plasmons using both finite element and boundary element methods.

Surface plasmons are computed by solving the **plasmonic eigenvalue problems (PEP)**.

The following cases are implemented:
- Smooth particle.
- Camembert in Cartesian and Euler coordinates.
- Poké-ball.
- Smooth particle perturbed by a straight corner. 


# Description

	fem/

Discretization of the PEP using finite element methods. The discretization is implemented in Python and relies on the following packages:
- `gmsh`: mesh generation.
- `fenics` or `fenicsx`: finite element assembly.
- `petsc4py`: sparse linear algebra.
- `slepc`: sparse eigensolver.

Most scripts rely on modules available in the repository [python-scientific-computing-utils](https://github.com/fmonteghetti/python-scientific-computing-utils). These modules are identified by the `_utils` suffix.


	bem/

Discretization of the PEP using high-order Nyström boundary element methods. Implementation in Julia.