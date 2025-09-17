# cpp_progs
Showing cpp programming skill
1. oop_example/oop_example1.cc
Sample of Object Oriented Programming showing:
i) Class construction
ii) Inheritance
iii) Polymorphism
iv) Abstraction
v) Encapsulation

2. oop_example/nbody
a) Solving numerically Gravitational N-Body problem (system of coupled differential equations of 2nd order)
   using method Runge-Kutta 4.
b) Generate initial conditions of position and velocity that produce analytical solutions
   (e.g. trajectories that move along conic sections (for example, elliptical trajectory) for
   different cases:
   Case 1) Three bodies, with arbitrary masses, on the vertices of equilateral triangle.
   Case 2) Three bodies, with arbitrary masses, aligned, configuration 321.
   Case 3) Three bodies, with arbitrary masses, aligned, configuration 231.
   Case 4) Three bodies, with arbitrary masses, aligned, configuration 213.
   Case 5) Four bodies on the vertices of an isosceles trapezoid, configuration nu = 2.
   Case 6) Four bodies on the vertices of an isosceles trapezoid, configuration nu = 1.
   Case 7) N bodies on vertices of a regular polygon (N = 6) with equal mass.

Note.- The program generates nbody_output.json file that contains trajectories of each of the bodies
       with respect to center of mass of the system.
       Plotting those trajectories you will see that each of the bodies describe an elliptical
       trajectory.

3. codewars_progs
Solving eight problems from codewars website

4. optimization_with_mkl
Implementation of standard geometry optimization algorithm: BFGS using MKL
(among other algorithms).

5. optimization_with_CUDA
Implementation of standard geometry optimization algorithm: BFGS using cuBLAS, cuSolver
(among other algorithms).

