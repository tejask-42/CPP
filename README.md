# Numerical Techniques Bootcamp

This repository contains C++ implementations of various numerical techniques commonly used in computational physics. The focus is on solving equations, numerical integration, and differential equations.

## Topics Covered

### Week 1
* [**Newton-Raphson Method**](https://www.geeksforgeeks.org/newton-raphson-method/): Solving for roots of equations.
* [**Secant Method**](https://www.geeksforgeeks.org/secant-method-of-numerical-analysis/): Root-finding without derivatives.
* [**Euler's Method**](https://tutorial.math.lamar.edu/classes/de/eulersmethod.aspx): Solving ordinary differential equations (ODEs).
* [**Runge-Kutta 4th Order Method**](https://math.libretexts.org/Courses/Monroe_Community_College/MTH_225_Differential_Equations/03%3A_Numerical_Methods/3.03%3A_The_Runge-Kutta_Method): High-accuracy ODE solution.
* **Solving Linear Systems**:
    * [**Jacobi Method**](https://en.wikipedia.org/wiki/Jacobi_method)
    * [**Gauss-Seidel Method**](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)

### Week 2
* [**Solving Second-Order ODEs**](https://kyleniemeyer.github.io/ME373-book/content/second-order/numerical-methods.html): Numerical solutions using Euler and RK4 methods, including a comparison with exact solutions.  
* [**Van der Pol Oscillator**](https://en.wikipedia.org/wiki/Van_der_Pol_oscillator): Implementation of the oscillator equations and animation using Euler and RK4 methods.  
* [**Lorenz Attractor**](https://en.wikipedia.org/wiki/Lorenz_system): Plotting chaotic trajectories in 2D using Euler and RK4 methods, showcasing sensitivity to initial conditions.
  
## How to Use
Clone this repository, compile, and run the executables:
```bash
git clone https://github.com/tejask-42/Numerical-Techniques-Bootcamp.git
g++ filename.cpp -o output
./output
