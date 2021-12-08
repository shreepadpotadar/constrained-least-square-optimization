# Constrained Least Square Optimization
This code implements the numerical algorithm to solve the simplest random constrained least-square optimization.


### The following software packages must be installed to run the code on Linux/MacOS
```
Eigen C++ open source library
C++ compiler
Cmake
```


### To create a C++ executable program run the following commands on the terminal:
```
cd main_folder
make
```
This creates C++ executable file `constrained-lso` inside the folder `build`


### To use this program the run the following command, 
`build/constrained-lso`
This will prompt a message asking the user to input M, N and Sigma values




### The following papers are used to implement this code repo.
[1] **Gene H. Golub and Urs von Matt**, *Quadratically constrained least squares and quadratic problems*. Numer. Math. 59, 561-580 (1991).

[2] **Yan V. Fyodorov and Rashel Tublin**, *Counting stationary points of the loss function in the simplest constrained least-square optimization*. arXiv:1911.12452v2 [math.PR], 1 Oct 2020.

[3] **W. Gander**. *Least Squares with a Quadratic Constraint*. Numer. Math. 36, 291-307 (1981).

[4] **Christopher C. Paige and Michael A. Saunders**, *LSQR: An Algorithm for Sparse Linear Equations
and Sparse Least Squares*. ACM Trans. Math. Softw. 8, 2, pages 43-71 (June 1982).

