Heat Equation via a Crank-Nicolson scheme
============================

The heat equations in 1D and 2D can be expressed as:

$$ \\frac{ \\partial Q}{\\partial t} = \\frac{ \\partial}{\\partial x} \\left( D \\frac{ \\partial Q}{ \\partial x} \\right) $$

$$ \\frac{ \\partial Q}{\\partial t} = \\frac{ \\partial}{\\partial y} \\left( D \\frac{ \\partial Q}{ \\partial y} \\right)+\\frac{ \\partial}{\\partial y} \\left( D \\frac{ \\partial Q}{ \\partial y} \\right) $$

For 1D and constant coefficient D, using a finite differencing method we can obtain a stable algoritm (unlike for the wave equation):

$$ \\frac{Q^{n+1}_j - Q^{n}_j}{\\Delta t} = D \\left( \\frac{Q^n_{j+1} - 2 Q^n_j + Q^n_{j-1}}{\\Delta x^2} \\right) $$

However in order to observe features of scale \\(\\lambda \gg \\Delta x\\) with a stable result, the number of time steps needed is unfeasibly large. As usual we will have to be smarter. 

The above method is called fully explicit, if instead we evaulate the RHS at the time step \\(t_{n+1}\\) we create a fully implicit method:

$$ \\frac{Q^{n+1}_j - Q^{n}_j}{\\Delta t} = D \\left( \\frac{Q^{n+1}_{j+1} - 2 Q^{n+1}_j + Q^{n+1}_{j-1}}{\\Delta x^2} \\right) $$

This scheme is unconditionally stable yet first order in time and second order in space. We can form a method which is second order in both space and time and unconditionally stable by forming the average of the explicit and implicit schemes. This is the Crank-Nicolson scheme:

$$ \\frac{Q^{n+1}_j - Q^{n}_j}{\\Delta t} = \\frac{D}{2} \\left( \\frac{Q^{n+1}_{j+1} - 2 Q^{n+1}_j + Q^{n+1}_{j-1} + Q^n_{j+1} - 2 Q^n_j + Q^n_{j-1} }{\\Delta x^2} \\right) $$

We now have a suitable algorithm for solving the heat equation. But it would seem it requires knowledge of \\(Q\\) at later time steps. However this notion can be dispelled by writing the above in a matrix equation form:

$$ \\begin{bmatrix} & \\vdots & \\vdots & \\vdots & \\\\ \\cdots & 1+s & -s/2 & 0 & \\cdots \\\\ \\cdots & -s/2 & 1+s & -s/2 & \\cdots \\\\ \\cdots & 0 & 1+s & -s/2 & \\cdots \\\\ & \\vdots & \\vdots & \\vdots & \\end{bmatrix} \\begin{bmatrix} \\vdots \\\\ Q^{n+1}_{j+1} \\\\ Q^{n+1}_{j} \\\\ Q^{n+1}_{j-1} \\\\ \\vdots \\end{bmatrix} = \\begin{bmatrix} & \\vdots & \\vdots & \\vdots & \\\\ \\cdots & 1-s & s/2 & 0 & \\cdots \\\\ \\cdots & s/2 & 1-s & s/2 & \\cdots \\\\ \\cdots & 0 & 1-s & s/2 & \\cdots \\\\ & \\vdots & \\vdots & \\vdots & \\end{bmatrix} \\begin{bmatrix} \\vdots \\\\ Q^{n}_{j+1} \\\\ Q^{n}_{j} \\\\ Q^{n}_{j-1} \\\\ \\vdots \\end{bmatrix} + \\begin{bmatrix}  \\\\  \\\\ b.c.s \\\\  \\\\ \\end{bmatrix}$$