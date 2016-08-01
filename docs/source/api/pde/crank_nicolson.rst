Heat Equation via a Crank-Nicolson scheme
============================

The heat equations in 1D and 2D can be expressed as:

$$ \\frac{ \\partial Q}{\\partial t} = \\frac{ \\partial}{\\partial x} \\left( D \\frac{ \\partial Q}{ \\partial x} \\right) $$

$$ \\frac{ \\partial Q}{\\partial t} = \\frac{ \\partial}{\\partial y} \\left( D \\frac{ \\partial Q}{ \\partial y} \\right)+\\frac{ \\partial}{\\partial y} \\left( D \\frac{ \\partial Q}{ \\partial y} \\right) $$

For 1D and constant coefficient D, using a finite differencing method we can obtain a stable algoritm (unlike for the wave equation):

$$ \\frac{Q^{n+1}_j - Q^{n}_j}{\\Delta t} = D \\left( \\frac{Q^n_{j+1} - 2 Q^n_j + Q^n_{j-1}}{\\Delta x^2} \\right) $$

However in order to observe features of scale \\(\\lambda \\gg \\Delta x\\) with a stable result, the number of time steps needed is unfeasibly large. As usual we will have to be smarter. 

The above method is called fully explicit, if instead we evaulate the RHS at the time step \\(t_{n+1}\\) we create a fully implicit method:

$$ \\frac{Q^{n+1}_j - Q^{n}_j}{\\Delta t} = D \\left( \\frac{Q^{n+1}_{j+1} - 2 Q^{n+1}_j + Q^{n+1}_{j-1}}{\\Delta x^2} \\right) $$

This scheme is unconditionally stable yet first order in time and second order in space. We can form a method which is second order in both space and time and unconditionally stable by forming the average of the explicit and implicit schemes. This is the Crank-Nicolson scheme:

$$ \\frac{Q^{n+1}_j - Q^{n}_j}{\\Delta t} = \\frac{D}{2} \\left( \\frac{Q^{n+1}_{j+1} - 2 Q^{n+1}_j + Q^{n+1}_{j-1} + Q^n_{j+1} - 2 Q^n_j + Q^n_{j-1} }{\\Delta x^2} \\right) $$

We now have a suitable algorithm for solving the heat equation. But it would seem it requires knowledge of \\(Q\\) at later time steps. However this notion can be dispelled by writing the above in a matrix equation form:

$$ \\begin{bmatrix} & \\vdots & \\vdots & \\vdots & \\\\ \\cdots & 1+s & -s/2 & 0 & \\cdots \\\\ \\cdots & -s/2 & 1+s & -s/2 & \\cdots \\\\ \\cdots & 0 & -s/2 & 1+s & \\cdots \\\\ & \\vdots & \\vdots & \\vdots & \\end{bmatrix} \\begin{bmatrix} \\vdots \\\\ Q^{n+1}_{j+1} \\\\ Q^{n+1}_{j} \\\\ Q^{n+1}_{j-1} \\\\ \\vdots \\end{bmatrix} = \\begin{bmatrix} & \\vdots & \\vdots & \\vdots & \\\\ \\cdots & 1-s & s/2 & 0 & \\cdots \\\\ \\cdots & s/2 & 1-s & s/2 & \\cdots \\\\ \\cdots & 0 & s/2 & 1-s & \\cdots \\\\ & \\vdots & \\vdots & \\vdots & \\end{bmatrix} \\begin{bmatrix} \\vdots \\\\ Q^{n}_{j+1} \\\\ Q^{n}_{j} \\\\ Q^{n}_{j-1} \\\\ \\vdots \\end{bmatrix} + \\begin{bmatrix} \\vdots \\\\  \\\\ b.c.s \\\\  \\\\ \\vdots \\end{bmatrix}$$

Where \\(s = \\frac{D \\Delta t}{\\Delta x^2} \\)

Hence the matrix equation \\(Ax = B \\) must be solved where \\(A\\) is a tridiagonal matrix. Here we can use SciPy's solve_banded function to solve the above equation and advance one time step for all the points on the spatial grid.

Solving for the diffusion of a Gaussian we can compare to the analytic solution, the heat kernel:

$$ Q(x,t) \\propto \\frac{1}{\\sqrt{2 \\pi (\\sigma_0^2 + 2Dt)}} \\exp \\left( -\\frac{(x-x_0)^2}{2(\\sigma_0^2 + 2Dt)}  \\right) $$

In the below video, the red outline shows the analytic solution and the black solid line shows the Crank-Nicolson result

.. raw:: html

 <video width="640" height="480" controls>
   <source src="https://raw.githubusercontent.com/PyCav/PyCav-Library/master/docs/source/api/pde/kernel.mp4" type="video/mp4">
 Your browser does not support the video tag.
 </video> 

Variable Coefficient
^^^^^^^^^^^^^^^^^

If \\(D\\) is a function of position then \\(s\\) needs to be evaulated at the spatial point. The form of this finite difference can be seen in finite difference methods page under introductory documentation.

Argument list
^^^^^^^^^^^^^^^^

CN_diffusion_equation(T_0, D, x_list, dx, N_t, s = 0.25, wall_T = [0.0,0.0,0.0,0.0])

   This function performs the Crank-Nicolson scheme for 1D and 2D problems to solve the inital value problem for the heat equation.

   **Parameters:**

   *T_0: numpy array*

   In 1D, an N element numpy array containing the intial values of T at the spatial grid points. In 2D, a NxM array is needed where N is the number of x grid points, M the number of y grid points. This array needs to be in "matrix indexing" rather than "Cartesian indexing" i.e. the first index (the rows) correspond to x values and the second index (the columns) correspond to y values. If using numpy.meshgrid, matrix indexing can be ensured by using the indexing='ij' keyword arg.
   
   *D: function*
   
   In 1D, must take a numpy array argument containing spatial coords and return a numpy array of equal length giving the value of the diffusivity at the given positions e.g.

   .. code-block:: python
   
    def D(x):
      return 0.5+0.5*x

   In 2D, must take a pair of floats of the x and y coords and return a float of the diffusivity at that point e.g.

   .. code-block:: python
   
    def D(x,y):
      return 0.5+0.5*(x-0.5)**2+0.5*(y-0.5)**2

   *x_list: numpy array / list of numpy array*

   In 1D, an N element numpy array of equally spaced points in space (creating using numpy linspace or arange is advised) at which the wave will be evaluated. In 2D, a list containing two numpy arrays of length N and M respectively. These correspond to the x and y spatial grids. e.g.

   .. code-block:: python
   
    dx = 0.01
    x = dx*np.arange(201)
    y = dx*np.arange(101)
    T,t = pde.CN_diffusion_equation(T_0, D, [x,y], dx, N_t)

   *dx: float*

   Must give the spacing between points in the x array (and y array for 2D)
   
   *N_t: integer*
   
   Number of time steps taken

   *s: float*
   
   This is used to set the time step via \\(\\Delta t = s \\Delta x**2 \\). Although the scheme is stable for any size \\(\\Delta t \\) in order to ensure accurate results s should set sufficiently low. Generally of order \\(1/D_{max} \\) is advisable.
 
   *wall_T: list of floats*

   A list of 2 or 4 floats (for 1D or 2D) containing the fixed T values for the boundaries.

   **Returns:**

   A N x N_t numpy array, N x M x N_t in 2D, which contains the approximated T at different times. A N_t element numpy array is also returned containing the time interval over which the simulation was run.