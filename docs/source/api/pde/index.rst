Partial Differential Equations
====================================

Introduction to the Partial Differential Equation module
--------------------------

This module contains functions for use in solving PDEs, namely the wave equation, the heat equation and the time dependent Schrödinger equation. These are solved using the Lax-Wendroff, Crank-Nicolson and Split Step Fourier methods respectively. Documentation explaining the use of each function and the algorithms used within will be presented in individual sections. A short introduction to finite difference methods will also be presented.

Introductory Documentation
-------

.. toctree::
   :maxdepth: 1

   finitedifference

In-depth Documentation
-------

.. toctree::
   :maxdepth: 1
   
   lax_wendroff
   crank_nicolson
   split_step


The list below summarises the functions, their input arguments and their outputs for quick reference for the informed user:

Functions
---------

LW_wave_equation(psi_0, x_list, dx, N_t, c, a = 1., bound_cond = 'periodic',init_grad = None, init_vel = None)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   This function performs the two-step Lax-Wendroff scheme for 1D problems and a Lax method for 2D problems to solve a flux-conservative form of the wave equation for variable wave speed, c. 

   **Parameters:**

   *psi_0: numpy array*

   In 1D, an N element numpy array containing the intial values of \\(\\psi\\) at the spatial grid points. In 2D, a NxM array is needed where N is the number of x grid points, M the number of y grid points. This array needs to be in "matrix indexing" rather than "Cartesian indexing" i.e. the first index (the rows) correspond to x values and the second index (the columns) correspond to y values. If using numpy.meshgrid, matrix indexing can be ensured by using the indexing='ij' keyword arg.

   *x_list: numpy array / list of numpy array*

   In 1D, an N element numpy array of equally spaced points in space (creating using numpy linspace or arange is advised) at which the wave will be evaluated. In 2D, a list containing two numpy arrays of length N and M respectively. These correspond to the x and y spatial grids. e.g.

   .. code-block:: python
   
    dx = 0.01
    x = dx*np.arange(201)
    y = dx*np.arange(101)
    psi_2d,t = pde.LW_wave_equation(psi_0_2d,[x,y],dx,N,c_2d)

   *dx: float*

   Must give the spacing between points in the x array (and y array for 2D)
   
   *N_t: integer*
   
   Number of time steps taken
   
   *c: function*
   
   In 1D, must take a numpy array argument containing spatial coords and return a numpy array of equal length giving the value of the wave speed at the given positions e.g.

   .. code-block:: python
   
    def c(x):
      return 0.5+0.5*x

   In 2D, must take a pair of numpy arrays containing the x and y coords and return a numpy meshgrid of the wave speeds at those points e.g.

   .. code-block:: python
   
    def c(x,y):
      XX,YY = np.meshgrid(x,y,indexing='ij')
      return 0.5+0.5*YY
   
   This gives a wavespeed that's only a function of y

   *a: float*
   
   The Courant number, for stability of the code this must be \\(\\leq 1\\) (look up Courant-Friedrichs-Lewy stability criterion for information on this). For lower a, the code is more stable but the time step is reduced so more time steps (N) are required to simulate the same time length 
   
   *bound_cond: string*
   
   Can be equal to 'fixed', 'reflective' and 'periodic' to impose those boundary conditions. For fixed, the wave must go to zero at the boundary. For reflective, the gradient parallel to the surface normal must vanish at the boundary. For periodic, the boundaries on opposite sides are set to be equal.

   *init_grad: function*

   A function which takes psi_0 as an argument and returns the gradient of the initial wave on the spatial grid. 1D example for a travelling Gaussian given below along with the init_vel example. For 2D, both \\(\\partial \\psi / \\partial x \\) and \\(\\partial \\psi / \\partial y \\) must be returned individually. For a 2D initially Gaussian wave:

   $$ \\psi_0 (x,y) = \\exp (- ((x - \\mu_x )^2+(y - \\mu_y )^2) / 2 \\sigma^2 ) \\to \\frac{ \\partial \\psi }{ \\partial x} = -(x- \\mu_x) \\psi_0 / \\sigma^2 $$

   .. code-block:: python

     def twoD_gaussian(XX,YY,mean,std):
      return np.exp(-((XX-mean[0])**2+(YY-mean[1])**2)/(2*std**2))

    def gradient_2d(x,y,mean,std):
      XX,YY = np.meshgrid(x,y, indexing='ij')
      def D(psi_0):
         dfdx = -(XX-mean[0])*twoD_gaussian(XX,YY,mean,std)/std**2
         dfdy = -(YY-mean[1])*twoD_gaussian(XX,YY,mean,std)/std**2
         return dfdx,dfdy
      return gradient_2d

   Here the init_grad argument would be set to gradient_2d(x,y,mean,std) so that the LW_wave_equation program recieves the function D. This removes the need for LW_wave_equation to know the values of mean and std. 

   If the default argument, None, is given then the initial gradient is estimated within the program using finite differencing. It is preferable to give the program a init_grad function when there exists an analytic form.

   *init_vel: function*

   A function which takes psi_0 as an argument and returns the velocity (\\(\\partial \\psi / \\partial t \\)) of the initial wave on the spatial grid. 1D example for a travelling Gaussian given below.

   If the default argument, None, is given then the initial velocity is set to zero at all points.

   Having defined the variables; x, dx, N_t, mean and std:

   .. code-block:: python
   
    def oneD_gaussian(x,mean,std):
      return np.exp(-((x-mean)**2)/(2*std**2))

    def gradient_1d(x,mean,std):
      def D(psi_0):
         return -(x-mean)*oneD_gaussian(x,mean,std)/std**2
      return D

    def velocity_1d(x,mean,std):
      def V(psi_0):
         return -c(x)*(x-mean)*oneD_gaussian(x,mean,std)/std**2
      return V

    psi_1d,t = pde.LW_wave_equation(oneD_gaussian(x,mean,std),x,dx,N_t,c, 
            init_vel = velocity_1d(x,mean,std), init_grad = gradient_1d(x,mean,std),
            bound_cond = 'reflective')
 
   **Returns:**

   A N x N_t numpy array, N x M x N_t in 2D, which contains the approximated wave at different times. A N_t element numpy array is also returned containing the time interval over which the simulation was run.

CN_diffusion_equation(T_0, D, x_list, dx, N_list, s = 0.25, wall_T = [0.0,0.0,0.0,0.0])
^^^^^^^^^^^^^^^^^^^^^^^^^^

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


split_step_schrodinger(psi_0, dx, dt, V, N, x_0 = 0., k_0 = None, m = 1.0, non_linear = False)
^^^^^^^^^^^^^^^^^^^^^^^^^^

   This function performs the split-step Fourier method to solve the 1D time-dependent Schrödinger equation for a given potential

   **Parameters:**

   *psi_0: numpy array*

   In 1D, an N element numpy array containing the intial values of \\(\\psi\\) at the spatial grid points. In 2D, a NxM array is needed where N is the number of x grid points, M the number of y grid points. This array needs to be in "matrix indexing" rather than "Cartesian indexing" i.e. the first index (the rows) correspond to x values and the second index (the columns) correspond to y values. If using numpy.meshgrid, matrix indexing can be ensured by using the indexing='ij' keyword arg.

   *dx: float*

   Must give the spacing between points in the x array

   *dt: float*

   Gives the time step taken within the split-step algorithm. This needs to be small to reduce the size of numerical errors (try 0.01 as a safe starting value)

   *V: function*

   Pass a function which takes a numpy array argument containing spatial coords and returns the potential at that point e.g.

   .. code-block:: python

    def V(x):
      V_x = np.zeros_like(x)
      a = 0.5
      x_mid = (x.max()+x.min())/2.
      V_x = -a**2*(1/np.cosh(a*(x-x_mid)))**2
      return V_x

   If non_linear = True then the potential function must now take an additional argument which is equal to the spatial wavefunction at the current time step e.g.

   .. code-block:: python

    def V(x,psi):
      V_x = np.zeros_like(x)
      V_x = -200.*np.absolute(psi)**2+0.05*x**2
      return V_x

   *N_t: integer*
   
   Number of time steps taken

   *x_0: float*

   Give the starting position of the spatial grid

   *k_0: float*

   Gives the starting position of the momentum space grid. If none is given then k_0 is set to \\(-\\pi/ \\Delta x \\) as it can be shown that this exactly satisfies the Nyquist limit.

   *m: float*

   The mass of the particle (default value of 1.0)

   *non_linear: boolean*

   Set to True if investigating the non-linear Schrödinger equation. Default is False
   
   **Returns:**

   Two N x N_t numpy arrays which contain the approximated real space and momentum space wavefunctions at different times. A N element numpy array is also returned containing the k space interval used.
