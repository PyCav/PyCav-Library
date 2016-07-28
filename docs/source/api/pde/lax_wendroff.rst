Wave Equation via Lax/Lax-Wendroff schemes
===================

The wave equations in 1D and 2D can be expressed as (for constant wave speed):

$$ \\frac{ \\partial^2 \\psi}{\\partial t^2} = c^2 \\frac{ \\partial^2 \\psi}{\\partial x^2} $$

$$ \\frac{ \\partial^2 \\psi}{\\partial t^2} = c^2 ( \\frac{ \\partial^2 \\psi}{\\partial x^2} + \\frac{ \\partial^2 \\psi}{\\partial y^2} )$$

If we re-express the 1D wave equation in a flux-conservative form, \\(\\frac{\\partial u}{\\partial t} = -\\frac{\\partial F}{\\partial x} \\) (which allows for the use of established numerical methods) then we obtain:

$$ u = \begin{array} r \\\ s \end{array}, r \equiv c \\frac{\\partial \\psi}{\\partial x}, s \equiv \\frac{\\partial \\psi}{\\partial t} $$

LW_wave_equation(psi_0, x_list, dx, N_t, c, a = 1., bound_cond = 'periodic',init_grad = None, init_vel = None)

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
   
   *a: float*
   
   The Courant number, for stability of the code this must be \\(\\leq 1\\) (look up Courant-Friedrichs-Lewy stability criterion for information on this). For lower a, the code is more stable but the time step is reduced so more time steps (N) are required to simulate the same time length 
   
   *bound_cond: string*
   
   Can be equal to 'fixed', 'reflective' and 'periodic' to impose those boundary conditions. For fixed, the wave must go to zero at the boundary. For reflective, the gradient parallel to the surface normal must vanish at the boundary. For periodic, the boundaries on opposite sides are set to be equal.

   *init_grad: function*

   A function which takes psi_0 as an argument and returns the gradient of the initial wave on the spatial grid. 1D example for a travelling Gaussian given below along with the init_vel example. For 2D, both \\(\\partial \\psi / \\partial x \\) and \\(\\partial \\psi / \\partial y \\) must be returned individually. For a 2D initially Gaussian wave:

   $$ \\psi_0 = \\exp (- ((x - \\mu_x )^2+(y - \\mu_y )^2) / 2 \\sigma^2 ) \\to \\frac{ \\partial \\psi }{ \\partial x} = -(x- \\mu_x) \\psi_0 / \\sigma^2 $$

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