Wave Equation via Lax/Lax-Wendroff schemes
===================

The wave equations in 1D and 2D can be expressed as (for constant wave speed):

$$ \\frac{ \\partial^2 \\psi}{\\partial t^2} = c^2 \\frac{ \\partial^2 \\psi}{\\partial x^2} $$

$$ \\frac{ \\partial^2 \\psi}{\\partial t^2} = c^2 \\left( \\frac{ \\partial^2 \\psi}{\\partial x^2} + \\frac{ \\partial^2 \\psi}{\\partial y^2} \\right) $$

A large class of inital value PDEs can be case into a flux-conservative form. In 1D:

$$ \\frac{\\partial u}{\\partial t} = -\\frac{\\partial F}{\\partial x} $$

If we re-express the 1D wave equation in a flux-conservative form (which allows for the use of established numerical methods) then we obtain:

$$ u = \\begin{bmatrix}r \\\\ s\\end{bmatrix}, r \\equiv c \\frac{\\partial \\psi}{\\partial x}, s \\equiv \\frac{\\partial \\psi}{\\partial t} $$

$$ F = \\begin{bmatrix}0 & -c \\\\ -c & 0\\end{bmatrix}\\cdot \\begin{bmatrix}r \\\\ s\\end{bmatrix} $$

We now have a first order differential equation to solve of a given inital value problem.

Unfortunately a simple finite difference method is unconditionally unstable. What a shame. 

For the 1D wave equation we shall use the two-step Lax Wendroff scheme. This includes evaulating u at half steps, using this to find the half step fluxes and using a properly centered expression to perform the full time step. Using the notation introduced in the finite difference method documentation page:

$$ u^{n+1/2}_{j+1/2} = \\frac{1}{2} (u_{j+1}^n+u_j^n) - \\frac{\\Delta t}{2 \\Delta x}(F^n_{j+1} - F^n_j) $$

$$ u^{n+1}_{j} = u_j^n - \\frac{\\Delta t}{\\Delta x}(F^{n+1/2}_{j+1/2} - F^{n+1/2}_{j-1/2}) $$ 

The stability of the scheme is parameterised by the Courant number, \\(\\alpha = c \\Delta t / \\Delta x \\), and the criteria for stability is that \\(\\alpha \\leq 1 \\). This ensures \\(c \\Delta t \\leq \\Delta x \\), hence information between spatial points can not have been communicated before the next time step is calculated.

Boundary conditions can be set to be periodic, fixed or reflective. These are defined in our flux-conservative formulation as:

Fixed (noting \\(s = \\frac{\\partial \\psi}{\\partial t} \\)):

$$ s(x = x_{min}) = s(x = x_{max}) = 0 $$

Reflective (noting \\(r = c \\frac{\\partial \\psi}{\\partial x} \\)):

$$ r(x = x_{min}) = r(x = x_{max}) = 0 $$

Periodic:

$$ \\psi_N = \\psi_1, \\psi_0 = \\psi_{N-1} $$

Now if we take the wave speed to be a function of position, \\( c(x) \\), then when the flux, \\(F\\), is evaulated at grid points then the wave speed at that point must be used e.g

$$ F^n_{j+1} = - \\begin{bmatrix}c(x_{j+1}) s^n_{j+1} \\\\ c(x_{j+1}) r^n_{j+1} \\end{bmatrix} $$

As an example of its use, here is a gaussian disturbance reflecting from a wall with fixed boundary conditions:

.. raw:: html

 <video width="640" height="480" controls>
   <source src="https://raw.githubusercontent.com/PyCav/PyCav-Library/master/docs/source/api/pde/wibble.mp4" type="video/mp4">
 Your browser does not support the video tag.
 </video> 

2D Lax Scheme
^^^^^^^^^^^^^^^

In 2D, the flux conservative form includes an additional flux:

$$ \\frac{\\partial u}{\\partial t} = -\\frac{\\partial F_x}{\\partial x}-\\frac{\\partial F_y}{\\partial y} $$

We define our new \\(u\\) vector as the following:

$$ u = \\begin{bmatrix}r \\\\ l \\\\ s\\end{bmatrix}, r \\equiv c \\frac{\\partial \\psi}{\\partial x}, l \\equiv c \\frac{\\partial \\psi}{\\partial y}, s \\equiv \\frac{\\partial \\psi}{\\partial t} $$

$$ \\frac{\\partial r}{\\partial t} = \\frac{\\partial}{\\partial x} (c s)$$

$$ \\frac{\\partial l}{\\partial t} = \\frac{\\partial}{\\partial y} (c s)$$

$$ \\frac{\\partial s}{\\partial t} = \\frac{\\partial}{\\partial x}(c r) + \\frac{\\partial}{\\partial y} (c l)$$

Hence the fluxes are given by:

$$ F_x = \\begin{bmatrix}0 & 0 & -c \\\\ 0 & 0 & 0 \\\\ -c & 0 & 0\\end{bmatrix}\\cdot \\begin{bmatrix}r \\\\ l \\\\ s\\end{bmatrix} $$

$$ F_y = \\begin{bmatrix}0 & 0 & 0 \\\\ 0 & 0 & -c \\\\ 0 & -c & 0\\end{bmatrix}\\cdot \\begin{bmatrix}r \\\\ l \\\\ s\\end{bmatrix} $$

Plugging in the definitions of \\(r\\), \\(l\\) and \\(s\\) to the expression for \\(\\dot{s}\\), the 2D wave equation is recovered.

Using a 2D Lax scheme, the components of \\(u\\) can be calculated by:

$$ u^{n+1}_{j,l} = \\frac{1}{4} (u^n_{j+1,l} + u^n_{j-1,l} + u^n_{j,l+1} + u^n_{j,l-1}) - \\frac{\\Delta t}{2 \\Delta}(F^n_{x,j+1,l}-F^n_{x,j-1,l}+F^n_{y,j,l+1}-F^n_{y,j,l-1}) $$

Where \\(\\Delta = \\Delta x = \\Delta y \\). The fluxs labelled with \\(x,y\\) refer to terms within the \\(x,y\\) partial derivatives in the flux conservative expressions above. e.g. \\(F_x\\) for \\(r\\) is \\(c s\\) while \\(F_y\\) is zero.

Boundaries conditions are dealt with in a similar way to the 1D case. But now \\(l\\) vanishes at \\(y_{min}\\) and \\(y_{max}\\) for reflective boundary conditions.

For 2D cases, the Courant number must now be less than \\( 2^{-1/2} \\). It should be noted that the definition of \\(c\\) in the Courant condition is replaced with the maximal wave speed when the wave speed is allowed to vary with position.

As an example of its use, here is a gaussian disturbance within a container with periodic boundary conditions:

.. raw:: html

 <video width="640" height="480" controls>
   <source src="https://raw.githubusercontent.com/PyCav/PyCav-Library/master/docs/source/api/pde/wibble2.mp4" type="video/mp4">
 Your browser does not support the video tag.
 </video> 

Form of the wave equation for spatially varying wave speed
^^^^^^^^^^^^^^^^^^

A important distinction should be made about the form of the wave equation. There are two possible forms of the wave equation for a variable wave speed, in 1D these are:

$$ \\frac{ \\partial^2 \\psi}{\\partial t^2} = c(x)^2 \\frac{ \\partial^2 \\psi}{\\partial x^2} $$

$$ \\frac{ \\partial^2 \\psi}{\\partial t^2} = \\frac{ \\partial }{\\partial x} \\left( c(x)^2 \\frac{\\partial \\psi}{\\partial x} \\right) $$

In the form we have cast the wave equation we are solving for the second of these equations. This describes systems such as surface waves on a fluid. The first equation follows from the electro-magnetic Maxwell equations in 1D.

*It should be noted when giving positional depedent wavespeeds with discontinuities, this will not give the familiar reflection and transmission results as the additional boundary conditions at the discontinuties are not included*

Argument list
^^^^^^^^^^^^

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