Time-Dependent Schrödinger equation via the Split-Step Fourier method
==============================

Writing the Schrödinger equation in the form (in units where \\(\\hbar = 1 \\):
$$ \\frac{ \\partial \\psi}{\\partial t} = i \\mathcal{L} \\psi+i \\mathcal{N} \\psi $$

Where for the TDSE:
$$ \\mathcal{L} = \\frac{1}{2m} \\frac{\\partial^2}{\\partial x^2},  \\mathcal{N} = -V(x) $$

The time evolution operator is then given by:

$$ \\psi (x,t_0 + t) = \\exp (i (\\mathcal{L} + \\mathcal{N}) t) \\psi (x,t_0) $$

We can split the exponential if \\([ \\mathcal{L}, \\mathcal{N} ] = 0 \\). This is not necessarily true but for a small time interval, \\( \\Delta t \\), commutativity can be assumed with an error of order \\( \\Delta t^2 \\).

By first neglecting \\(\\mathcal{L} \\) in time interval \\([t_0,t_0+ \\Delta t/2]\\) we are left with an ODE with a solution of the form:

$$ \\psi (x,t_0 + \\Delta t / 2) = \\exp (i \\Delta t \\mathcal{N} /2) \\psi (x,t_0) $$

Now neglecting \\(\\mathcal{N} \\), moving to momentum space \\(\\mathcal{L} \\) is simply multiplication. Hence in the full time interval \\(\\Delta t\\):

$$  \\tilde{ \\psi } (k,t_0 + \\Delta t) = \\exp (i \\Delta t \\mathcal{F} ( \\mathcal{L})) \\tilde{ \\psi }(k,t_0) = \\exp (-i \\Delta t k^2 / 2 m) \\tilde{ \\psi }(k,t_0) $$

For the initial \\( \\tilde{ \\psi }(k,t_0)\\) we use the Fourier transform of the time half step result we found first. Finally we must perform an additional spatial domain time half step to recover the split step approximation to time evolution operator for \\(\\mathcal{L} + \\mathcal{N}\\) by \\( \\Delta t\\).

In full, the process is the following:

$$ \\psi (x,t_0+ \\Delta t) = \\exp (i \\Delta t \\mathcal{N} /2) \\mathcal{F}^{-1}( \\exp (i \\Delta t \\mathcal{F}( \\mathcal{L})) \ \\mathcal{F} (\\exp(i \\Delta t \\mathcal{N} /2) \ \\psi(x,t_0))) $$

We will be using Fast Fourier Transforms (FFTs) from the SciPy library so need to take into consideration the discrete nature of our input.

The basic argument behind this is to match the continuous Fourier transform pair \\( \\psi(x,t) \\leftrightarrow \\tilde{ \\psi} (k,t)\\) to a discrete approximation, \\( \\psi(x_n,t) \\leftrightarrow \\tilde{ \\psi} (k_m,t)\\). Here we use n and m to index x and k:

Starting with a continuous Fourier transform, we can form the discrete approximation:

$$ \\tilde \\psi (k,t) = \\frac{1}{ \\sqrt{2 \\pi}} \\int^\\infty _\\infty \\psi (x,t) e^{-ikx} dx \\to \\tilde \\psi (k_m,t) \\approx \\frac{ \\Delta x}{ \\sqrt{2 \\pi}} \\sum^{N-1}_{n=0} \\psi(x_n,t) e^{-ik_mx_n} $$

Comparing these to the discrete Fourier transform definitions we find the discrete Fourier transform pair:

$$ \\frac{ \\Delta x}{ \\sqrt{2 \\pi}} \\psi(x_n,t) e^{-ik_0x_n} \\leftrightarrow \\tilde \\psi (k_m,t) e^{im \\Delta k x_0} $$

Where \\( \\Delta k = 2 \\pi / (N \\Delta x) \\)

Note that just as we have limited the range of x above, we have here limited the range of k as well. This means that high-frequency components of the signal will be lost in our approximation. The Nyquist sampling theorem tells us that this is an unavoidable consequence of choosing discrete steps in space.

*It should be noted that the wavefunctions reaching the boundaries should be avoided. The spatial domain should be large enough and the inital wavefunction should be of a suitable form e.g. Gaussian wavepackets*

Gaussian wavepackets
^^^^^^^^^^^^^^^^^

A Gaussian wavepackets can be used to investigate the quantum mechanical evolution of particles as they are well localized in both real and momentum space (i.e. they are minimum uncertainty states). However a Gaussian wavepacket will spread over time so a smart choice of a initial wavefunction will help the wavepacket last long enough to observe long times. The width of the probability density increases with time as follows:

For a initial wavepacket of the form:

$$ \\psi (x, t = t_0) \\propto \\exp \\left(- \\frac{(x-x_0)^2}{4 \\sigma_0^2 } + i k_0 x \\right) $$

$$ \\sigma^2 (t) = \\sigma^2_0 + \\frac{1}{\\sigma^2_0} \\frac{t^2}{4 m^2} $$

Hence the width of the packet is minimized at the end of the simulation when

$$ \\frac{d \\sigma}{\\sigma_0} = 0 \\to \\sigma^2_0 = \\frac{t}{2m} $$

Below is an example of setting up such a wavepacket:

.. code-block:: python

 def oneD_gaussian(x,mean,std,k0):
     return np.exp(-((x-mean)**2)/(4*std**2)+ 1j*x*k0)/np.sqrt(4*np.pi*std**2)

 dt = 0.01
 N_t = 2000

 p0 = 2.0
 d = np.sqrt(N_t*dt/2.)

 psi_0 = oneD_gaussian(x,x.max()-10*d,d,-p0)

Non-Linear Schrödinger
^^^^^^^^^^^^^^^^^^^^

The non-linear Schrödinger equation includes a term which depends on the probability density of the wavefunction. This can be included by modifying our \\(\\mathcal{N}\\) operator:

$$ \\mathcal{N} = - V(x) + \\kappa \\left| \\psi (x,t) \\right| ^2 $$

Depending on the sign, this corresponds to a repulsive or attractive contact potential between particles described by the wavefunction.

Argument list
^^^^^^^^^^^^

split_step_schrodinger(psi_0, dx, dt, V, N_t, x_0 = 0., k_0 = None, m = 1.0, non_linear = False)

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
