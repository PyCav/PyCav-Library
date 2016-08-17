Numerov Method
====================================

In units where \\(\\hbar = 1 \\), the 1D time-independent Schrödinger equation can be expressed in the form:

$$ \\frac{d^2 \\psi}{dx^2} = -2m\\left(E - V(x) \\right) \\psi = -g(x) \\psi $$


This differential equation can be solved numerically via Numerov's method (see_ pages 10 - 11). For a 1D spatial grid, the wavefunction at the \\(n+1\\)th point along the x-axis can be approximated by:

.. _see: http://www.fisica.uniud.it/~giannozz/Corsi/MQ/LectureNotes/mq-cap1.pdf

$$ \\psi_{n+1} = \\frac{(12-10f_n) \\ \\psi_n-f_{n-1}\\psi_{n-1}}{f_{n+1}} $$

where:
$$ f_n \\equiv \\left( 1 + \\frac{\\delta x^2}{12}g_n \\right), \ \ \ \ \ \ \ g_n = 2m(E-V(x_n)) $$

Hence to start Numerov's method we require \\(\\psi_0\\) and \\(\\psi_1\\), in other words \\(\\psi(x = x_{min})\\) and \\(\\psi(x = x_{min} + \\delta x)\\), where \\(\\delta x\\) is the step size.

When investigating bound states we require \\(\\psi( x \\to \\pm \\infty) = 0\\). However, we cannot consider an infinite domain. Instead we must choose a large enough domain that setting \\(\\psi( x = x_{min}) = 0\\) is a good approximation (and similarly for \\(x_{max}\\)).

With Numerov's Method in place, the shooting method can be used to find the energy eigenstates. It goes as follows:

1. Setting \\(\\psi_0 = 0\\) approximately satisfies the boundary condition that the wavefunctions must vanish at the boundary.
2. Since the Schrödinger Equation is linear and homogeneous we are free to set \\(\\psi_1\\) to any non-zero constant as multiplying by a constant does not affect the solution. In this case we set \\(\\psi_1 = \\delta x\\).
3. Using the Numerov algorithm, \\(\\psi(x)\\) can be found. Exponential growth near \\(x_{max}\\) is observed if the input energy is not near a energy eigenvalue

Argument list
^^^^^^^^

numerov(x,dx,V,E,initial_values,params)

   This function performs a Numerov integration at the given energy within the given domain and returns the un-normalised wavefunction evaluates it over the whole domain.

   **Parameters:**

   *x: numpy array*

   An N element numpy array of equally spaced points (creating using numpy linspace is advised) at which the wavefunction will be evaluated

   *dx: float*

   The spacing between points in the x array
   
   *V: function*
   
   Function which takes x as an argument and returns the value of potential at that point, V(x)
   
   *E: float*
   
   The energy of the time independent wavefunction. This will give exponential growth if E does not correspond to a bound or free state energy eigenvalue.
   
   *inital_value: list*
   
   A list with 2 elements. The first of these is the boundary condition for the wavefunction at the lowest value spatial coordinate. The second is used to initialise the next point in from this, which must be non-zero. Apart from this requirement, the second element only effects normalisation.
   
   *params: list*
   
   List which can be used within your code to hold various physical parameters. The first element must be equal to the particle mass, but apart from this the size of the list and the other parameters are not called
   
   **Returns:**

   A *numpy array* containing the approximated wavefunction evaluated at x

