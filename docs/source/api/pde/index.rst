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

numerov(x,dx,V,E,initial_values,params)
^^^^^^^^^^^^^^^^^^^

   This function performs a Numerov integration at the given energy within the given domain and returns the un-normalised wavefunction evaluated at the given spatial coords

   **Parameters:**

   *x: numpy array*

   An N element numpy array of equally spaced points in space (creating using numpy linspace is advised) at which the wavefunction will be evaluated

   *dx: float*

   Must give the spacing between points in the x array
   
   *V: function*
   
   Pass a function which takes x as an argument and returns the value of potential at that point, V(x)
   
   *E: float*
   
   The energy of the time independent wavefunction. This will give exponential growth if E does not correspond to a bound or free state energy eigenvalue.
   
   *inital_value: list*
   
   A list with 2 elements. The first of these is the boundary condition for the wavefunction at the lowest value spatial coordinate. The second of these is used to initialise the next point in from this, it must be non-zero but apart from this requirement it only effects normalisation (see in detail documentation)
   
   *params: list*
   
   List which can be used within your code to hold various physical parameter. The first element is required to be equal to the particle mass but apart from this size of the list and the other parameters are not called
   
   **Returns:**

   A *numpy array* containing the approximated wavefunction evaluated at x