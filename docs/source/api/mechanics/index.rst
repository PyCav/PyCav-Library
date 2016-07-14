Mechanics module API Reference
====================================

Introduction to the mechanics module
--------------------------

The mechanics module provides the infrastructure to create simple systems consisting of springs and particles, and handles collisions between them as well as the visualization. A quick example of their use:

.. code-block:: python
	
	from pycav.mechanics import *
	from vpython import *

	container = Container(1)
	system = System(collides = True, interacts = False, visualize = True, container = 	container)
	avg_speed = 1
	system.create_particles_in_container(number = 100, speed = avg_speed, radius = 0.03)
	system.run_for(10)

This creates a system consisting of a 100 particles which can collide with each other within a cubic container of unit size, which is then run for 10 seconds in the system's time.

Classes
-------

.. toctree::
   :maxdepth: 1

   particle
   spring
   container
   pointerarrow
   system

Functions
---------

no_force(pos, time)
^^^^^^^^^^^^^^^^^^^

   This is used as the default force applied to Particles. Returns no force, as name suggests.

   **Parameters:**

   *pos: numpy array*

   A 3 element numpy array giving the position of the object feeling the force

   *time: float*

   Time at which this force is felt
   
   **Returns:**

   A *numpy array* 0., 0., 0.


element_mult(vec_1, vec_2)
^^^^^^^^^^^^^^^^^^^^^^^^^^
   
   Performs element wise multiplication (computes the Hadamard product) between two 3 dimensional vectors.

   **Parameters:**

   *vec_1: numpy array*

   A 3 element numpy array which represents the first vector for which elements are multiplied

   *vec_2: numpy array*

   A 3 element numpy array which represents the second vector for which elements are multiplied

   **Returns:**

   The Hadamard product of the two vectors given as a 3-element numpy array


vector_from(arr=np.array([0, 0, 0]))
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   
   Creates a vpython vector from a numpy array.

   **Parameters:**

   *arr: numpy array*

   3 element numpy array which will be converted to into a vpython vector

   **Returns:**

   The numpy array converted into the equivalent vpython array

normalized(arr)
^^^^^^^^^^^^^^^

   Creates a normalised 1-D array for the 1-D array given.

   **Parameters:**

   *arr: numpy array*

   Array to be normalised

   **Returns:**

   A 1-D numpy array that is of unit length