Spring
======

The Spring class represents a spring. It is not tied to any visualization method by design. It connects two Particles together and applies a force F = kx to them, as per Hooke's Law. It also shares many visualization properties with the Particle class.

Functions
---------

__init__(particle_1, particle_2, k, l0=None, radius=0.5, color=None, alpha=1.)
^^^^^^^^^^^^^^^^^^
	Initialises the Spring object by supplying the 2 particles it connects and the 	value of the stiffness, k.
	
	**Parameters:**
	
	*particle_1: Particle*
	
	Particle on one end of spring
	
	*particle_2: Particle*
	
	Particle on other end of spring
	
	*k: float*
	
	The stiffness of the spring (F = kx)

	*l0: float*

	Original length of the spring
	
	*radius: float*
	
	Radius of spring.
	
	*color: array*
	
	Color of particle, given in form [R G B]
	
	*alpha: float*
	
	Alpha of particle, 1 is completely opaque, 0 is completely transparent, used in 	visualisation

force_on(particle, if_at=np.array([None]))
^^^^^^^^^^^^^^^^^^^
	Given an arbitary particle, gives the force on that particle. No force if the 	spring isn't connected to that particle.
	
	**Parameters:**
	
	*particle: Particle*
	
	Particle which feels the force
	
	*if_at: NumPy Array*
	
	If this parameter is used, this gives the force felt if the particle were at this position.

	**Returns:**

	A *numpy array* with 3 elements giving the vector force on the particle from the spring.


Properties
----------

particle_1
^^^^^
  *Particle*

  First particle that the spring is attached to.

particle_2
^^^^^
  *Particle*

  Second particle that the spring is attached to.

k
^^^^
	*float*

	The stiffness of the spring (F = kx)

l0
^^^^^
	*float*

	The original length of the spring.

color
^^^^^
  *array*

  Gives the color of the spring as an array, given in the form [R G B]. Each element of the array should be between 0 and 1.

radius
^^^^^^
  *float*

  Gives the radius of the spring.

alpha
^^^^^
  *float*

  Float between 0 and 1, giving the opacity of the spring.

pos
^^^
  *numpy array, read only*

  3 element array giving position of one end of the spring in 3D space.

axis
^
  *numpy array, read only*

  3 element array the axis, i.e. the vector showing the orientation and length of the spring.