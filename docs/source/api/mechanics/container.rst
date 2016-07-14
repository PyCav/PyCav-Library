Container
=========

Container is a class that represents a cubic Container that particles can be inside. It is not tied to any visualization method by design. 

Functions
---------

__init__(dimension, pos=None, color=None, alpha=0.3)
^^^^^^^^^^^^^^^^^^^^^^^^
	
	Initialises the Container object

	**Parameters:**

	*dimension: float*

  The dimension of the Container, i.e. the length, width, height of the Container

	*pos: numpy array*

	Position of centre of cube, default 0, 0, 0

	*color: array*

	Color of particle, given in form [R G B], default 1,1,1

	*alpha: float*

	Alpha of particle, 1 is completely opaque, 0 is completely transparent, used in visualisation

contains(particle)
^^^^^^^^^^^^^^
	
	Function which checks if a given particle is entirely inside the Container

	**Parameters:**

	*particle: Particle*

	Particle which is being checked to see if inside Container or not

	**Returns:**

	True if the Container contains the particle, and if not, returns the index of the axis along which the particle is outside the Container

Properties
----------

pos
^^^
  *numpy array*

  3 element array giving position of tail end of Container in 3D space.

dimension
^^^^^^^^^
	*float*

	The dimension of the Container, i.e. the length, width, height of the Container

axis
^^^^^
  *numpy array*

  3 element array the axis, i.e. the vector showing the orientation and length of the Container.

color
^^^^^
  *array*

  Gives the color of the Container as an array, given in the form [R G B]. Each element of the array should be between 0 and 1.

alpha
^^^^^
  *float*

  Float between 0 and 1, giving the opacity of the Container.

surface_area
^^^^^^^^^^^
	*float*

	Float giving the surface area of the Container.