PointerArrow
============

PointerArrow is a class that represents an arrow/pointer. It is not tied to any visualization method by design.

Functions
--------

__init__(pos, axis, shaftwidth=1, color=None, alpha=1)
^^^^^^^^^^^^^^^^^^^^^

	Initialises the PointerArrow object

	**Parameters:**

	*pos: numpy array*

  Location of tail end of PointerArrow

	*axis: numpy array*

	A 3 dimensional vector giving the length and orientation of the pointer

	*shaftwidth: float*

  The width of the pointer's shaft

	*color: array*

	Color of PointerArrow, given in form [R G B], default 1, 1, 1

	*alpha: float*

  Alpha of PointerArrow, 1 is completely opaque, 0 is completely transparent, used in visualisation

Properties
----------

pos
^^^
  *numpy array*

  3 element array giving position of the centre of the PointerArrow in 3D space.


color
^^^^^
  *array*

  Gives the color of the PointerArrow as an array, given in the form [R G B]. Each element of the array should be between 0 and 1.

shaftwidth
^^^^^^
  *float*

  Gives the shaft width of the PointerArrow.

alpha
^^^^^
  *float*

  Float between 0 and 1, giving the opacity of the PointerArrow.