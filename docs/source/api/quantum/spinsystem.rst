SpinSystem
==========

The SpinSystem class represents a system of spins interacting with one another. 
Interactions with an external inhomogeneous magnetic field can be included as 
an optional argument (note that the gyromagnetic ratios are unity).

Functions
---------
__init__(spins, couplings, B_field = None)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	Initialises the spin system object.

    **Parameters**

    *spins: list or numpy array*

    Spin quantum numbers of the particles (integer or half-integer)
    e.g. [0.5, 1.0]

    *couplings: list or numpy array*

    NxN list for a system of N spins. The element J[i,j] gives the
    coupling strength between spins i and j. Positive coupling favours
    spin alignment.
    e.g. [[0,1],[1,0]]

    *B_field: list, optional*

    List of 3-dimensional vectors, [B_x, B_y, B_z], representing the 
    magnetic flux density at each spin.
    e.g. [[3.0, 0, 0], [1.5, 0, 0]]

create_Hamiltonian()
^^^^^^^^^^^^^^^^^^^^
	Create a Hamiltonian with spin-pairing given by the coupling array. 
	A linear coupling between spins and the magnetic field is included if 
	a B_field was used to initialise the SpinSystem.
    The Hamiltonian acts on Kronecker products of the spins' states.

    This function is called on initialisation and gives the instance of
    SpinSystem the attribute *H*, the Hamiltonian, as a numpy array.

get_energies()
^^^^^^^^^^^^^^
	Calculates the energy levels and eigenstates of the Hamiltonian using
	numpy.linalg.eigh. They are then stored as attributes *energies* and 
	*states* in order of increasing energy, both are numpy arrays.

count_multiplicities(tolerance = 0.0001)
^^^^^^^^^^^^^^^^^^^^^^
	Counts the multiplicities of the energy eigenvalues, generating a list
	of multiplicities and list of non-repeated energies. These lists are 
	subsequently stored as attributes *multiplicities* and *reduced_energies*.

	**Parameters:**

	*tolerance: float, optional*

	Energy levels within *tolerance* of one another are considered degenerate.

Attributes
----------
    *The following attributes are created on initialisation*

    *N: int*

    Number of spins in the system.


    *J: numpy array*

    NxN array of coupling strengths. J[i,j] is the coupling strength
    between spins i and j.

    *Sx: list of numpy arrays*

    A list of x-angular momentum operators corresponding to the
    Hilbert spaces of each spin.

    *Sy and Sz: list of numpy arrays*

    See *Sx*.

    *I: list of numpy arrays*

    A list of the identity operators corresponding to the Hilbert spaces
    of each spin.

    *B: numpy array*

    The magnetic flux densities at each spin in the system. B[i,j] is the
    jth component of the magnetic field at the ith spin.

    *H: numpy array*

    The Hamiltonian of the system.

    *The remaining attributes are not created on initialisation*

    *energies: numpy array*

    Energy eigenvalues of the system, repeated according to degeneracy and
    stored in order of increasing energy.
    *Created by get_energies() and count_multiplicities()*

    *states: numpy array*

    Energy eigenstates, stored in order of increasing energy.
    *Created by get_energies() and count_multiplicities()*

    *multiplicities: list*

    List of the multiplicities of the energy eigenstates in order of
    increasing energy.

    *reduced_eneries: list*

    List of the non-repeated energy levels. e.g. an 8-fold degenerate state
    will feature only once in this list.










  