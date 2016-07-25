Perturbation Analysis
====================================

Using the following notation for our perturbation analysis:

Total Hamiltonian:
$$ \\hat{H} = \\hat{H}^{(0)}+\\hat{H}'$$

First order energy shift:
$$ \\Delta E_n^{(1)} = \\langle n^{(0)} | \\hat{H}' | n^{(0)} \\rangle$$

First order perturbed wavefunctions:
$$ | n^{(1)} \\rangle \\approx | n^{(0)} \\rangle +  \\sum_{m \\neq n} | m^{(0)} \\rangle \\frac{\\langle m^{(0)} | \\hat{H}' | n^{(0)} \\rangle}{E_n^{(0)} - E_m^{(0)}} $$

Evaluating inner products are done using an integration over space i.e.
$$ \\langle n^{(0)} | \\hat{H}' | n^{(0)} \\rangle = \\int_{x_{min}}^{x_{max}} \\psi_{n}^{(0)}(x) \\hat{H}'(x) \\psi_{n}^{(0)}(x) dx $$

These are calculated using SciPy's quad integration function.

For first_order_wf, \\(I_{mn} = \\langle m^{(0)} | \\hat{H}' | n^{(0)} \\rangle / (E_n^{(0)}-E_m^{(0)})\\) is calculated for m values around n until \\(I_{mn} < \\epsilon\\), where \\(\\epsilon\\) is the given tolerance. Two seperate iterations are run for even and odd values of m around n as the form of the perturbation may cause inner products to vanish for certain configurations of even/odd wavefunctions. 

If return_list is set to True, the list of m values used in the sum is returned along with the corresponding \\(I_{mn}\\) values.

The perturbed wavefunction is calculated within the function using the following sum:
$$ | n^{(1)} \\rangle \\approx | n^{(0)} \\rangle + \\sum_{k \\in k_{list}} I_{kn} \ | k^{(0)} \\rangle $$
This is returned as a function of position i.e.

.. code-block:: python
  
   x = np.linspace(-10.,10.,100)  
   perturb_wf = first_order_wf(n,H,unperturb_wf,unperturb_erg,params)
   perturb_wf_x = perturb_wf(x)

Here perturb_wf_x is a numpy array containing the perturbed wavefunction evaluated at all the points in x
   

first_order_energy_sft(n,H,unperturb_wf,params,limits = [-np.inf,np.inf])
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   
   Works out the 1st order energy shift from time independent perturbation theory for a given unperturbed system and the applied perturbation

   **Parameters:**

   *n: int or list*

   Principal quantum number (or list of these) which labels the unperturbed wavefunctions
   
   *H: function*

   The applied perturbation as a function of position 
   
   *unperturb_wf: function*

   A function which is passed params and n and returns a function of position e.g.
   
   .. code-block:: python
   
    def unperturb_wf(params,n):
      a = params[0]
      m = n+1
      def psi_n(x):
        return np.sqrt(2./a)*np.sin((m*np.pi/a)*(x+a/2.))
      return psi_n
      
   Here psi_n(x) is the returned function of position
  
   *params: list*
   
   List which can be used within your code to hold various physical parameter used by unperturb_wf and other functions (see later)
  
   *limits: list*
  
   List containing the integration limits for the inner product of the wavefunctions and the perturbation. Default is the whole space. For hard wall potentials adjust these limits to respect the boundaries

   **Returns:**

   A list or float depending on the input argument n containing the first energy shifts to the unperturbed wavefunctions for the given quantum numbers in n

first_order_wf(n,H,unperturb_wf,unperturb_erg,params,tolerance = 0.01, limits = [-np.inf,np.inf], return_list = False)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   Calculates the 1st order perturbed wavefunction for a given unperturbed system and the applied perturbation. The system is defined by its known unperturbed wavefunctions and energies. This function takes similar parameters as first_order_energy_sft (see above) so only new parameters will be defined.

   **Unique Parameters:**

   *unperturb_erg: function*

   A function which is passed params and n and returns the energy of the n-th unperturbed eigenstate e.g. for a harmonic oscillator

   .. code-block:: python
  
    def unperturb_erg(params,n):
      return (n+0.5)*params[1]

   where params[1] contains the angular frequency (for hbar = 1)

   *tolerance: float*

   The value below which terms in the 1st order wavefunction sum are ignored

   *return_list: float*

   Set to True if you require the perturbation sum prefactors and values of the principal quantum numbers of the unperturbed wavefunctions

   **Returns:**

   A function of position which corresponds to the 1st order perturbed wavefunction and if return_list = True, copies of the principal quantum number lists and the sum prefactors list which were used to calculate the resultant perturbed wavefunction