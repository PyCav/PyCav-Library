Time-Dependent Schrödinger equation via the Split-Step Fourier method:
==============================

Writing the Schrodinger equation in the form:
$$ \\frac{ \\partial \\psi}{\\partial t} = i \\mathcal{L} \\psi+i \\mathcal{N} \\psi $$

Where for the TDSE:
$$ \\mathcal{L} = -\\frac{1}{2} \\frac{\\partial^2}{\\partial x^2},  \\mathcal{N} = -V(x) $$

By first neglecting \\(\\mathcal{L} \\) in time interval \\([t_0,t_0+ \\Delta t/2]\\) we are left with an ODE with a solution of the form:

$$ \\psi (x,t_0 + \\Delta t) = \\exp (i \\Delta t \\mathcal{N} /2) \\psi (x,t_0) $$

Now neglecting \\(\\mathcal{N} \\), moving to momentum space \\(\\mathcal{L} \\) is simply multiplication. Hence in the full time interval $\Delta t$:

$$  \\tilde{ \\psi } (k,t_0 + \\Delta t) = \\exp (i \\Delta t \\mathcal{F} ( \\mathcal{L})) \\tilde{ \\psi }(k,t_0) = \\exp (-i \\Delta t k^2) \\tilde{ \\psi }(k,t_0) $$

For the initial \\( \\tilde{ \\psi }(k,t_0)\\) we use the Fourier transform of the time half step result we found first. Finally we must perform an additional spatial domain time half step to recover the split step approximation to time evolution by \\( \\Delta t\\).

In full, the process is the following:

$$ \\psi (x,t_0+ \\Delta t) = \\exp (i \\Delta t \\mathcal{N} /2) \\mathcal{F}^{-1}( \\exp (i \\Delta t \\mathcal{F}( \\mathcal{L})) \ \\mathcal{F} (\\exp(i \\Delta t \\mathcal{N} /2) \ \\psi(x,t_0))) $$

We will be using Fast Fourier Transforms (FFTs) from the SciPy library so need to take into consideration the discrete nature of our input.

The basic argument behind this is to match the continuous Fourier transform pair \\( \\psi(x,t) \\leftrightarrow \\tilde{ \\psi} (k,t)$ to a discrete approximation, \\( \\psi(x_n,t) \\leftrightarrow \\tilde{ \\psi} (k_m,t)\\).

$$ \\tilde \\psi (k,t) = \\frac{1}{ \\sqrt{2 \\pi}} \\int^\\infty _\\infty \\psi (x,t) e^{-ikx} dx \\to \\tilde \\psi (k_m,t) \\approx \\frac{ \\Delta x}{ \\sqrt{2 \\pi}} \\sum^{N-1}_{n=0} \\psi(x_n,t) e^{-ik_mx_n} $$

Comparing these to the discrete Fourier transform definitions we find the discrete Fourier transform pair:

$$ \\frac{ \\Delta x}{ \\sqrt{2 \\pi}} \\psi(x_n,t) e^{-ik_0x_n} \\leftrightarrow \\tilde \\psi (k_m,t) e^{im \\Delta k x_0} $$

Where \\( \\Delta k = 2 \\pi / (N \\Delta x) \\)