Installing PyCav
================

PyCav requires Python 2.7 or ≥ 3.5.

There are a few ways to install PyCav.
* The simplest and recommended way is to run the following command in Command Prompt (Windows) or Terminal(Mac/Linux), with pip installed:

.. code-block:: bash

    $ pip install pycav

  This way, PyCav and all of its dependencies will be installed.
  If you want to visualize the simulations, vpython also needs to be installed with

.. code-block:: bash

    $ pip install vpython

* If you don't want to use the previous command, you can get the latest release as a tarball from `from PyPI <https://pypi.python.org/>`_. With this method, however, you will not get the dependencies automatically installed, so make sure that you have:
  * numpy
  * matplotlib
  * scipy
  * vpython (If you want to visualize simulations)
  installed when you want to use PyCav.