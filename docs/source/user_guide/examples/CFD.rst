CFD
===

Lets do a more complicated example. Lets use the results from a Computational
Fluid Dynamics Simulation to calculate some optical properties. We will not go
into the details on setting up a CFD simulation, but the following result will
be using results from an OpenFOAM. 

For this example an additional module to parse the CFD results data  will have to be installed ``pyvista`` please follow the folowing lines. Skip next code lines if you already have ``pyvista``.

.. code:: console

    $ pip install pyvista


Lets now load the `*.foam` file. Please refer to hyperlink:`OpenFOAM documentation <https://www.openfoam.com>`_.

.. code:: python

