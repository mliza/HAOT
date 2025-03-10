OpenFOAM
========
We will not delve into the details of setting up a CFD simulation; instead, the following example will use results generated with OpenFOAM. This example is using turbulent flow in equilibrium.

For this example, an additional module, ``pyvista``, is required to parse the CFD results data. Please follow the instructions in the next section to install ``pyvista``. If you already have ``pyvista`` installed, you can skip the next code lines.

.. code:: console

    $ pip install pyvista


Next, lets load the `*.foam` file. For additional details, please refer to the `OpenFOAM documentation <https://www.openfoam.com>`_ and the `PyVista documentation <https://pyvista.org/>`_. The following lines of code will load the foam data using ``pyvista``. 

.. code:: python

    import pyvista as pv

    # Loading results data
    reader = pv.POpenFOAMReader('results.foam')
    mesh = reader.read()
    internal_mesh = mesh['internalMesh']

Lets now use the ``HAOT`` package to process the results

.. code:: python

    import haot

    # Compute index of refraction as a function of temperature and density
    index_of_refraction = haot.index_of_refraction_density_temperature(
                                            internal_mesh['T'],
                                            internal_mesh['rho'],
                                            'Air', 633)

    # Compute Kerl polarizability for air as a function of temperature
    kerl_polarizability = haot.kerl_polarizability_temperature(
                                            internal_mesh['T'],
                                            'Air', 633)

    # Compute permittivity of the medium using dilute index of refraction
    permittivity_dilute = haot.permittivity_material(index_of_refraction['dilute'])

    # Compute electric susceptibility using dilute index of refraction
    susceptibility_dilute = haot.electric_susceptibility(index_of_refraction['dilute'])

.. toctree::
   :maxdepth: 1
   :caption: Contour Plots:

   openFoamContour
