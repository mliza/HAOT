OpenFOAM
========
We will not delve into the details of setting up a CFD simulation; instead, the following example will use results generated with OpenFOAM.


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


Now, let's plot the dilute index of refraction using ``pyvista``.

.. code:: python

    # Add the index of refraction to the mesh
    internal_mesh.cell_data['n'] = index_of_refraction['dilute']

    plotter = pv.Plotter(window_size=[1800, 900])
    plotter.view_xy()
    plotter.add_mesh(internal_mesh, scalars='n', cmap='turbo',
                     reset_camera='True', show_scalar_bar=False)
    plotter.set_background('white')
    plotter.camera.zoom(2.0)

    plotter.add_scalar_bar(
        title='Dilute Index of refraction',
        title_font_size=22,
        label_font_size=18,
        bold=True,
        position_x=0.02,
        position_y=0.6,
        width=0.3,
        n_labels=8,
        height=0.1,
        vertical=False,
        fmt=""
    )

    plotter.show()

.. image:: images/index_of_refraction_dilute.png

Now, let's plot the Kerl polarizability using ``pyvista``.

.. code:: python

    # Add polarizability to the mesh
    internal_mesh.cell_data['pol'] = kerl_polarizability

    plotter = pv.Plotter(window_size=[1800, 900])
    plotter.view_xy()
    plotter.add_mesh(internal_mesh, scalars='pol', cmap='turbo',
                     reset_camera='True', show_scalar_bar=False)
    plotter.set_background('white')
    plotter.camera.zoom(2.0)

    plotter.add_scalar_bar(
        title='Polarizability',
        title_font_size=22,
        label_font_size=18,
        bold=True,
        position_x=0.02,
        position_y=0.6,
        width=0.3,
        n_labels=8,
        height=0.1,
        vertical=False,
        fmt=""
    )

    plotter.show()

.. image:: images/polarizability_kerl.png


Now, let's plot the permittivity of the medium using ``pyvista``.

.. code:: python

    # Add Permittivity constant to the mesh
    internal_mesh.cell_data['permittivity_dilute'] = permittivity_dilute

    plotter = pv.Plotter(window_size=[1800, 900])
    plotter.view_xy()
    plotter.add_mesh(internal_mesh, scalars='permittivity_dilute', cmap='turbo',
                     reset_camera='True', show_scalar_bar=False)
    plotter.set_background('white')
    plotter.camera.zoom(2.0)

    plotter.add_scalar_bar(
        title='Permittivity',
        title_font_size=22,
        label_font_size=18,
        bold=True,
        position_x=0.02,
        position_y=0.6,
        width=0.3,
        n_labels=8,
        height=0.1,
        vertical=False,
        fmt=""
    )

    plotter.show()

.. image:: images/permittivity.png

Now, let's plot the electric susceptibility using ``pyvista``.

.. code:: python

    # Add Electric Susceptibility constant to the mesh
    internal_mesh.cell_data['susceptibility_dilute'] = susceptibility_dilute

    plotter = pv.Plotter(window_size=[1800, 900])
    plotter.view_xy()
    plotter.add_mesh(internal_mesh, scalars='susceptibility_dilute', cmap='turbo',
                     reset_camera='True', show_scalar_bar=False)
    plotter.set_background('white')
    plotter.camera.zoom(2.0)

    plotter.add_scalar_bar(
        title='Susceptibility',
        title_font_size=22,
        label_font_size=18,
        bold=True,
        position_x=0.02,
        position_y=0.6,
        width=0.3,
        n_labels=8,
        height=0.1,
        vertical=False,
        fmt=""
    )

    plotter.show()

.. image:: images/susceptibility.png
