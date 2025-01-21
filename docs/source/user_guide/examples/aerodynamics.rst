Aerodynamics
============
Let pick a problem from the [Anderson]_ textbook.

Consider a normal shock wave in air. The upstream conditions are given by
:math:`M_1 = 3, P_1 = 101 [kPa]` and :math:`\rho_1=1.23 [kg/m3]`. Calculate the
downstream values :math:`M_2, P_2, rho_2`.

Open your favorite python instance, for this case I will open IPython and in
IPython I will add the following lines.

.. code:: python

    import haot

    downstream = haot.normal_shock_relations(3.0)

    dowstream # Print dictionary values

    >> {'mach_2': 0.4751909633114914,
        'pressure_r': 10.333333333333334,
        'temperature_r': 2.6790123456790123,
        'density_r': 3.857142857142857,
        'pressure_s': 0.32834388819073684}

The output is a dictionary, but the variable names may not be very intuitive. To gather more information, you can use Python's help() function or consult the documentation in the modules section under the :ref:`Aerodynamics module<Module aerodynamics target>`. However, for this example, we will rely solely on standard Python techniques.

.. code:: python

    help(haot.normal_shock_relations)

    >> Calculates normal shock relations

        Parameters:
            mach_1: pre-shock mach number
            adiabatic_indx: adiabatic index, 1.4 (default)

        Returns:
            dict: A dictionary containing:
                - mach_2: post-shock mach number
                - pressure_r: pressure ratio (post-shock / pre-shock)
                - temperature_r: temperature ratio (post-shock / pre-shock)
                - density_r: density ratio (post-shock / pre-shock)
                - pressure_s: stagnation pressure ratio (post-shock / pre-shock)

        Reference:
            Normal Shock Wave - NASA (https://www.grc.nasa.gov/www/k-12/airplane/normal.html)


Now we can calculate the downstream values by multiplying the upstream values
to the results provided by the module

.. code:: python

    print(downstream['pressure_r'] * 101)

    >> 1043.6666666666667

    print(downstream['density_r'] * 1.23)

    >> 4.744285714285715

    print(downstream['mach_2'])

    >> 0.4751909633114914 


.. [Anderson] Modern Compressible Flow: With Historical Perspective (4th), Aderson J., ISBN: 978-1260471441
