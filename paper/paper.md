---
title: 'HAOT: A Python package for hypersonic aero-optics analysis'
tags:
  - Python
  - Hypersonics 
  - Aerodynamics
  - Optics
authors:
  - name: Martin E. Liza
    orcid: 0009-0000-2231-667X
    affiliation: 1
affiliations:
  - index: 1
    name: The University of Arizona
date: 13 December 2019
bibliography: paper.bib
---

# Summary

Hypersonic flows present a unique challenges due to the complex interplay of fluid dynamics, chemical reactions, and optical phenomena. As a signal from a Light Detection and Ranging (LiDAR) travels through a hypersonic flow field, the beam would be affected by the flow, this can lead to errors on targeting and detection measurements.

`HAOT` is a Hypersonic Aerodynamics Optics Tools Python package developed to calculate different aerodynamic properties in a hypersonic medium. Its source code is available on [GitHub](https://github.com/mliza/HAOT), the documentation is available on [Read the Docs](https://haot.readthedocs.io/en/latest/) and an example on the usage of the package is given on the GitHub repo under the example folder. 

# Statement of Need

Many techniques used to calculate optical properties are scattered across various papers, but there is no centralized repository containing all these calculations. Furthermore, some of these calculations require spectroscopy constants, which are often unclear or inconsistently presented in the literature. This package includes a constants module that provides and documents numerous spectroscopy constants for diatomic molecules. 

This package has been used by [@Liza2023]. In this work he used results from
Computational Fluid Dynamics (CFD) to calculate index optic properties of
interest.

# Algorithms
The `HAOT` package, contains five modules:

    - Aerodynamics
    - Optics
    - Quantum Mechanics
    - Constants
    - Conversions

Each module can be imported independently. The [documentation](https://haot.readthedocs.io/en/latest/) explains he functions in each module as well as their usage. Docstrings were used, so the function prototypes and usage are also available in an interactive Python session. Results from these algorithms were compared with literature and a unit test was developed and it is located under test. 

This section provide some of the capabilities of the packages but not all of them. For instance the package can calculate some compressible flow properties such as isentropic, normal shock, and oblique shock relations. Please refer to the [documentation](https://haot.readthedocs.io/en/latest/) to see the complete list of available functions.

# Results 

The equation below was introduced by [@Smith1953], and it is crude approximation
for the change in the index of refraction as a function of altitude.

$$ n(h) \approx 1 + \frac{K_1}{T(h)} \left( p(h) + K_2\frac{e(h)}{T(h)} \right) \label{eq:atmosphericIndex} $$

Where: $K_1$ and $K_2$ are constants, $T$ is the temperature as a function of altitude, $p$ is pressure as a function of altitude, and $e(h)$ is the partial pressure of water vapor.

Results for this equation are provided in the figure below. As expected the
index of refraction and the atmospheric density have a close relationship.

\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{atmosphericOptics.png}
    \caption{Atmospheric index of refraction for dry air. \label{fig:atmosphericIndexOfRefraction}}
\end{figure}

The Gladstone-Dale constant is an important constant used to calculate the
index of refraction, for a dilute gas, the index of refraction can be
approximated as:

$$ n -1 = \frac{1}{2\epsilon_0}\sum\limits_{i=1}^N \alpha_i N_i $$

Where: $n$ is the index of refraction, $\epsilon_0$ is the dielectric constant
in vacuum, $\alpha_i$ is the species polarizability, and $N_i$ is the species
partial fraction.

Figure \ref{fig:speciesGladstoneDale}
Results for this equation are provided in the figure below. This particular
results required the use of a Computational Fluid Dynamics (CFD) tool, SU2 [@Maier2021], [@Maier2023a], to calculate the fluid properties used by the `HAOT` tool. 
\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{3C_speciesGladstoneDale.png}
    \caption{Species Gladstone-Dale constants for a five-species gas.\label{fig:speciesGladstoneDale}}
\end{figure}


\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{kerlPolarizability_O2_633nm.png}
    \caption{Polarizability of $O2$ using Kerl's method
    [@Kerl1992].\label{fig:kerlPolarizability}}
\end{figure}


\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{boltzmannDistribution_N2.png}
    \caption{Boltzmann Distribution for $N2$.\label{fig:boltzmannDistribution}}
\end{figure}
    




# Acknowledgements
The author gratefully thank Kyle Hanquist, who supported with the tool
verification.

# References
