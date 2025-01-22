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
date: 17 December 2019
bibliography: paper.bib
---

# Summary

Hypersonic flows present unique challenges due to the complex interplay of fluid dynamics, chemical reactions, and optical phenomena. As a signal from a Light Detection and Ranging (LiDAR) system travels through a hypersonic flow field, the beam is affected by the flow, potentially leading to errors in targeting and detection measurements.

`HAOT` is a Hypersonic Aerodynamics Optics Tools Python package developed to calculate different aerodynamic properties in a hypersonic medium. Its source code is available on [GitHub](https://github.com/mliza/HAOT), the documentation is available on [Read the Docs](https://haot.readthedocs.io/en/latest/) and an example on the usage of the package is given on the GitHub repo under the example folder. 

# Statement of Need

Many techniques used to calculate optical properties are scattered across various papers, but there is no centralized repository containing all these calculations. Furthermore, some of these calculations require spectroscopy constants, which are often unclear or inconsistently presented in the literature. This package includes a constants module that provides and documents numerous spectroscopy constants for diatomic molecules. 

This package was used by [@Liza2023]. In this work, they used results from
Computational Fluid Dynamics (CFD) to calculate optical properties of
interest.

# Algorithms
The `HAOT` package, contains five modules:

    - Aerodynamics
    - Optics
    - Quantum Mechanics
    - Constants
    - Conversions

Each module can be imported independently. The [documentation](https://haot.readthedocs.io/en/latest/) explains the functions in each module as well as their usage. Docstrings are include, so the function prototypes and usage can also be accessible in an interactive Python session. Results from these algorithms were compared with the literature, and a unit test was developed, which is located under the unit_test directory. 

This section highlights some of the capabilities of the package but not all of them. For instance, the package can calculate various compressible flow properties, such as isentropic, normal shock, and oblique shock relations. Please refer to the [documentation](https://haot.readthedocs.io/en/latest) for a complete list of available functions.

# Results 

Equation \ref{eq:indexAtmosphere} was introduced by [@Smith1953], and it provides an approximation for the index of refraction as a function of atmospheric altitude.

\begin{equation}\label{eq:indexAtmosphere}
n(h) \approx 1 + \frac{K_1}{T(h)} \left( p(h) + K_2\frac{e(h)}{T(h)} \right) 
\end{equation}

Where: $K_1$ and $K_2$ are constants, $T(h)$ is the temperature as a function of altitude, $p(h)$ is pressure as a function of altitude, and $e(h)$ is the partial pressure of water vapor.

Results for Equation \ref{eq:indexAtmosphere} are provided in Figure \ref{fig:atmosphericIndexOfRefraction}. This approximation is a useful way of analyzing the region in which the index of refraction has the greatest impact on a seeker's performance. As expected, the critical region is between $20~\mathrm{[km]}$ and $30~\mathrm{[km]}$ above sea level, which is a region where a seeker's performance is particularly important.

\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{atmosphericOptics.png}
    \caption{Atmospheric index of refraction for dry air. \label{fig:atmosphericIndexOfRefraction}}
\end{figure}

The Gladstone-Dale constant is an important constant used to calculate the index of refraction (equation \ref{eq:indexGD}). For a dilute gas, the index of refraction can be approximated as:

\begin{equation}\label{eq:indexGD}
n - 1 = \sum\limits_{s = 1}^{N} K_s \rho_s
\end{equation}

Where: $n$ is the index of refraction, $K_s$ is the species' Gladstone-Dale
constants, and $\rho_s$ is the species density. 

Figure \ref{fig:speciesGladstoneDale} shows the Gladstone-Dale constants for a five-species gas. The species density was calculated using the CFD tool, SU2 [@Maier2021], [@Maier2023a]. These densities were then loaded into the `HAOT` tool to calculate Figure \ref{fig:speciesGladstoneDale}.
\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{3C_speciesGladstoneDale.png}
    \caption{Species Gladstone-Dale constants for a five-species gas.\label{fig:speciesGladstoneDale}}
\end{figure}
Another way of calculating the index of refraction is using the polarizability, as shown in equation \ref{eq:indexPolarizability}.
\begin{equation}\label{eq:indexPolarizability}
    n - 1 = \frac{1}{2\epsilon_0}\sum\limits_{s = 1}\alpha_s N_s
\end{equation}
Where: $\epsilon_0$ is the dielectric constant in vacuum, $\alpha_s$ is the species' polarizability constant, and $N_s$ is the partial species mass fraction.

Figure \ref{fig:kerlPolarizability} uses the extrapolation (equation \ref{eq:kerlExtrapolation}) developed by [@Kerl1992], based on the work of [@Hohm1986].

\begin{equation}\label{eq:kerlExtrapolation}
\alpha(\omega, T) = \frac{\alpha(0,0)}{1 - \left(\frac{\omega}{\omega_0}\right)^2} \left( 1 + bT + cT^2\right)
\end{equation}
Where: $\alpha$ is the polarizability as a function of the laser's frequency $\omega$ and temperature $T$. $b$ and $c$ are extrapolation constants, all of which are implemented in the `HAOT` package, and $\omega_0$ is the oscillation frequency of the diatomic molecule.

Figure \ref{fig:kerlPolarizability} shows the change of polarizability as a
function of temperature.
\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{kerlPolarizability_O2_633nm.png}
    \caption{Polarizability of $O2$ using Kerl's method for a laser of
    $633[nm]$.\label{fig:kerlPolarizability}}
\end{figure}

The results presented so far highlight some of the capabilities of the aerodynamics and optics modules. Figure \ref{fig:boltzmannDistribution} shows the results of the Boltzmann Distribution for $N_2$ (in the quantum mechanics module). Some aero-optics calculations [@Buldakov2000], [@Tropina2018] required a Boltzmann Distribution to calculate the polarizability and index of refraction for diatomic molecules.

\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{boltzmannDistribution_N2.png}
    \caption{Boltzmann Distribution for $N2$.\label{fig:boltzmannDistribution}}
\end{figure}

These tools can be applied to the results from a CFD code fairly easy results of calculating the index of refraction, Gladstone Dale Constant and dielectric medium are below. A low speed compressible Large Eddys Simulation was used with OpenFOAM to perform results below, the goal of this is not to go into very especific details of the CFD ssetup but to showcase the ablity of the HAOT tool. 

# References
