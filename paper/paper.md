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
date: 30 January 2025
bibliography: paper.bib
---

# Summary

Hypersonic flows present unique challenges due to the complex interplay of fluid dynamics, chemical reactions, and optical phenomena. As a signal from a Light Detection and Ranging (LiDAR) system travels through a hypersonic flow field, the beam is affected by the flow, potentially leading to errors in targeting and detection measurements.

`HAOT` is a Hypersonic Aerodynamics Optics Tools Python package developed to calculate different aerodynamic properties in a hypersonic medium. Its source code is available on [GitHub](https://github.com/mliza/HAOT), the documentation is available on [Read the Docs](https://haot.readthedocs.io/en/latest/) and examples on the usage of the package are given on the documentation site under usage. These examples cover some easy examples such as the basic functionality of the tool and also cover examples using results from Computational Fluid Dynamics (CFD) simulations. 

# Statement of Need

Many techniques for calculating optical properties are dispersed across various research papers, with no centralized repository consolidating these methods. Additionally, some of these calculations depend on spectroscopy constants, which are often inconsistently presented or unclear in the literature. To address this issue, this package includes a dedicated constants module that provides and documents numerous spectroscopy constants for diatomic molecules, ensuring consistency and accessibility.

`HAOT` was developed to facilitate the computation of aerodynamic quantities of interest from CFD results in a streamlined manner. The package has been employed in research, including the work by [@Liza2023], where it was used to investigate nonequilibrium effects on aero-optics in hypersonic flows. By integrating relevant optical and aerodynamic calculations, `HAOT` severs as a valuable tool for researchers serves as a valuable tool for researchers studying the complex interactions between fluid dynamics and optical wave propagation in extreme flow conditions. 

# Algorithms
The `HAOT` package, contains five modules:

- Aerodynamics
- Optics
- Quantum Mechanics
- Constants
- Conversions


Each module can be imported independently, allowing for flexible usage
depending on the specific computational needs. The [documentation](https://haot.readthedocs.io/en/latest/) provides detailed explanations of the functions within each module, along with their proper usage. Additionally, docstrings are included, making function prototypes and usage accessible within an interactive Python session through built-in help functions.

To ensure accuracy and reliability, results obtained from these algorithms have been validated against published literature. A comprehensive suite of unit tests has been developed to detect potential issues and maintain code integrity.

# Results 

This section highlights some of the package’s capabilities but does not encompass all of them. For example `HAOT` can compute various compressible flow properties, including including isentropic flow, normal shock relations, and oblique shock relations. For a complete list of available functions and their descriptions, please refer to the [documentation](https://haot.readthedocs.io/en/latest).

Equation \ref{eq:indexAtmosphere}, introduced by [@Smith1953], provides an approximation for the index of refraction as a function of atmospheric altitude.

\begin{equation}\label{eq:indexAtmosphere}
n(h) \approx 1 + \frac{K_1}{T(h)} \left( p(h) + K_2\frac{e(h)}{T(h)} \right) 
\end{equation}

Where: $K_1$ and $K_2$ are empirically determined constants, $T(h)$ is the temperature as a function of altitude, $p(h)$ is the atmospheric pressure, and $e(h)$ is the partial pressure of water vapor. These constants are documented in the constants module.

The results obtained using Equation \ref{eq:indexAtmosphere} are illustrated in Figure \ref{fig:atmosphericIndexOfRefraction}. This approximation provides valuable insight into the altitudes where variations in the index of refraction most significantly affect a seeker's performance. As expected, the critical region lies between $20~\mathrm{[km]}$ and $30~\mathrm{[km]}$ above sea level, a range where atmospheric effects play a crucial role in optical tracking. 

\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{atmosphericOptics.png}
    \caption{Atmospheric index of refraction for dry air. \label{fig:atmosphericIndexOfRefraction}}
\end{figure}

The Gladstone-Dale constant is a key parameter used to calculate the index of
refraction as given by Equation \ref{eq:indexGD}. This constant is directly related to the fluid’s density and, in the case of hypersonic flows, must account for the contributions of multiple species due to the high-temperature chemical nonequilibrium effects present in such flows.
 
 For a dilute gas, the index of refraction can be approximated as:

\begin{equation}\label{eq:indexGD}
n - 1 = \sum\limits_{s = 1}^{N} K_s \rho_s
\end{equation}

Where: $n$ is the index of refraction, $K_s$ represents the Gladstone-Dale
constant for species $s$, and $\rho_s$ is the corresponding species density. 

Figure \ref{fig:speciesGladstoneDale} presents the Gladstone-Dale constants for a five-species gas. The species densities were computed using the CFD tool, SU2 [@Maier2021], [@Maier2023a]. These computed densities were then imported into the `HAOT` tool to generate Figure \ref{fig:speciesGladstoneDale}.
\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{3C_speciesGladstoneDale.png}
    \caption{Species Gladstone-Dale constants for a five-species gas.\label{fig:speciesGladstoneDale}}
\end{figure}

An alternative method for calculating the index of refraction involves using polarizability, as expressed in Equation \ref{eq:indexPolarizability}:
\begin{equation}\label{eq:indexPolarizability}
    n - 1 = \frac{1}{2\epsilon_0}\sum\limits_{s = 1}\alpha_s N_s
\end{equation}
Where: $\epsilon_0$ is the dielectric constant in vacuum, $\alpha_s$ is the polarizability constant for species $s$, and $N_s$ represents the number density of species $s$.

Figure \ref{fig:kerlPolarizability} presents results obtained using an extrapolation method (Equation \ref{eq:kerlExtrapolation}) developed by [@Kerl1992], based on the work of [@Hohm1986]. This method provides a temperature-dependent correction to the static polarizability, improving accuracy for high-temperature conditions relevant to hypersonic flows and aero-optical calculations.
\begin{equation}\label{eq:kerlExtrapolation}
\alpha(\omega, T) = \frac{\alpha(0,0)}{1 - \left(\frac{\omega}{\omega_0}\right)^2} \left( 1 + bT + cT^2\right)
\end{equation}
Where: $\alpha(\omega, T)$ represents the polarizability as a function of
laser's frequency $\omega$ and temperature $T$. The parameters $b$ and $c$ are
empirical extrapolation constants, and $\omega_0$ is the characteristic oscillation frequency of the
diatomic molecule. These values are implemented in the `HAOT` package, allowing for easy computation of temperature-dependent polarizability.

Figure \ref{fig:kerlPolarizability} illustrates the variation of polarizability with temperature for oxygen $(O_2)$ at a laser wavelength of $633\, nm$, demonstrating the effectiveness of Kerl’s method in modeling optical properties under varying thermodynamic conditions.
\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{kerlPolarizability_O2_633nm.png}
    \caption{Polarizability of $O2$ using Kerl's method for a laser of
    $633[nm]$.\label{fig:kerlPolarizability}}
\end{figure}

The results presented thus far highlight some of the key capabilities of the aerodynamics and optics modules within `HAOT``. Additionally, Figure \ref{fig:boltzmannDistribution} presents the Boltzmann distribution for nitrogen $(N_2)$, implemented in the quantum mechanics module. Several aero-optics studies [@Buldakov2000], [@Tropina2018] have incorporated Boltzmann distribution calculations to determine the polarizability and index of refraction for diatomic molecules, demonstrating its significance in hypersonic flow modeling.
\begin{figure}[h!]
    \centering
    \includegraphics[width=0.8\textwidth]{boltzmannDistribution_N2.png}
    \caption{Boltzmann Distribution for $N2$.\label{fig:boltzmannDistribution}}
\end{figure}


These tools can be readily applied to the results obtained from a CFD simulation. Below, we present calculations for the index of refraction, Gladstone-Dale constant, and dielectric properties using the `HAOT` package. A low-speed, compressible Large Eddy Simulation (LES) was performed using OpenFOAM to generate the results shown below. The primary objective here is not to provide a detailed description of the CFD setup, but rather to demonstrate the capabilities of the `HAOT` tool in processing and analyzing flow-field optical properties.

\begin{figure}[h!]
    \centering
    \includegraphics[width=1.0\textwidth]{index_of_refraction.png}
    \caption{Computed index of refraction.\label{fig:contourIndex}}
\end{figure}

\begin{figure}[h!]
    \centering
    \includegraphics[width=1.0\textwidth]{kerl_polarizability.png}
    \caption{Computed polarizability.\label{fig:contourPolarizability}}
\end{figure}

\begin{figure}[h!]
    \centering
    \includegraphics[width=1.0\textwidth]{dielectric.png}
    \caption{Computed dielectric medium constant.\label{fig:contourDielectric}}
\end{figure}

# Acknowledgements

The author acknowlege the support from Dr. Kyle Hanquist at the University of
Arizona and Dr.Ozgur Tumuklu at Rensselaer Polytechnic Institue.

# References
