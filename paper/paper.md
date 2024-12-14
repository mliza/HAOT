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

Hypersonic flows present a unique challenges due to the complex interplay of
fluid dynamics, chemical reactions, and optical phenomena. As a signal from a
Light Detection and Ranging (LiDAR) travels through a hypersonic flow field,
the beam would be affected by the flow. 

`HAOT` is a Hypersonic Aerodynamics Optics Tools Python package developed to
calculate the index of refraction of a hypersonic medium. Its source code is
available on [GitHub](https://github.com/mliza/HAOT), the documentation is
available on
[Read the Docs](https://haot.readthedocs.io/en/latest/) and example on the
usage of the package is given on the GitHub repo under the example folder. 

# Statement of Need
Many techniques used to calculate optical properties are scatter in papers but
there is not a local repo containing all this calculations, furthermore some of
these calculations require the use of spectroscopy constants, which have been
properly documented and added to the package.

# Algorithms
The `HAOT` package, contains five  modules:

- Modules:
    - Aerodynamics
    - Optics
    - Quantum Mechanics
    - Constants
    - Conversions

Each module can be imported independed
and the documentation explains in detail what each module does. Furthremore,
docstrings have been added to the function and the description of each function
can be seeing in an interactive python session. 

Equation below was introduced by [@Smith1953], and it is a good approximation
for the change of the index of refraction as a function of altitude. 

$$ n(h) \approx 1 + \frac{K_1}{T(h)} \left( p(h) + K_2\frac{e(h)}{T(h)} \right) \label{eq:atmosphericIndex} $$

Where: $K_1$ and $K_2$ are constants, $T$ is the temperature as a function of
altitude, $p$ is pressure as a function of altitude, and $e(h)$ is partial pressure
of water vapor.

![Atmospheric index of refraction for dry air.\label{fig:atmosphericIndexOfRefraction}](atmosphericOptics.png)

Equation below shows the equation used to calculate the dilute index of refraction. 

$$ n -1 = \rho \sum\limits_{s =1}^N K_s \rho_s \label{eq:diluteIndexOfRefraction}$$
Where:  $\rho_s$
is the species density, $\rho$ is the flow's density, and $K_s$ is the specie's
Gladstone-Dale constant.

![Species Gladstone-Dale constants for a 5 species gas.\label{fig:indexOfRefraction5Species}](3C_speciesGladstoneDale.png)

![Index of Refraction for a 5 species gas.\label{fig:indexOfRefraction5Species}](3C_refractionIndex.png)


and referenced from text using \autoref{fig:diluteIndexOfRefraction}.





A more extensive work showing the results of this pacakge was done by [@Liza2023].

# Acknowledgements

# References
