---
title: 'chombo-discharge: An AMR code for gas discharge simulations in complex geometries'
tags:
  - Python
  - astronomy
  - dynamics
  - galactic dynamics
  - milky way
authors:
  - name: Robert Marskar
    orcid: 0000-0003-1706-9736
    affiliation: 1
affiliations:
 - name: SINTEF Energy Research, Norway
   index: 1
date: 5 January 2023
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

chombo-discharge is a parallelized Cartesian 2D and 3D adaptive code for simulating low-temperature gas discharges in complex geometries. 
Such discharges occur when electrons accelerate in strong electric fields and ionize the gas, and further evolution is affected by residual space charges.
Streamers, for example, are filamentary plasma dominated by space charge effects.
They are the natural precursors to leader, sparks, and lightning.

Gas discharge modeling involves simulations over multiple scales in time and space.
chombo-discharge reduces the cost of such simulations by using Cartesian Adaptive Mesh Refinement (AMR).
It also provides support for multi-material complex geometries (gas phase, electrodes, and solid dielectrics) through an embedded boundary (EB) formulation.
Geometries are represented as implicit functions, and can be created using constructive solid geometry.
Conversion of surface meshes to implicit functions is also supported.
Under the hood, chombo-discharge uses Chombo [@ebchombo] for the AMR and EB infrastructure, and is parallelized using MPI.
However, chombo-discharge supplies all numerical solvers.

chombo-discharge uses a solver-centered modular design where larger applications are developed by coupling numerical solvers in the chombo-discharge base code, using C++ interfaces.
Many solvers already exist in chombo-discharge, all of which are parallelized and compatible with EBs and AMR:

* Advection-diffusion-reaction solvers.
* Helmholtz equation solvers, using geometric multigrid. 
* An electrostatic solver (with support for discontinuous coefficients).
* Kinetic Monte Carlo chemistry solvers. 
* Radiative transfer solvers (continuum and Monte Carlo).
* Various particle solvers, e.g. for Monte Carlo radiative transfer, tracer particles, and microscopic drift-diffusion.
* ODE solvers defined over volume or surface meshes.

All solvers exist as stand-alone applications, and many of them are also coupled through more complex physics applications that aim at resolving different types of discharge phenomena (e.g. statistical inception models, or particle and fluid models of streamer discharges).
The interaction of these solvers occurs through a common AMR "core", which can also use dual grids where e.g. fluid and particle kernels are load-balanced separately.
Depending on their needs, users can therefore enter the framework at several levels. 
E.g. they need to learn interfaces when using existing applications (e.g., streamer models); use C++ APIs if developing new physics applications, or use the EB-AMR infrastructure if contributing with new solvers.

# Statement of need

There is already a number of discharge simulation codes currently available.
Commercial codes used for simulating discharges include COMSOL (see e.g. @mcplas), ANSYS Fluent (see e.g. @Niknezhad_2021), and PLASIMO [@plasimo].
Another example is Afivo-streamer [@Teunissen:2017], which is also open source and uses Cartesian AMR.

While chombo-discharge is not the only open-source discharge simulation code, it has a number of unique features.
In particular, chombo-discharge supports complex geometries, and is therefore useful in many discharge-related science applications.
Support for AMR is also important, as AMR is a virtual requirement in many 3D applications (certainly the ones involving filamentary plasma).
The code is also quite performant, and its design pattern permits a flexible and extensible approach to numerically solving various discharge-related problems, even when these end up requiring many thousands of CPU cores.
Originally, the code was written for studying pre-breakdown discharges in high-voltage equipment [@Marskar:2019], but over time it has been adapted in order to fit a wider category of discharge-related problems.
Two science examples are given in \autoref{fig:surface} (HV technology) and \autoref{fig:corona} (pulsed discharge).
Lightning initiation investigations (from hydrometeors), plasma medicine, and plasma-assisted combustion are other examples where chombo-discharge could potentially be used.

![2D Surface discharges over complex surfaces [@Meyer:2022]. Top: Electrode (shaded region) and dielectric (profiled surface). Bottom left panel: Electric field magnitude. Bottom right panel: Plasma density. \label{fig:surface}](figures/SquareEvolution.pdf)

![Streamer discharge tree simulation in full 3D using Kinetic Monte Carlo. The two figures show the same discharge, viewed from the side and from the bottom. \label{fig:corona}](figures/ItoKMC.pdf)

# Acknowledgements

The development of chombo-discharge was partially achieved through funding by the Research Council of Norway through grants 245422, 319930/E20, and 321449.

# References