---
title: 'chombo-discharge: An adaptive code for gas discharge simulations in complex geometries'
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

chombo-discharge is a 2D and 3D adaptive code for simulating low-temperature gas discharges in complex geometries. 
Such discharges occur when electrons accelerate in strong electric fields and ionize the gas, and further evolution is also determined by residual space charge.
Streamers, for example, are filamentary plasma dominated by space charge effects.
They are the natural precursors to leader, sparks, and lightning.

Gas discharge modeling involves simulations over multiple scales in time and space.
chombo-discharge reduces the cost of such simulations by using Cartesian Adaptive Mesh Refinement (AMR).
It also provides support for multi-material complex geometries (gas phase, electrodes, and solid dielectrics) through an embedded boundary (EB) formulation.
Geometries can be constructed using either constructive solid geometry (CSG) or imported as surface meshes. 
The code uses a solver-centered modular design where larger applications are developed by strongly or weakly coupling numerical solvers shipped with the chombo-discharge base code, using existing C++ interfaces.
Several such applications are included in the chombo-discharge code. 
Depending on their needs, users can therefore enter the framework at several levels. 
E.g. they only need to learn textual interfaces interfaces when using existing applications; use solver C++ APIs when developing new applications, or use the EB AMR infrastructure when contributing entirely new solvers.

Several solvers exist in chombo-discharge, all of which are parallelized and compatible with EBs and AMR:

* Advection-diffusion-reaction solvers.
* Helmholtz equation solvers, using geometric multigrid. 
* Electrostatic solvers (with support for discontinuous coefficients).
* Kinetic Monte Carlo solvers. 
* Radiative transfer solvers (continuum and Monte Carlo solvers).
* Various particle solvers, e.g. for Monte Carlo radiative transfer, tracer particles, Brownian walkers, or kinetic Particle-In-Cell.
* Volumetric and cut-cell ODE solvers.

All solvers exist as stand-alone applications, and many of them are also coupled through more complex physics applications that aim at resolving different types of discharge phenomena (e.g. statistical inception models or streamer discharges).
The interaction of these solvers occurs through a common AMR core, which can also use dual grids where e.g. fluid and particle kernels are load-balanced separately.
chombo-discharge uses the Chombo infrastructure for the AMR and EB infrastructure, and is parallelized using MPI.
However, chombo-discharge supplies most solver and infrastructure discretizations.

# Statement of need

There is already a number of discharge simulation codes currently available.
Commercial codes used for simulating discharges include COMSOL, ANSYS Fluent, PLASIMO, VSIM/VPIC, and Vizglow.
Several non-commerical codes also exist, such as Afivo-streamer [@Teunissen:2017], which is also open source and uses Cartesian AMR.

While chombo-discharge is not the only open-source discharge simulation code, it has a number of unique features.
In particular, chombo-discharge supports complex geometries, and is therefore useful in many discharge-related science applications.
Lightning initiation investigations (from hydrometeors), high-voltage technology, plasma medicine, and plasma-assisted combustion are typical examples.
The code is also quite performant, and its design pattern permits a flexible and extensible approach to numerically solving various discharge problems, even when these require many thousands of CPU cores.

Originally, the code was written for studying pre-breakdown phenomena in high-voltage equipment [@Marskar:2019], but over time it has been adapted to fit a wider category of discharge problems.
Some examples science examples are given in \autoref{fig:surface} (HV technology), \autoref{fig:corona} (nanosecond pulsed discharges), and \autoref{fig:sprite} (upper atmospheric lightning). 

![2D Surface discharges over complex surfaces [@Meyer:2022]. Top: Electrode (shaded region) and dielectric (profiled surface). Bottom left panel: Electric field distribution. Bottom right panel: Plasma density. \label{fig:surface}](figures/SquareEvolution.pdf)

![Streamer corona in a nanosecond pulsed discharge in atmospheric air. The two figures show the same discharge, viewed from the side and from the bottom. \label{fig:corona}](figures/ItoKMC.pdf)

![Sprite discharges in the upper atmosphere.\label{fig:sprite}](figures/sprite3d_3.png)

# Acknowledgements

The development of chombo-discharge was partially achieved through funding by the Research Council of Norway through grants 245422, 319930/E20, and 321449.

# References