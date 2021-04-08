/**
\mainpage What is PlasmaC?

<table bgcolor="#FFFFFF" align="top">
<tr>
<td valign="top">
PlasmaC is a scalable computer code for two- and three-dimensional low-temperature fluid plasmas in complex geometries. We currently support

  * Multi-material Poisson solvers
  * Diffusive radiate transport
  * Convection-diffusion reaction systems

PlasmaC is build on top of Chombo, and therefore features

  * Cut-cell representation of geometries
  * Multi-material geometries (electrodes and insulators)
  * Patch-based adaptive mesh refinement
  * Excellent scalability to hundres of thousands of cores

Please see the [gallery] (@ref gallery) for a selection of results

To get started, go [here](@ref doxy-contents)

You can obtain the PlasmaC from here:

      git clone ssh://git@git.code.sintef.no/~robertm/chombo-streamer

If you do not have MFChombo installed, it can be pulled from here:

      git clone ssh://git@git.code.sintef.no/~robertm/mf-chombo
</td>
<td valign="top">
\image html spacecharge_black.png
\image html mechshaft_charge_red_white.png
</td>
</table>
