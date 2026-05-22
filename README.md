[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.595402.svg)](https://doi.org/10.5281/zenodo.595402)
[![Documentation Status](https://readthedocs.org/projects/mizuroute/badge/?version=main)](https://mizuroute.readthedocs.io/en/main/?badge=main)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

# mizuRoute

## Overview

MizuRoute is a community river model for simulating streamflow through user-defined, catchment-based vector river network. It takes runoff or lateral inflow from from a hydrologic, land surface model, or other data sources and routes water through rivers, lakes, and reservoirs.
MizuRoute is designed for multiple applications, including stand-alone hydrologic and water-resources studies, regional and continental river-routing experiments, model intercomparison studies, and coupled Earth system modeling. 
MizuRoute supports parallel computing with Passing Interface (MPI), enabling multi-decdal to a century simulations over large and fine scale river networks such as regional, continental or global domains.
We are also intereseted in expanding modeling capability, for example, tracers, streamtemperature , data assimilations etc. as optional features, in addition to streamflow estimations.

## Applications and use cases

mizuRoute can be used in several ways:

2. **Hydrologic and water-resources modeling**  
   mizuRoute can be applied to regional or continental river networks for hydrologic studies, practical water-resource planning such as climate-impact assessments, and hydrologic model evaluation via observed streamflow comparison.

3. **Coupled Earth system modeling**  
   mizuRoute can be used as a river routing component in Earth System Model to study feedbacks from river to terrestial hydrologic system and discharge into oscean. 

4. **Method development and community research**  
   mizuRoute is intended to support development and testing of routing methods, lake/reservoir representations, network datasets, coupling strategies, and hydraulic or hydrologic process improvements.

## Getting started

mizuRoute is commonly used in two configurations:

1. **Stand-alone model**  
   Users provide runoff or lateral inflow time series and a river network to simulate streamflow, river storage, lake/reservoir volume, and related routing diagnostics. This is how majorty of users use.
   Please refer to [getting started stand-alone mode](docs/starting_stand-alone_mode.md) 

2. **ESM/land model coupled mode** 
   MizuRoute runs as a river model component in Community Earth System Model (CESM). 
   Currently, mizuRoute is coupled to Community Terrestrial Systems Model ([CTSM](https://github.com/ESCOMP/CTSM)), the land-model component of Community Earth System Model(CESM). 
   Please refer to [getting started CTSM coupled mode](docs/starting_ctsm_couple_mode.md) 
   
Most new users interested in rive routing applications should begin with [getting started stand-alone mode](docs/starting_stand-alone_mode.md) and [mizuRoute documentation](http://mizuRoute.readthedocs.io/)
Users interested in coupled Earth system simulations should consult the CESM/CTSM documentation in addition to the mizuRoute documentation.

## Community model and contributions

mizuRoute is a multi-application community river routing model. It supports stand-alone hydrologic and water-resources applications as well as coupled Earth system modeling applications such as CESM/CTSM. 
The project welcomes contributions that improve model physics, numerical methods, input/output workflows, testing, documentation, examples, and usability across these applications.

For small changes, such as documentation updates or minor bug fixes, contributors may open pull request directly. For larger changes, such as new routing schemes, changes to input/output formats, data-structure refactoring, or modifications that may affect CESM/CTSM coupling, contributors are encouraged to open a GitHub issue first. This allows maintainers and interested users to discuss the proposed augmentation, identify possible impacts on existing workflows, coordinate with ongoing development, and agree on a review strategy before implementation.

A suggested contribution process is:

1. Open an issue describing the proposed change, motivation, expected application, and likely impacts.
2. Discuss the proposed design with maintainers and interested community members.
3. Develop the change in a feature branch in your fork.
4. Add or update tests, examples, and documentation.
5. Open a pull request for review.

Current development interests include but are not limited to, improving support for stand-alone and coupled applications, refining routing methods, improving lake/reservoir and lateral inflow representations, expanding tests and examples, improving documentation, and modernizing the code structure to better support multiple routing options.

## Prerequisite 
To effectively work with mizuRoute, users and developers are desired to have the following working commands. 
 1. Linux commands.
 2. Geographic Information System to develop and visualize river network data
 3. python or similar other languages to analyze/visualize data, prepare for the input
 4. Fortran (if a user desires to change the codes)


## Funding
 - U.S. Army Corps of Engineers
 - Bureau of Reclamation
 - NASA’s Advanced Information Systems Technology program
 - National Science Foundation


LICENSE
-------
mizuRoute is [licensed](LICENSE.txt) under Apache 2.0
