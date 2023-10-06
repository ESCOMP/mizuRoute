# Welcome to mizuRoute webpages

xxxxx

## History

Predecessor to mizuRoute is a river channel routing model coupled with the TopNet hydrologic Model (Clark et al., 2008), which has been used in the National Institute of Water and Atmospheric Research (NIWA) in New Zealand. 
The wave tracking algorithm employed in this routing model was originally developed by Derek Goring in NIWA, and detailed in Goring (1994). 
When Martyn Clark transitioned from NIWA to NCAR in 2010, and embarked on developing a new hydrologic model, SUMMA, the need for a river routing model became apparent. 

Subsequently the wave tracking subroutine in the TopNet model was extracted and underwent refactoring, ultimately enabling it to run over the Continental US with USGS catchment-river Geospatial Fabric data. 
During that time period, US hydrology projection studies (e.g. CMIP3 and CMIP5 hydrology projection), which were our focal point of our research, were utilizing the Lohmann’s Impulse Response Function (IRF) routing (Lohmann et al., 1996) to produce streamflow projections at selected gauges throughout the Continental US. 

This prompted us to implement the same Lohman’s routing scheme, allowing us to emulate both Lohman’s routing and Goring’s wave tracking routing. Although the routing code was capable of running using USGS Geospatial Fabric over multi-decades, it was slow using a single core. 
Consequently, we opted to document the routing development work as a publication in the Geoscientific Model Development (GMD) journal. 
Given the journal’s requirement for the model’s name to be included in the paper’s title, we initiated discussions in our group’s Slack channel. 
At a certain juncture, Andy Wood suggested, "How about “mizuRoute”? since I (my last name can be divided into two parts and the first part is “mizu” ) led the paper and became increasingly involved in the code development. 
As It turns out, “mizu” translates to “water” in Japanese, and the suggestion received anonymous approval. 
Somewhat importantly, It is not an acronym (which can be challenging to recall in its entirety) and is concise enough to remember.

Since then, mizuRoute has undergone several rounds of refactoring to enhance both efficiency and usability. 
Initially, the river network connectivity computation was improved through a traditional pair coding exercise conducted over one weekend by Martyn and me, and subsequently, a hybrid parallelization implementation was introduced. 
This enhancement has empowered us to produce an ensemble of projected flows at numerous critical gauge locations, apply it to the global river network, then couple it to Earth System Model (i.e., Community Earth System Model).


## Documentations

Documentation is available [online](https://mizuroute.readthedocs.io/en/latest/). More scientific information should be referred to [Publications](#publication)


## Code

Code is available at [github](https://github.com/ESCOMP/mizuRoute)

For use for smaller domain, use main branch (only use OpenMP for parallel run)

For use for high resolution or large domain, use cesm-coupling branch (require MPI library)


## Contributions

We value community efforts on scientific work. This includes model development. We welcome any contributions to help mizuRoute to improve further. Please contact xxx.


## Publications(#publication)

Gharari, S., Vanderkelen, I., Tefs, A., Mizukami, N., Lawrence, D., and Clark, M. P.: A Flexible Multi-Scale Framework to Simulate Lakes 1 and Reservoirs in Earth System Models, Water Resour. Res., in review, https://doi.org/10.1002/essoar.10510902.1, 2022.

Cortés-Salazar, N., Vásquez, N., Mizukami, N., Mendoza, P. A., and Vargas, X.: To what extent does river routing matter in hydrological modeling?, Hydrol. Earth Syst. Sci., 27, 3505–3524, https://doi.org/10.5194/hess-27-3505-2023, 2023.

Vanderkelen, I., Gharari, S., Mizukami, N., Clark, M. P., Lawrence, D. M., Swenson, S., Pokhrel, Y., Hanasaki, N., van Griensven, A., and Thiery, W.: Evaluating a reservoir parametrization in the vector-based global routing model mizuRoute (v2.0.1) for Earth system model coupling, Geosci. Model Dev., 15, 4163–4192, https://doi.org/10.5194/gmd-15-4163-2022, 2022.

Li, Z., Gao, S., Chen, M., Gourley, J., Mizukami, N., and Hong, Y.: CREST-VEC: a framework towards more accurate and realistic flood simulation across scales, Geosci. Model Dev., 15, 6181–6196, https://doi.org/10.5194/gmd-15-6181-2022, 2022.

Mizukami, N., Clark, M. P., Gharari, S., Kluzek, E., Pan, M., Lin, P., et al. (2021). A vector-based river routing model for Earth System Models: Parallelization and global applications. Journal of Advances in Modeling Earth Systems, 13, e2020MS002434. https://doi.org/10.1029/2020MS002434.

Mizukami, N., Clark, M. P., Sampson, K., Nijssen, B., Mao, Y., McMillan, H., Viger, R. J., Markstrom, S. L., Hay, L. E., Woods, R., Arnold, J. R., and Brekke, L. D.: mizuRoute version 1: a river network routing tool for a continental domain water resources applications, Geosci. Model Dev., 9, 2223–2238, https://doi.org/10.5194/gmd-9-2223-2016, 2016.


## Authors of the page

Naoki Mizukami
