# Welcome to mizuRoute webpages

xxxxx

## History

Predecessor of mizuRoute is a river channel routing model coupled with TopNet hydrologic Model (Clark et al., 2008), which is used in National Institute of Water and Atmospheric Research (NIWA) in New Zealand. 
The wave tracking algorithm used in this routing model was originally developed by Derek Goring in NIWA, and described in Goring (1994). 
When Martyn Clark moved to NCAR from NIWA in 2010, and began developing a new hydrologic model, SUMMA, a need for a river routing model was realized. 

Then the wave tracking subroutine in the TopNet model was extracted and went through refactoring, and eventually enabled it to run over Continental US with USGS catchment-river Geospatial Fabric data. 
Meanwhile, US hydrology projection researches (e.g. CMIP3, CMIP5 hydrology etc.), which our research focuses on, was using the Lohmann’s Impulse Response Function (IRF) routing (Lohmann et al., 1996) to produce streamflow projections at selected gauges across the Continental US. 

This motivated us to implement the same Lohman’s routing scheme so we could mimic Lohman’s routing as well as Goring’s wave tracking routing. The routing code was able to run using USGS Geospatial Fabric over multi-decades, though quite slow using a single core. 
Then we decided to document the routing development work as a publication in the Geoscientific Model Development (GMD) journal. Since GMD journal requires the name of the model in the paper title, we started discussing it in the group slack channel. At some point, Andy Wood suggested, "How about “mizuRoute”? because I (my last name can break into two parts and the first part is “mizu” ) was leading the paper and also started working on the code more and more. It turns out that “mizu” is “water” in Japanese, and everyone liked it. It is not an acronym (which can be hard to remember its full name) and short enough.

Since then mizuRoute went through several refactoring for efficiency and usability; first, river network connectivity computation by traditional pair coding exercise done one weekend, and hybrid parallelization implementation. 
This allows us to produce an ensemble of projected flows at many key gauge locations, as well as run it over the global river network.


## Contributions

We value community efforts on scientific work. This includes model development. We welcome any contributions to help mizuRoute to improve further. Please contact xxx.


## Publication

Gharari, S., Vanderkelen, I., Tefs, A., Mizukami, N., Lawrence, D., and Clark, M. P.: A Flexible Multi-Scale Framework to Simulate Lakes 1 and Reservoirs in Earth System Models, Water Resour. Res., in review, https://doi.org/10.1002/essoar.10510902.1, 2022.

Cortés-Salazar, N., Vásquez, N., Mizukami, N., Mendoza, P., and Vargas, X.: To what extent does river routing matter in hydrological modeling?, Hydrol. Earth Syst. Sci. Discuss. [preprint], https://doi.org/10.5194/hess-2022-338, accepted, 2023.

Vanderkelen, I., Gharari, S., Mizukami, N., Clark, M. P., Lawrence, D. M., Swenson, S., Pokhrel, Y., Hanasaki, N., van Griensven, A., and Thiery, W.: Evaluating a reservoir parametrization in the vector-based global routing model mizuRoute (v2.0.1) for Earth system model coupling, Geosci. Model Dev., 15, 4163–4192, https://doi.org/10.5194/gmd-15-4163-2022, 2022.

Li, Z., Gao, S., Chen, M., Gourley, J., Mizukami, N., and Hong, Y.: CREST-VEC: a framework towards more accurate and realistic flood simulation across scales, Geosci. Model Dev., 15, 6181–6196, https://doi.org/10.5194/gmd-15-6181-2022, 2022.

Mizukami, N., Clark, M. P., Gharari, S., Kluzek, E., Pan, M., Lin, P., et al. (2021). A vector-based river routing model for Earth System Models: Parallelization and global applications. Journal of Advances in Modeling Earth Systems, 13, e2020MS002434. https://doi.org/10.1029/2020MS002434.

Mizukami, N., Clark, M. P., Sampson, K., Nijssen, B., Mao, Y., McMillan, H., Viger, R. J., Markstrom, S. L., Hay, L. E., Woods, R., Arnold, J. R., and Brekke, L. D.: mizuRoute version 1: a river network routing tool for a continental domain water resources applications, Geosci. Model Dev., 9, 2223–2238, https://doi.org/10.5194/gmd-9-2223-2016, 2016.

