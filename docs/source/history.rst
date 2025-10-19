.. _rst_History:

History
==============

Predecessor to mizuRoute is a river channel routing model coupled with the TopNet hydrologic Model (:ref:`Clark et al., 2008 <Clark2008>`), which has been used in the National Institute of Water and Atmospheric Research (NIWA) in New Zealand.
The wave tracking algorithm employed in this routing model was originally developed by Derek Goring in NIWA, and detailed in :ref:`Goring (1994) <Goring1994>`.
When Martyn Clark transitioned from NIWA to NCAR in 2010, and embarked on developing a new hydrologic model, SUMMA, the need for a river routing model became apparent.

Subsequently the wave tracking subroutine in the TopNet model was extracted and underwent refactoring, ultimately enabling it to run over the Continental US with USGS catchment-river Geospatial Fabric data.
During that time period, US hydrology projection studies (e.g. CMIP3 and CMIP5 hydrology projection), which were our focal point of our research, were utilizing the Lohmann’s Impulse Response Function (IRF) routing (Lohmann et al., 1996) to produce streamflow projections at selected gauges throughout the Continental US.

This prompted us to implement the same Lohman’s routing scheme, allowing us to emulate both Lohman’s routing and Goring’s wave tracking routing.
Although the routing code was capable of running using USGS Geospatial Fabric over multi-decades, it was slow using a single core. Consequently, we opted to document the routing development work as a publication in the Geoscientific Model Development (GMD) journal.
Given the journal’s requirement for the model’s name to be included in the paper’s title, we initiated discussions in our group’s Slack channel.
At a certain juncture, Andy Wood suggested, “How about “mizuRoute”? since I (my last name can be divided into two parts and the first part is “mizu”, 水 in kanji) led the paper and became increasingly involved in the code development.
As It turns out, “mizu” translates to “water” in Japanese, and the suggestion received anonymous approval.
Somewhat importantly, It is not an acronym (which can be challenging to recall in its entirety) and is concise enough to remember.

Since then, mizuRoute has undergone several rounds of refactoring to enhance both efficiency and usability.
Initially, the river network connectivity computation was improved through a traditional pair coding exercise conducted over one weekend by Martyn and me, and subsequently, a hybrid parallelization implementation was introduced.
This enhancement has empowered us to produce an ensemble of projected flows at numerous critical gauge locations, apply it to the global river network, then couple it to Earth System Model (i.e., Community Earth System Model).

At the same time, lakes are added in river-only network, and new capability of simulating lake water balance was added by Shervan Gharari (University of Saskatchewan), and Inne vanderkelen at VUB worked with Shervan to implement paramertric reservoir release models.


**References**

.. _Clark2008:

Clark, M.P., Rupp, D.E., Woods, R.A., Zheng, X., Ibbitt, R.P., Slater, A.G., Schmidt, J. and Uddstrom, M.J. (2008).
Hydrological data assimilation with the ensemble Kalman filter: Use of streamflow observations to update states in a distributed hydrological model.
Advances in water resources, 31(10), pp.1309-1324.
https://doi.org/10.1016/j.advwatres.2008.06.005

.. _Goring1994:

Goring, D. G. (1994)
Kinematic shocks and monoclinal waves in the Waimakariri, a steep, braided, gravel-bed river,
Proceedings of the International Symposium on waves: Physical and numerical modelling, University of British Columbia, Vancouver, Canada, 336–345.
