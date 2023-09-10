## Welcome to mizuRoute webpages 

This website provides informal information on mizuRoute, and link to technical  

### Descriptions

xxxxxx

### History

Predecessor of mizuRoute is a river channel routing model coupled with TopNet hydrologic Model (Clark et al., 2008), used in National Institute of Water and Atmospheric Research (NIWA) in New Zealand. 
The wave tracking algorithm used in the routing model was originally developed by Derek Goring in NIWA, and described in Goring (1994). When Martyn Clark moved to NCAR from NIWA in 2010, and focused more on developing new hydrologic model, SUMMMA, 
a need for a river routing model was realized. Since then a routing model component in the TopNet model was extracted and went through refactoring, and eventually enabled it to run over Continental US domain with USGS' catchment-river geospatial fabric data. 
Meanwhile, US hydrolgy projection researches (e.g., CMIP3, CMIP5 etc.), which our research focuses on, have been used Lohmann's Impulse Response Funtion (IRF) routing (reference) to produce at streamflow projections at selected gauges across the CONUS. 
This motivated us to implement the same Lohman's routing scheme so we could mimic Lohman's routing as well as Goring's wave tracking routing. The routing code was able to run using USGS Geospatial Fabric over multi-decades first and then we decided to document the routing development work as a publication in Geoscientific Model Development journal.
Since GMD journal requires the name of the model to be included, we discussed it in slack. 
At somepoint, Andy Woond suggested "mizuRoute", because I (my last name include "mizu" as short name) was leading the paper and somewhow woring on the code more and more. It turns out that "mizu" is "water" in Japanese, and everyone liked it. It is not acronyum, which can be hard to remember (for especially non-English speaker).
 
Since then mizuRoute went through a few refactoring for efficiency and usability; first, river network connectivity computation, and hybrid parallelization implementation. This allows us to produce ensemble of projected flows at many key gauge locations, as well as run over global domain.


```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```
