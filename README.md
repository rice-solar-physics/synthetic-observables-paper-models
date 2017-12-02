# synthesized_timelag_maps
Synthesize AIA Emission from several active regions and Compute Timelag Maps for Channel Pairs

The main set of ARs we will be looking at come from [Warren et al. (2011)](http://iopscience.iop.org/article/10.1088/0004-637X/759/2/141/meta)
In particular, we'll focus on 
* NOAA 1158

Note that right now, we are only looking at NOAA 1158. If we expand to more ARs, we can create subdirectories
for each active region. For each distinct AR geometry, we'll build a base field that all of the different 
heating configurations will inherit.

For each AR, we'll run a series of four heating models:
* cooling only
* high frequency nanoflares
* intermediate frequency nanoflares
* low frequency nanoflares

