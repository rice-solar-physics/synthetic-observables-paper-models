# Paper Ideas

Here are some ideas for this paper and other papers to follow based on this project of computing timelags from simulated AIA intensities over a whole AR.

## Parameters survey of nanoflare heating frequencies

* Use EBTEL to simulate many loops in AR NOAA 1158
* Survey several different heating models
  * Low frequency
  * Intermediate frequency
  * High frequency
  * Randomly occurring nanoflares (maximally decoherent)
  * Single pulse, cooling only (maximally coherent)
* Compute AIA intensities in all 6 channels
* Compute timelags for all 15 channel pairs
* Compute emission measure slopes
* Compute these observables from the observed intensities as well.
* Ask two primary questions:
  * What are the observable signatures of these different heating frequencies?
  * How do the modeled observables compare to the observations?
* Compare modeled timelags between heating types
* "By-eye" comparison between models and observations
* Train random forest classifier on model data and apply to observations
  * How do results change when we add the emission measure slope?
  * How do the timelag only results compare to the emission measure slopes?

### Color Palette

* For our three different heating functions, use the first three colors in the `colorblind` palette in Seaborn
* For different elements in the response function plot, use colorbrewer `Set1`

## Bundle Heating Model

* Bundle heating model that connects heating frequency to field strength
* Still constrained by AR flux
* Run AR model using HYDRAD, ~1000-2000 loops
* Look at both electron and ion heating (?)
* Viability of this type of mechanism
* May help to constrain how the heating is connected to the field
* Compare EM and timelag observables with observations
* Companion to Paper 1

## Investigation of how viewing angle affects observables

* Shorter more succinct paper
* HYDRAD simulation of dipole AR or use AR from previous study
* Compute timelags and emission measure slopes for a single heating frequency
* Look at multiple viewing angles
  * Disk center
  * Slightly off center
  * Off limb
  * Arcade
* Look at distributions of observables for varying LOS angles
* Are they different? How?
* Would we come to different conclusions about the heating for a different LOS?

## Classifying timelags in a large survey of active regions

* Look at all 15 ARs from Warren et al. (2011)
* Run AR simulations with HYDRAD for ~1000 loops
* Survey ~3 heating frequencies and compute timelags in all 15 pairs
* Compute timelags in all 15 ARs
* Train random forest classifier on model timelags
* Classify observed timelags in all ARs
* Survey of heating frequencies across many ARs
* Systematic assessment of the frequency in the context of the our models, accounting for all dimensions of our results
* Good testbed for analysis of AIA data at scale

## JOSS paper describing the synthesizAR Python package

* describe code
* minimal scientific results
* peer review
* gives a way for others to cite the code, a DOI!
* review of the software, not the results