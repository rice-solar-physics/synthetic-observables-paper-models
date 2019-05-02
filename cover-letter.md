---
# spell-checker: disable
documentclass: scrartcl
title: Cover Letter
subtitle: "*Understanding Heating in Active Region Cores through Machine Learning I. Numerical Modeling and Predicted Observables*"
author:
- W.T. Barnes
- S.J. Bradshaw
- N.M. Viall
bibliography: [paper/references.bib]
# spell-checker: enable
---
We are grateful to the referee for providing a thorough review of the manuscript and for making several very useful suggestions, all of which have greatly improved this manuscript. Below, we respond to each of the referee's specific concerns. Additionally, we note these changes in the text as appropriate using the "track changes" functionality provided by the AAS\TeX\ document class. A summary of these changes is included on the last page of the corrected manuscript as well.

# Specific Concerns

The comments of the referee are noted in *italics*. Our response to each concern is directly below in the indented block of text. The numbering scheme is the same as in the original referee report.

1. *In the introduction, the authors write that "nanoflare is now synonymous with any impulsive energy release and is not specific to any particular physical mechanism." While some authors assume this, it is far from being universally true in the community, and many would take it to be synonymous with magnetic reconnection.*

      > We have removed this sentence altogether. While we do not model any specific physical heating mechanism in this paper, this sentence contradicts our later explanation for why the heating rate should scale with the waiting time between events. See also point 4.

2. *In Section 2, the authors state that they trace 5000 field lines. Usually magnetic field extrapolations do not fix the number this way, so it would be worth explaining why the authors have done so.*

      > We choose to trace 5000 field lines in order to make the emission as volume filling as possible while also balancing computation time. In order to carry out this step, we place seed points at the lower boundary in areas of strong positive polarity. We have added a sentence in the text to make this procedure more clear. The magnetic field extrapolation calculation is independent of the field line tracing step.

3. *The authors have chosen to deposit all of the energy into electrons. Would the results differ in any significant way should they deposit energy into ions instead (or some mixture)?*

      > While likely an important consideration in the exploration of coronal heating, the effect of ion heating on these observables is outside the scope of this paper. However, we investigated the effect of energizing the ions versus the electrons in two previous papers [@barnes_inference_2016;@barnes_inference_2016-1]. We found that preferentially heating the electrons versus the ions affects the presence of "hot" plasma, but does not impact cooling behavior. As the emission measure slope and the time lag primarily probe the cooling and draining phase of the plasma evolution, we do not expect that preferentially heating the ions will significantly change our results. We have added a note to this effect in section 2.2. The one exception to this could be the presence of the negative 131 \AA\ time lags as these are likely a probe of much hotter plasma that is created before the two species have had time to equilibrate.

4. *In this paper, nanoflare is stated to mean heating by any mechanism, but the authors explicitly state that motivation for the heating model is the reconnection scenario in Parker 1988. This is a contradiction.*

      > The contradictory statement in the introduction has been removed. See point 1.

5. *Parameterizing the waiting time between heating events in terms of the cooling time seems odd. Can this be physically justified? Naively, one might expect shorter loops to reconnect more often (e.g. Tiwari et al 2017, APJ, 843, L20).*

      > We choose to define the heating frequency in terms of a cooling time because the cooling time scales approximately with the loop length [see appendix of @cargill_active_2014]. Thus, when we say "low frequency," we mean that every loop, independent of the loop length, will, on average, undergo a complete cooling and draining cycle before being reheated. Under this definition, it is indeed true that shorter loops reconnect more often. Taking the absolute waiting time, $t_\textup{wait}$, as a proxy for the reconnection timescale, we find that, under this definition of the heating frequency, shorter loops do reconnect more often as $t_\textup{wait}\propto L$ for a given value of $\varepsilon$ (see Equation 1). This is also addressed in the second to last paragraph of the discussion section.

6. *The authors have shown the derived slopes of the EM distribution, but it's not clear that the peak temperatures are consistent with observations. Figure 3 suggests that the temperature might peak relatively low (2 MK or so). It would be worthwhile to show a few example distributions, or maybe a map of the peak temperature.*

      > Distribution of peak temperature in the EM, estimating emission from this AR; are we putting in enough energy? check flux out of our AR by multiplying DEM and loss function and integrating over temperature

7. *The scaling law derived at the bottom of page 12 may be dubious for a few reasons. The exponent $\alpha = -1/2$ is not a value that appears in @rosner_dynamics_1978. In fact, @rosner_dynamics_1978 give $\alpha = 0$ over the range $5.75 < \log T < 6.3$, and $\alpha = -2/3$ from $6.3 < \log T < 7.0$, which will obviously give a different scaling law. The second issue is that the scaling law from @cargill_implications_1994 may not be appropriate: @cargill_implications_1994 notes that the approximation is only to zeroth order, assuming static radiative cooling, whereas the @bradshaw_cooling_2010 scaling explicitly includes coronal draining. This appears incompatible.*

      > We have removed the reference to @rosner_dynamics_1978. The reference for $\alpha=-1/2$ should be @cargill_implications_1994. The inconsistency in the scaling law formation is a valid concern. We have corrected the scaling law to explicitly use $\ell=2$ and make a note that consideration of enthalpy-driven cooling is likely to alter this scaling relation, particularly for longer loops.

8. *In Figure 7 and accompanying text, it is difficult to evaluate how realistic the time-lag maps are compared to the observations by @viall_survey_2017. The colorbar and spatial scale differ from that paper. It would be worthwhile to use the same scales to facilitate comparison.*

      > The limits on the colorbar in Figure 7 have been changed to $\pm6000$ s to be consistent with @viall_survey_2017. While the field-of-view of our simulated images is slightly different compared to that of @viall_survey_2017, this is not critical as we do not wish to make detailed comparisons between the simulation and observation here as this will be covered in much more detail in Paper II. See also point 11.

9. *The authors claim on page 14 that the 94-335 time lags are consistent with @viall_survey_2017. It's not clear that any of the simulated cases are similar to the observations. @viall_survey_2017 have time lags up to their limit of -6000 s, and a significant amount of 0 lag. The simulations here do not reproduce either of these features. Similarly, the authors should comment on the discrepancies between the simulated and observed 193-171 pair.*

      > Our only intention here is to say that in our "cooling" heating scenario, we find negative time lags of similar magnitude to @viall_survey_2017 in the longest loops in the periphery of the AR. We have clarified the language in this paragraph to make our intentions more explicit. Additionally, we note that observed 94-335 time lags have far fewer positive time lags than our simulations. We have also added a note that there are far more zero time lags in the observations than in our simulations. This discrepancy in the number of zero time lags is addressed in the last paragraph of this same section. We do not attempt to make detailed qualitative or quantitative comparisons between observed and simulated results here. This is reserved for Paper 2 (see point 11).

10. *The third conclusion in the summary is contradicted by Figure 9. Examining 94-335, 94-131 or 335-211, for example, there does not appear to be any trend towards a uniform distribution between high and low frequencies.*

      > We have corrected the third point of our conclusions to only state that the distribution becomes increasingly broad with increasing frequency as the tendency toward a uniform distribution is not demonstrated here. We mean to convey that in the limit of steady heating ($t_\textup{wait}\to0$), variations in the intensity will be dominated by noise which leads to uniform distribution in the time lags [@viall_signatures_2016]. We have added a clarification to this effect in the discussion section.

11. *Can the authors provide any conclusions about the heating frequency? The observed time lag and EM maps appear closest to the high frequency case, though there are certainly differences*

      > The comparison of our model to observations in the context of the heating frequency is a very valid concern and one that we will address in the second paper in this series which is in its final stages of preparation. The primary purpose of this paper (Paper I) is to explain our forward modeling code and show our simulated emission measure slopes and time lags. In Paper II, we will show the time lags and emission measure slopes calculated from observations of this same active region. We will then use machine learning, specifically a random forest classifier, in order to quantitatively evaluate which heating frequency matches the observation "best" given all 15 time lags, 15 cross-correlations, and the emission measure slope. We expect to submit this paper to The Astrophysical Journal within a week of the submission of the Paper 1 revisions.

# Minor Concerns

* *"Nanoflares" is misspelled on page 2, first paragraph.*

    > This has been corrected. Several other minor spelling mistakes have also been corrected.

* *Version 9 of CHIANTI has been released. If the authors are so inclined, it might be worth redoing the calculations as a check, particularly since there are changes to ionization equilibria for iron.*

    > While we appreciate the referee's attention to detail regarding the importance of atomic data in this study, redoing these calculations with CHIANTI v9 would essentially mean repeating the entire study from the very beginning. Furthermore, we suspect that the major changes in this latest version of the database (i.e. the inclusion of dielectronic recombination and autoionization in the level population calculation) are not likely to significantly affect our results. In light of this and the significant effort that this change would require, we will opt not to incorporate CHIANTI v9 in this paper.

\section{References}
