---
title: "Bayesian MR Dissertation Notes"
output: html_document
date: "2024-09-29"
---


# TwoSampleMR




# Literature Review

Mendelian randomisation (MR) is a methodology intended to support causal inference from observational data. It applies the principles of instrumental variable (IV) analysis to genetic data: naturally occurring genetic variants - “instruments” – with a known association to an exposure of interest are chosen, and by comparing the association of those same genetic instruments to an outcome of interest, genetic data can be used to investigate causal links between exposures and outcomes. In theory, provided that the assumptions of IV analysis are met, random assignment of genetic variants from parents to offspring during meiosis can create a form of natural experiment, analogous to randomisation during a clinical trial – both measured and unmeasured confounders should be distributed evenly between the groups created, allowing valid causal inference after other sources of bias and random variation are accounted for.

Three key assumptions of IV analysis must be met:
1.	Relevance – the genetic variant must be associated with the exposure of interest
2.	Independence – the genetic variant is independent of confounders of the relationship between exposure and outcome
3.	Exclusion restriction – the genetic variant must not be associated with the outcome except via the exposure
If these assumptions are satisfied, the “causal effect” of the exposure on the outcome can be estimated by the Wald ratio, i.e. by dividing the co-efficient of gene-outcome association by the co-efficient of gene-exposure, giving a numerical measure of strength of causal exposure-outcome association.
 
Figure X. Taken from Burgess et al 2016 (DAG)




In practice, only the relevance assumption can be directly tested and proven, as independence and exclusion restriction depend on all possible confounders of the exposure-outcome association, both measured and unmeasured. Threats to the independence assumption will vary depending on the population, exposure and outcome being studied. Exclusion restriction is a particularly universal issue in MR, due to so-called (horizontal) genetic pleiotropy, where a single genetic variant may have multiple “pleiotropic” effects – i.e. it may influence several traits simultaneously. Such pleiotropic effects may be unknown and open unmeasured causal pathways between a genetic instrument and the outcome, separate to the path involving the exposure of interest, thus potentially biasing MR estimates of the association between exposure and outcome.2

Although not possible to prove exclusion restriction for any MR study, several methods attempt to produce exposure-outcome causal effect estimates which are robust to violations of this assumption. A common approach is the Weighted Median Estimator (WME) method, proposed by Bowden et al 2. 

[ In WME analysis, several genetic instruments are used to estimate the exposure-outcome causal effect; each instrument is known to be associated with the exposure of interest, but a proportion of these instruments may be invalid due to unknown pleiotropic genetic effects. Any genetic instrument causally linked an outcome via multiple pleiotropic causal pathways would be expected to exhibit a less consistent gene-outcome association than if only a single pathway mediated the gene-outcome relationship; this would be reflected in larger variance in causal effect estimates derived from invalid/pleiotropic genetic instruments. WME therefore assigns a weight to each genetic instrument’s estimate of the causal effect according to the inverse of the variance of the estimate; these weighted effect estimates are used to construct a cumulative distribution function for probability of true causal effect size across the range of estimated values.

Causal effect estimates from each instrument are ordered by size, then used to create a cumulative distribution function for probability of true causal effect size. The relative contribution of each instrument’s effect estimate to the probability distribution is weighted according to the inverse of the variance of the estimate. Genetic instruments whose causal effect estimates exhibit a large variance, which would be expected would therefore contribute less to ]

In WME analysis, several genetic instruments are used to estimate the exposure-outcome causal effect. Any instrument linked to an outcome via multiple pleiotropic causal pathways will exhibit a less consistent gene-outcome association than a relationship mediated by a single pathway; this results in larger variance in causal estimates derived from invalid/pleiotropic genetic instruments. WME therefore assigns a weight to each genetic instrument’s causal estimate according to the inverse of its variance, then constructs a cumulative distribution function for probability of true causal effect size across the range of estimated values, before taking the 50th percentile of this distribution as a “weighted median estimate” of the true causal effect, theoretically producing consistent causal estimates even if up to 50% of the included information comes from invalid instruments.

# Methods

Precision and Consistency of Causal Estimates – Simulation

The precision and consistency of causal estimates will first be quantified for each methodology using simulated datasets where characteristics of genetic instruments are known, including pleiotropy of effects. 

Data will be generated from models and parameters based on those previously published in the original exposition of the weighted median estimator method 2. Three scenarios will be examined in turn:

1.	Balanced pleiotropy, InSIDE assumption satisfied
2.	Directional (positive) pleiotropy, InSIDE assumption satisfied
3.	Directional (positive) pleiotropy via a confounder, InSIDE assumption not satisfied

For each set of assumptions, 25 candidate genetic instruments will be simulated, and causal estimates will be calculated using both weighted median and MR-Hevo methods. Furthermore, the effects on the performance of each method in the two-sample MR setting will be examined for differing combinations of the following variables across 10,000 simulated datasets:

•	Proportion of invalid genetic instruments: this will be simulated in turn at 0.1, 0.2 and 0.3 for each scenario considered.

•	Number of participants: 10,000 versus 20,000

•	Value of causal effect: null (β = 0) and positive (β = 0.1)

For each scenario and combination of variables, 10,000 datasets will be generated and used to calculate the mean causal estimate, mean standard error of the causal estimate, and the ??power of the 95% confidence interval to reject the null hypothesis??(/positive report rate/similar)?? generated by both the weighted median and MR-Hevo methods. These estimates will initially be compared by tabulation. 

In addition, estimates will be graphically compared by plotting the coefficients of regression of outcome on each instrument against coefficients of regression of exposure of each instrument. Plotting genetic associations with outcome versus genetic association with exposure in this way for both WME and MR-Hevo allows comparison of their causal effect estimates, represented as the gradient of a line of best fit through the points generated by each method.
 
Effect on Conclusions Drawn – Re-Analysis of Published Data

The potential impact of MR-Hevo casual effect estimation on real-world conclusions will be investigated through re-analysis of published data. Studies for re-analysis will be identified which have used WME to make causal inferences in a two-sample MR setting. 

To determine the upper-bound of impact that MR-Hevo methods could plausibly have on causal effect estimation versus WME in the field of MR, the search for studies to re-analyse will prioritise studies which have been most highly cited. Highly cited studies demonstrably have directly influenced the largest numbers of other academic works; it is also likely that they will also have the greatest indirect influences on the field, via the downstream impact on subsequent studies performed. In addition, highly cited studies have arguably been exposed to the greatest cumulative scrutiny, and would be expected to be free of significant methodological flaws which might otherwise impact validity of study conclusions. If even one highly-cited MR study were found to have differing conclusions after re-analysis using the MR-Hevo method, then this could represent a meaningful impact on the field.

On this basis, potential MR studies for re-analysis will first be identified through a forward citation search, using the Scopus platform (***?REF***) to identify papers referencing the original WME exposition paper. Returned studies will then be ranked by the number of times they themselves have been cited. The top 10 most cited papers meeting following criteria will be re-analysed:
•	Inclusion criteria:
o	Reporting an original, two-sample MR study
o	Study populations reported sufficiently to determine ancestry and presence/degree of participants overlapping between studies
o	Reporting on ≥20 human genetic instruments relating to a defined exposure
o	Genetic instruments reported sufficiently to allow re-analysis, i.e. to include:
•	Details of effect/non-effect alleles
•	Regression coefficient estimates and standard errors for both exposure and outcome available as reported in manuscript OR accessible via MR-Base
o	Uses Weighted Median Estimator causal estimate methods
•	Exclusion criteria: 
o	Methodology paper, review article, editorial or letter only; does not report original two-sample MR study as the primary focus
o	English full-text not accessible

The above sample of 10 highly-cited studies will then be re-analysed using MR-Hevo in place of WME, and the resulting conclusion regarding the main study question will be recorded. The proportion of studies where conclusions of main research questions differ between MR-Hevo and MWE analyses will then be quantified, and the potential reasons for and implications of such differences will be discussed

