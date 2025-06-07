---
title: "Causal Effect Estimation in Mendelian Randomisation Studies - Evaluating a Novel Bayesian Approach To Genetic Pleiotropy Versus Established Weighted Median Methodology"
author: "B233241"
date: "September 2024 - June 2025"
bibliography: 
  - Lit_Rev_Refs.bib
  - grateful-refs.bib
csl: vancouver-superscript.csl
link-citations: TRUE
urlcolor: blue
acronyms:
  loa_title: ""
  include_unused: FALSE
  insert_loa: FALSE
  insert_links: TRUE
  id_prefix: "acronyms_"
  sorting: "alphabetical"
  non_existing: "key"
  style: "long-short"
  fromfile: ./_acronyms.yml
---

\newpage

## Acknowledgements {-}

I would like to acknowledge


## Contributions {-}

Mine others

## Statement of Originality {-}

I confirm that all work is my own except where indicated, that all sources are clearly referenced....


## Word Count {-}

Word count: 
6647

\newpage

<!--chapter:end:index.Rmd-->

# Introduction and Background




## Introduction to Mendelian Randomisation (MR)

Epidemiology is the study of determinants and distribution of disease across populations; a common epidemiological study aim is therefore to seek evidence as to whether a given exposure (e.g. cigarette smoking) may cause a given outcome (e.g. lung cancer) [@coggon_chapter_2003]. Logistics limit experimental interventions across large groups, so insights into associations between exposures and outcomes are gleaned from observational data of people in the population of interest. Comparing health outcomes between individuals with different levels of a particular exposure may highlight potential links, e.g. higher cancer incidence in those who smoke more is consistent with a causal role for cigarettes in carcinogenesis [@coggon_chapter_2003]. 

However, correlation does not prove causation. A key epidemiological challenge is accounting for so-called "confounding" factors; these are other variables, associated with both the exposure and the outcome of interest, which represent an alternative causal explanation for any exposure-outcome links observed[@martens_instrumental_2006]. If smokers also drink more alcohol than non-smokers, then an observed link between smoking and increased cancer risk could plausibly be caused by increased alcohol exposure, either partially or entirely. Another potential issue with observational data is "reverse causation", where the presumed outcome is in fact a cause of the exposure; this might be the case if a cancer diagnosis drove individuals to drink and smoke more, and data were collected without respect to exposure timings.

\acr{MR} is a methodology intended to support causal inference from observational data. It applies the principles of \acr{IV} analysis to genetic data, performing a type of natural experiment often likened to a \acr{RCT} [@hernan_instruments_2006]. 

In a properly conducted \acr{RCT}, causality can be inferred due to a randomisation process being used as an "instrument" to allocate different levels of exposures to different experimental groups. If groups are randomly allocated, any confounding variables which might otherwise influence exposure-outcome relationships should be evenly distributed between groups, whether these confounders are known or not. As such, there should be no systematic differences between individuals from different groups in the exposure of interest - that is, there should be no bias [@stel_instrumental_2013]. Statistical methods can quantify the probability that any observed outcome differences could have occurred by chance, and thereafter any outcome differences can be interpreted as caused by exposure differences. As allocation and receipt of exposures is known to precede outcome measurements, reverse causality is impossible.

In \acr{MR}, naturally occurring genetic variants - “genetic instruments” – are chosen based on their known association to an exposure of interest. Provided that assumptions of \acr{IV} analysis are met, random assignment of alleles (i.e. variants of a given gene) from parents to offspring during meiosis creates randomisation analagous to that performed for an \acr{RCT} – both measured and unmeasured confounders should be distributed evenly between the groups created, allowing valid causal inference after other sources of bias and random variation are accounted for [@davies_reading_2018].

## Causal Effect Estimation in MR

At its simplest, the relationship between two continuous variables - an exposure $X$ and outcome $Y$ - can be represented as a linear model:

\begin{equation} 
Y = \alpha + \beta X + \epsilon
\end{equation}

where $\alpha$ represents all non-$X$ determinants of $Y$, $\beta$ is the causal effect of $X$ on $Y$ and $\epsilon$ is an error term. The $\beta$ term is a numerical measure of strength of causal exposure-outcome association, where:

  - $\beta = 0$ implies no causal link between exposure and outcome
  - $\beta > 0$ implies $X$ causes $Y$
  - $\beta < 0$ implies $X$ prevents $Y$ 
  
To estimate a causal effect using a genetic variant in an \acr{IV} analysis, three key assumptions must be met [@lousdal_introduction_2018]:

1.	Relevance – the genetic variant must be associated with the exposure of interest
2.	Independence – the genetic variant is independent of confounders of the relationship between exposure and outcome
3.	Exclusion restriction – the genetic variant must not be associated with the outcome except via the exposure

These assumptions are represented graphically in Figure \@ref(fig:DAG-assumptions-plot). 





![(\#fig:DAG-assumptions-plot)Causal diagram illustrating the relationships between genetic instrument _G_, exposure _X_, outcome _Y_ and confounders of the exposure-outcome relationship _U_ in Mendelian randomisation studies. Blue text & crosses represent key assumptions to ensure valid inference of causal effect of _X_ on _Y_ using _G_ as an instrumental variable. Red text represents violations of these assumptions that may lead to invalid inference through opening of alternate causal pathways. Greek characters represent the key parameters/association coefficients to be estimated. Adapted from Burgess et al 2016[@burgess_sensitivity_2016]](2_Intro_Background_files/figure-latex/DAG-assumptions-plot-1.pdf) 


Typically, \acr{MR} studies estimate causal effect using a set of several genetic instruments; the causal effect estimate derived from the $jth$ instrument is denoted $\hat{\beta}_j$. Each estimate $\hat{\beta}_j$ acknowledges there will be specific effects on the observed values of exposure and outcome given the presence of that specific genetic variable $G_j$ under study, i.e. $\hat{\beta}_j$ is based on the instrument-conditioned exposure ${X|G_j}$ and outcome ${Y|G_j}$. These observed values of exposure and outcome can be described by their own linear models:

\begin{equation} 
X|G_j = \gamma_0 + \gamma_j G_j + \epsilon_{X_j}
\end{equation}

\begin{equation} 
Y|G_j = \Gamma_0 + \Gamma_j G_j + \epsilon_{Y_j}
\end{equation}

where, for exposure and outcome respectively:

  - $\gamma_0$ and $\Gamma_0$ reflect base values without influence of the genetic variant
  - $\gamma_j$ and $\Gamma_j$ are coefficients of association with the genetic variant, representing the extent to which an effect allele of $G_j$ will perturb the value of $X$ or $Y$ versus the non-effect allele
  - $\epsilon_{X_j}$ and $\epsilon_{Y_j}$ are error terms, containing contributions from confounders of the exposure-outcome relationship ($U$ in the causal diagram), and all genetic variants except $G_j$.

It can be shown that a simple causal effect estimate for the exposure on the outcome can be obtained from a single genetic instrument by the Wald method, dividing the coefficient of gene-outcome association by the coefficient of gene-exposure association, i.e.:

\begin{equation} 
\hat{\beta}_j = \frac {\hat{\Gamma}_j} {\hat{\gamma}_j}
\end{equation}

Each instrument may be valid or invalid, depending on it meeting the above assumptions. The overall causal effect estimate $\hat{\beta}$ from any given \acr{MR} method will typically seek to pool effect estimates from several instruments so as to minimise effects of any invalid instruments included, e.g. by removing/down-weighting contributions of genetic instruments which violate one or more assumptions. This is equivalent to plotting all estimated coefficients of gene-outcome association ($\bar{\Gamma}$) versus all estimated coefficients of gene-exposure association ($\bar{\gamma}$) for the set of instruments, then using the gradient of a regression line through the points as the causal effect estimate $\hat{\beta}$; picking an \acr{MR} methodology is analogous to choosing the method to draw the line of best fit (Figure \@ref(fig:Gamma-gamma-plot)). For binary outcomes, the causal effect estimate can be converted to an odds ratio (OR) through exponentiation, i.e.:

\begin{equation} 
OR = e^{\hat{\beta}}
\end{equation}

![(\#fig:Gamma-gamma-plot)Simulated MR Study on 10,000 individuals using 25 genetic instruments, of which 30% are invalid (red points) and introduce directional pleiotropic effects. The true value of the exposure-outcome causal effect is 0.25 (grey line, causal effect represented by gradient). Regression using an unajusted least-squares linear model (light blue line) results in a biased estimate in the positive direction due to the influence of the invalid instruments. Using the Weighted Median Estimator method (pink line) attenuates the effects of the invalid instruments, resulting in an estimate closer to the true value. Adapted from Bowden et al 2016[@bowden_consistent_2016]](2_Intro_Background_files/figure-latex/Gamma-gamma-plot-1.pdf) 

## Violations to Assumptions

In practice, only the relevance assumption can be directly tested and proven. Typically, genetic variants for \acr{MR} studies are selected as instruments based on Genome Wide Association Studies (GWAS), which quantify associations between small genetic variations - known as \acr{SNP}s - and various phenotypes. Association between genetic variants and a phenotypes representing exposures of interest can be partly assured by selection using an appropriate genome-wide significance level (e.g. $p < 10 ^{-8}$). Statistical testing can also quantify the gene-exposure relationship; commonly used measures include the $r^2$ statistic, representing the proportion of variance in the exposure explained by the genotype, and the related $F$-statistic, which additionally accounts for the the sample size under investigation [@richmond_mendelian_2022]. An $F$-statistic of $\ge$ 10 is generally considered to represent a strong enough gene-exposure association to consider a genetic instrument for use [@martens_instrumental_2006].

The assumptions of independence and exclusion restriction depend on all possible confounders of the exposure-outcome association, both measured and unmeasured; as such, these can never be proven absolutely. Various methods have been proposed to quantify and account for violations of these two additional assumptions, including the weighted median estimator, described below [@bowden_consistent_2016].
 
The main methods to avoid violations of the independence assumption relate to appropriate selection of populations studied to avoid confounding due to ancestry or population stratification. For example, in two-sample MR studies, where gene-exposure and gene-outcome coefficients are estimated from two separate GWAS studies, it is recommended to select GWAS studies performed in similar population groups (e.g. both in Western Europeans). This practice helps avoid spurious exposure-outcome associations being generated by confounding due to underlying differences in allele frequency, baseline disease risks etc between ancestrally different populations [@richmond_mendelian_2022].

Exclusion restriction is a particularly universal issue in MR, due to so-called (horizontal) genetic pleiotropy, where a single genetic variant may have multiple “pleiotropic” effects – i.e. it may influence several traits simultaneously. Such pleiotropic effects may be unknown and open unmeasured causal pathways between a genetic instrument and the outcome (Figure \@ref(fig:DAG-assumptions-plot)), thus potentially biasing \acr{MR} estimates of the association between exposure and outcome. As pleiotropy influences outcome separate to the path involving the exposure of interest, the term "direct effects" is also used[@hemani_evaluating_2018]. Where pleiotropic effects are in both positive and negative directions with a mean of zero - "balanced pleiotropy" - then they only add noise to causal effect estimation [@morrison_mendelian_2020]. By contrast, "directional pleiotropy", where the mean of pleiotropic effects is non-zero, may introduce bias [@bowden_consistent_2016].

If such an additional causal pathway acts between gene $G$ and outcome $Y$ via a confounding factor $U$, then the magnitude of direct/overall effects of $G$ on $Y$ will correlate with the effects of $G$ on $X$ (i.e. $\Gamma \propto \gamma$), and "correlated pleiotropy" is present. If an additional causal pathway acts directly between gene $G$ and outcome $Y$ independent of both exposure $X$ and confounders $U$, this results in "uncorrelated pleiotropy" (Figure \@ref(fig:DAG-assumptions-plot)). Both correlated and uncorrelated pleiotropy can introduce bias which distorts the estimate of the true causal effect. In general, correlated pleiotropy is more challenging to account for; several MR methods explicitly require an additional assumption of \acr{InSIDE}, i.e no correlated pleiotropy to be present [@grant_bayesian_2024].

## Weighted Median Estimator (WME)

A common approach to produce exposure-outcome causal effect estimates robust to violations of the exclusion restriction assumption is the \acr{WME}  method, proposed by Bowden et al [@bowden_consistent_2016]. 

In \acr{WME} analysis, several genetic instruments are used to estimate the exposure-outcome causal effect $\hat{\beta}$. Each instrument is known to be associated with the exposure of interest, but an unknown proportion of these instruments may be invalid due to pleiotropic genetic effects. Any instrument linked to an outcome via multiple pleiotropic causal pathways will exhibit a less consistent gene-outcome association than a relationship mediated by a single pathway; this results in larger variance in causal estimates derived from invalid/pleiotropic genetic instruments versus estimates from valid instruments.

\acr{WME} therefore assigns a weight to each genetic instrument’s estimate of the causal effect according to the inverse of the variance of the estimate; these weighted effect estimates are used to construct a cumulative distribution function for probability of true causal effect size across the range of estimated values. The 50th percentile of this distribution can then be taken as a “weighted median estimate” of the true causal effect, theoretically producing consistent causal estimates even if up to 50% of the included information comes from invalid instruments [@bowden_consistent_2016]. An example of \acr{WME} attenuating the effects of invalid instruments is shown in Figure \@ref(fig:Gamma-gamma-plot).

## Issues With WME CIs

\acr{WME} calculation methods are available via several prolific \acr{MR} tools: the R packages “MendelianRandomization”[@yavorska_mendelianrandomization_2017] and “TwoSampleMR”, and the MR-Base web platform[@hemani_mr-base_2018]. However, these implement the original authors’ suggested process of generating 95% confidence intervals for \acr{WME}, which deviates from accepted re-sampling methodology:

>“We found the bootstrap confidence interval…too conservative. However, the bootstrap standard error… gave more reasonable coverage using a normal approximation (estimate ±1.96 x standard error) to form a 95% confidence interval”[@bowden_consistent_2016]

This modification, explicitly aiming to boost estimate precision artificially, would be expected to lead to a high Type 1 error rate, which has been a growing concern in the field of late[@stender_reclaiming_2024]. The theoretical issues with this approach, and the fundamentals of bootstrapping in general, are covered in Appendix \@(ref:appendix-boot).

## MR-Hevo

MR-Hevo is an R package which uses more typical Bayesian methodology to estimate MR causal effects and corresponding 95% confidence intervals. It uses the probabilistic programming language, Stan, to directly sample the posterior probability distribution of pleiotropic effects on the outcome, rather than making untested assumptions about the shape of this distribution as current \acr{WME} implementations do[@mckeigue_inference_2024].

MR-Hevo incorporates several additional features which its creators claim further aid valid causal inference. Most MR methods can only account for one genetic variant per genetic locus (i.e. per location in the genome). If multiple variants exist at a given locus, generally only one can be selected as an instrument for further \acr{MR} analysis. HR-Hevo handles multiple instruments per genetic locus via scalar construction, essentially assigning a "score" to each locus based on the variant(s) present, thus incorporating more information than if closely grouped variants had been discarded [@mckeigue_inference_2024]. As MR-Hevo is based on a Bayesian approach, it generates estimates which incorporate relevant existing information generated by prior studies, increasing the amount of data informing each estimate. In this case, MR-Hevo specifies a prior probability distribution reflecting prior knowledge that most individual genetic instruments will have only small effects on complex traits[@park_estimation_2010; @piironen_sparsity_2017], further aiding biologically plausible inference regarding distribution of pleiotropic effects.

## Aims and Objectives

The main aim of this study will be to demonstrate if the \acr{WME} approach gives over-confident causal estimates in the presence of pleiotropy, and whether this issue is more correctly handled by the MR-Hevo Bayesian approach. This will be achieved through addressing the research questions and objectives as outlined below:


### Research Questions:

1. How does MR-Hevo perform versus the weighted median estimator when estimating causal effects in MR studies?
2. Do conclusions of existing MR studies using weighted median causal effect estimation change if MR-Hevo methods are used?

### Objectives:

1. Quantify the precision of MR-Hevo causal estimates for simulated data under differing sets of common assumptions, with reference to the weighted median estimator
2. Evaluate the consistency of MR-Hevo causal estimates for simulated data under differing sets of common assumptions, with reference to the weighted median estimator
3. Compare the conclusions drawn from MR-Hevo causal effect estimation versus the weighted median estimator on real-world data



\newpage

<!--chapter:end:2_Intro_Background.Rmd-->

# Methods



## Simulation Study

To evaluate the performance of MR-Hevo causal estimation relative to WME, the precision and consistency of both methods were quantified using simulated datasets with known parameter values.

### Data Simulation

To aid comparability with existing methods and literature, the simulation methodology of the original WME exposition was reproduced based on published models and parameters in Appendix 3 of its supplementary materials[@bowden_consistent_2016]. Full details of simulation reproduction, including code and validation of outputs, is presented in \@ref(appendix-sim). 

In brief, simulations were created based on three different scenarios, each representing a common set of assumptions about underlying data used for MR, and each increasingly challenging to the performance of any given MR causal estimation methodology:

1. Balanced pleiotropy, InSIDE assumption satisfied - A proportion of invalid genetic instruments are present and introduce pleiotropic effects uncorrelated with the instrument strength; these pleiotropic effects are equally likely to be positive as negative with a mean value = 0, thus introducing noise into the estimation of causal effect.

2. Directional pleiotropy, InSIDE assumption satisfied - A proportion of invalid genetic instruments are present and introduce pleiotropic effects uncorrelated with the instrument strength; these pleiotropic effects are positive only, with a mean value > 0, thus biasing the causal effect estimate in a positive direction.

3. Directional pleiotropy, InSIDE assumption not satisfied - A proportion of invalid genetic instruments are present and introduce pleiotropic effects correlated with the instrument strength through action via a confounder; these pleiotropic effects are positive only, with a mean value > 0, thus potentially biasing the causal effect estimate in a positive direction to an even greater extent than Scenario 2.

1,000 simulated datasets of participant-level data were generated for every combination of each scenario and each the following simulation parameters:

  - Proportion of invalid instruments: 0%, 10%, 20% or 30%
  - Number of participants: $n = 10,000$ or $n = 20,000$
  - Causal effect: null ($\beta = 0$) or positive ($\beta = 0.1$)
  
  
The same set of 25 simulated genetic instruments were used across all datasets, with the status of each as valid/invalid determined by random draw per instrument at the start of each simulation run of 1,000 datasets.

Genotypes were simulated as for a two-sample setting: where number of particpants was $n = 10,000$, 20,000 genotypes were simulated - 10,000 for the cohort used to estimate gene-exposure association ($\hat{\gamma}$), and a separate cohort of 10,000 used to estimate gene-outcome association ($\hat{\Gamma}$). Parameter values for effect allele frequency were not specified by Bowden et al, though initial testing showed values around 0.5 produced WME causal effect estimates closest to published values when other parameters were matched[@bowden_consistent_2016]. As such, effect allele frequencies were assigned per instrument from a uniform distribution between 0.4 to 0.6. Each effect allele freqency thus generated per instrument was then used as a probability to assign each simulated participant effect alleles for each instrument via two draws from a binomial distribution. 

### Analysis of Simulated Data

Each dataset generated was analysed using both WME and MR-Hevo methods, via functions from the `TwoSampleMR` and `mrhevo` packages, respectively[@hemani_mr-base_2018; @mckeigue_inference_2024]. Results were aggregated per group of 1,000 simulated datasets corresponding to a particular combination of scenario and parameter values. This resulted in one meta-analysis reported per combination of scenario/parameter values, each including 1,000 simulated MR studies using the same 25 genetic instruments in the same population. Aggregated measures for both WME and MR-Hevo per meta-analysis were mean causal effect estimate; mean standard error of the causal effect estimate; and causality report rate, i.e. percentage of simulated studies reported as showing a non-null causal effect, either by p-value <0.05 (WME), or by a 95% credible interval for causal effect estimate not including 0 (MR-Hevo).

Results of the above aggregations were tabulated as per Tables 2 and 3 of Bowden et al[@bowden_consistent_2016] to allow direct comparisons of both methods versus each other and versus the published characteristics of existing MR causal estimation methods.

## Re-Analysis of Published Data

To investigate the potential implications of any differences in performance between WME and MR-Hevo methods, a selection of published studies resporting causal effect estimates using the WME method was re-analysed. A sample size of 10 published studies was decided as a pragmatic compromise between the scope of this study and the need to check consistency of any observed differences. In the original Bowden et al simulation studies, the WME causal estimation method was  shown to generate a false-positive report rate of $\ge$ 30% with relatively minor violations of relevant assumptions[@bowden_consistent_2016]; therefore, even this relatively small sample of 10 studies might be expected to demonstrate differences between methods if the MR-Hevo approach is as appropriately conservative as its creators propose.

To estimate the upper bound of the potential impact of MR-Hevo versus existing WME methodology, studies were chosen for re-analysis based on their number of citations in the wider MR literature. Compared to studies with few or no citations, highly-cited studies would be expected to have a larger impact on their respective fields if their conclusions were to change. In addition, highly-cited works will typically have been submitted to more scrutiny than less-cited works - both during peer review whilst under consideration by journals likely to produce highly-cited works, and from the wider scientific community following the widespread dissemination evidenced by a high citation count. As such, it would be expected that highly-cited works are likely to be free of significant methodological flaws which may impede interpretation of any re-analysis.

### Citation Search

The Scopus search platform [@] was used on 15/04/2025 to retrieve all articles citing the original weighted median estimator exposition paper [@bowden_consistent_2016]. The articles returned were sorted by the number of times each article itself had been cited, and the resulting list was saved to RIS format in blocks of ten articles for upload into the Covidence evidence synthesis platform. Abstracts were screened by a single reviewer (B233241), starting with the most cited article and proceeding in descending order of citation count, against the following inclusion and exclusion criteria:

Inclusion criteria:

- Original two-sample MR study

- Able to determine samples’ ancestry sufficient to establish presence/potential degree of participant overlap between groups

- Reporting $\ge$ 20 human genetic instruments relating to exposure

- Reports details of effect/non-effect alleles

- Regression coefficients and standard errors and/or confidence intervals available for each genetic instrument used
  
- Uses Weighted Median Estimator

Exclusion criteria: 

- Methodology paper, review article, editorial or letter

- English full-text not accessible


Where eligibility could not be determined from abstract screening alone, full texts were retrieved and screened against the same criteria. Screening of abstracts and full texts was undertaken in blocks of ten articles, until the target of ten included studies for reanalysis had been reached. 

Where an article reported multiple exposure-outcome associations, data were only extracted for the association with the highest number of genetic instruments available, or else for the first reported association where several were based on the same number of instruments. Data were extracted from full texts of included studies using a standardised data collection template, which included publication details, citation count, primary study question, degree of participant overlap between groups, number/details of genetic instruments used, effect estimates/standard errors calculated, and conclusion regarding causality as determined by the weighted median estimator method. 

## Data Manipulation and Analysis

All simulations, data manipulations and data analyses were performed in R version 4.4.3 (2025-02-28 ucrt)[@base]. 

For the simulation study, full details of computation are available in Appendix \@ref(appendix-sim).

For citation search data, a standardised data collection form was Microsoft Excel [@microsoft_corporation_microsoft_2018] to create .csv files for subsequent analysis in R; Excel's "Get Data" function was also used to extract tables of genetic instruments where these were presented in non-csv format (e.g. pdf). 

Data cleaning for citation search data was primarily undertaken using the Tidyverse suite of R packages [@tidyverse]. A full list of packages used can be found in Appendix \@(ref:appendix-pkg). 

Data were manually screened at summary level and relevant features were extracted. Data were checked for completeness, consistency, duplicate values and plausibility. Data were transformed to an appropriate data type, and encoding of genetic variables was standardised into a single format. Missing values for association coefficients and \acr{SE}s were imputed as the mean value calculated per dataset. It was noted during early testing that MR-Hevo functions do not operate correctly when zero values are present in coefficients of genetic association or their standard errors; such zero values were therefore re-coded as an arbitrarily low value of $10^{-100}$.

## Ethical Approval

The protocol for this work has been reviewed and approved by the \acr{UMREG} at the University of Edinburgh, Ethics ID: UM241126. Due to the nature of the project, using simulated and publically available data only, no significant ethical issues were foreseen, and sponsorship was deemed unnecessary by the \acr{UMREG} reviewing panel.



<!--chapter:end:3_Methods.Rmd-->

# References

<div id="refs"></div>

<!-- export Zotero as Western encoding -->

\newpage

<!--chapter:end:8_References.Rmd-->

# (APPENDIX) Appendix {-} 

\newpage


# Appendix: List of Abbreviations {#appendix-acr}


\printacronyms


\newpage

# Appendix: Bootstrapping {#appendix-boot}

## Bootstrapping - General Method

The typical process for "bootstrap" generating an estimate, \acr{SE} and \acr{CI}s of a population parameter (e.g. population mean $\mu$) from a sample $x$ is as follows[@buscaglia_chapter_2020]:

1. A sample, $x$, of $n$ individuals is selected from a total population, $X$, of $N$ individuals
2. This sample $x$ is then treated as the "bootstrap population"; the empirical distribution of values in the $n$ individuals in the bootstrap population is taken to be broadly representative of the distribution of values in the underlying population $X$ of $N$ individuals
3. A "bootstrap sample", $x^*$,  is then obtained by re-sampling individuals from the bootstrap population with replacement $n$ times per bootstrap sample, i.e. the new bootstrap sample comprises $n$ sampled individuals, $x^*_1, x^*_2,...x^*_n$. As such, individuals from the original bootstrap population $x$ may contribute once, more than once or not at all to each bootstrap sample $x^*$.
4. A total of $k$ bootstrap samples are generated, $x^{*1}, x^{*2},...x^{*k}$, and the statistic of interest (e.g. sample mean $\bar{x}$) is estimated in each individual sample, $\bar{x}^{*i}$, giving the complete set of $\bar{x}^{*1}, \bar{x}^{*2},...\bar{x}^{*i}...\bar{x}^{*k}$.
5. The set of $k$ statistics are combined to form a "bootstrap distribution"; as expected from \acr{CLT}[@ross_chapter_2014], this is typically closer to a normal distribution than the underlying distribution of values in either the bootstrap population $x$ or the total population $X$. (See Figure \@ref(fig:prostate-vol) for an example of this)
6. The final values are derived as follows:

> - the parameter estimate (e.g. estimate of the true population mean, $\hat{\mu}$) is taken as the mean of the bootstrap distribution of $k$ estimates, $(\sum^k_{i = 1} \bar{x}) \div n$
> - the \acr{CI}s are taken as the values at the appropriate centiles at the edges of the sampling distribution, e.g. a 95% \acr{CI} would be generated using values at the 2.5th and 97.5th centiles
> - the \acr{SE} of the estimate is taken as the \acr{SD} of the sampling distribution, given by $\sqrt{\frac{1}{k - 1} \sum^k_{i = 1} (\bar{x_i} - \hat{\mu})^2}$



## Bootstrapping - Example: Prostate Volume

The above process is illustrated in \@ref(fig:prostate-vol). Data on prostate volume in 307 prostate cancer patients demonstrates a right-skewed distribution (A). An empirical distribution from a sample of 100 of these patients mirrors this right skew, and is used as the "bootstrap population" (B) for further re-sampling. As the bootstrap population is re-sampled more and more times, the "bootstrap distribution" of the sample means generated (C and D) gradually tends towards a normal distribution. The 95% \acr{CI} is given by the bounds defining the middle 95% of the bootstrap distribution of estimated means, as shown.

![(\#fig:prostate-vol)Histograms demonstrating distribution of prostate volumes in patients with prostatic cancer, taken from Cata et al 2011[@cata_blood_2011] via the R package `medicaldata`[@medicaldata]. A) Distribution from whole study population of 307 patients with non-missing data, exhibiting right-skew. B) Distribution from random sample of 100 patients, still exhibiting right-skew. C) Bootstrap distribution generated by re-sampling 1,000 bootstrap samples from the original sample of 100 patients, right-skew less apparent. D) Bootstrap distribution generated by re-sampling 100,000 bootstrap samples from the original sample of 100 patients, approaching normality. 95% confidence intervals are demonstrated in plots C and D by marking the 2.5th and 97.5th centiles.](9_Appendices_files/figure-latex/prostate-vol-1.pdf) 

\newpage

## Bootstrapping - Relevance to WME

In current implementations of \acr{WME}, the \acr{WME} estimate of the causal effect ($\hat{\beta}_{WME}$) is calculated as described in Bowden et al[@bowden_consistent_2016], and the 95% \acr{CI} is generated separately using bootstrapping, though notably not using the method described above. 

The bootstrapping process begins similarly, with re-sampling undertaken (a default of $k = 1000$ times) to generate $k$ bootstrap samples $x^{*1}, x^{*2},...x^{*k}$. Each individual bootstrap sample $x^{*i}$ is used to estimate the causal effect using the \acr{WME} method $\hat{\beta}^{*i}_{WME}$, and thus a bootstrap distribution of $k$ values of \acr{WME} is created, $\hat{\beta}^{*1}_{WME},  \hat{\beta}^{*2}_{WME}....\hat{\beta}^{*k}_{WME}$. 

At this stage, however, the bootstrap distribution is then assumed to be approximately normally distributed without verifying this assumption. The 95% \acr{CI} of the bootstrap estimate is then calculated as 1.96 \acr{SD}s of the bootstrap distribution either side of the mean estimate, i.e. $\hat{\beta}_{WME} \pm 1.96  \times SE$. This approach may be problematic for several reasons. 

Although \acr{CLT} leads us to expect that the bootstrap distribution will approach normality as the number of bootstrap iterations $k$ increases, the extent to which this occurs for a given $k$ may depend on the inital distribution of values in the population $X$, and so also on the distribution in the sample/bootstrap population $x$. If the true distribution of values is very non-normal, as may be the case for traits determined by complex genetic and environmental influences, it may take relatively more bootstrap iterations for the bootstrap distribution to become sufficiently normal to assume mean and \acr{SD} accurately describe it. 

Additionally, the bootstrap \acr{SE} is inversely proportional to the number of bootstrap iterations $k$, as opposed to the usual standard error (given by $SE = \frac{SD}{\sqrt{n}}$), which is inversely proportional to the square root of the sample size $n$. It is therefore possible to generate smaller \acr{SE}s by increasing the number of bootstrap samples obtained. This may lead to false confidence in estimates generated despite potential issues with initial sample $x$, e.g. if it too small, or sampled in such a way that it is not representative of the underlying population $X$. Although such issues are inherent to any bootstrapping approaches, the usual method of generating bootstrapped \acr{CI}s detailed above uses more information (i.e. using the entire bootstrap distribution) to generate these values than the parameter-based $estimate \pm 1.96 \times SE$ method (i.e. using approximate summary statistics to represent the distribution). The usual method of bootstrap \acr{CI} generation may therefore be expected to highlight any variation or uncertainty present more readily than the parameter-based approach; this would be represented as wider \acr{CI}s.

\newpage

# Appendix: Simulation Code {#appendix-sim}


## Generating Data and Models {#appendix-sim-gen}

The data generating model used was from Appendix 3 of Bowden et al [@bowden_consistent_2016]; the relevant section describing their model is reproduced below:

>_"..._

>\begin{equation} 
U_i = \sum^J_{j=1} \phi_jG_{ij} + \epsilon_i^U
\end{equation}


>\begin{equation} 
X_i = \sum^J_{j=1} \gamma_jG_{ij} + U_i + \epsilon_i^X
\end{equation}

>\begin{equation} 
Y_i = \sum^J_{j=1} \alpha_jG_{ij} + \beta X_i + U_i + \epsilon_i^Y
\end{equation}

>_for participants indexed by $i = 1, . . . , N$, and genetic instruments indexed by $j = 1, . . . , J$._

>_The error terms $\epsilon_i^U , \epsilon_i^X$ and $\epsilon_i^Y$ were each drawn independently from standard normal distributions. The genetic effects on the exposure γj are drawn from a uniform distribution between 0.03 and 0.1. Pleiotropic effects $\alpha_j$ and $\phi_j$ were set to zero if the genetic instrument was a valid instrumental variable. Otherwise (with probability 0.1, 0.2, or 0.3):_

>_1. In Scenario 1 (balanced pleiotropy, InSIDE satisfied), the $\alpha_j$ parameter was drawn from a uniform distribution between −0.2 and 0.2._

>_2. In Scenario 2 (directional pleiotropy, InSIDE satisfied), the $\alpha_j$ parameter was drawn from a uniform distribution between 0 and 0.2._ 

>_3. In Scenario 3 (directional pleiotropy, InSIDE not satisfied), the $\phi_j$ parameter was drawn from a uniform distribution between −0.2 and 0.2._


>_The causal effect of the exposure on the outcome was either $\beta X = 0$ (null causal effect) or $\beta X = 0.1$ (positive causal effect). A total of 10 000 simulated datasets were generated for sample sizes of N = 10 000 and 20 [sic] participants. Only the summary data, that is genetic associations with the exposure and with the outcome and their standard errors as estimated by univariate regression on the genetic instruments in turn, were used by the analysis methods. In the two-sample setting, data were generated on 2N participants, and genetic associations with the exposure were estimated in the first N participants, and genetic associations with the outcome in the second N participants."_ [@bowden_consistent_2016]

To reproduce this model, code was written in R to generate the relevant participant level data. First, a function (`get_simulated_MR_data`) was written which included parameters specified by Bowden et al, and also to allow testing of data simulation:





This initial simulation function generated data in the following format:




A function (`get_models`) was then written to create linear models from each dataset generated as per Bowden et al:





These models generated estimates of the coefficient of gene-exposure association (`coeff_G_X`), coefficient of gene-outcome association (`coeff_G_Y`), and the relevant standard errors of these estimates. The values of parameters inputted were also returned to aid in further testing of data/model generation, i.e. actual gene-exposure associations (`gamma`), pleiotropic effects of invalid instruments (`alpha`), additional pleiotropic effects when \acr{InSIDE} assumption not satified (`phi`), causal effect of exposure on outcome (`beta`) and the proportion of invalid genetic instruments with pleiotropic effects on the outcome (`prop_invalid`).




\newpage
## Testing Generation of Data and Models {#appendix-sim-test}

A series of test plots were used to verify that data were simulated as intended under the various conditions specified by input parameters. Test plots were not created for the parameters `n_participants`, `n_instruments` or `n_datasets`, as the functioning of these parameters could be readily inferred from the structure of the  datasets outputted, as above.

### Proportion of Invalid Instruments
\leavevmode\newline The `prop_invalid` parameter specifies the proportion of invalid genetic instruments simulated, i.e. the proportion of genetic instruments affecting the outcome via direct/pleiotropic effects, and thus not solely via the exposure of interest. If simulated correctly, increasing the value of `prop_invalid` should increase the number of instruments with pleiotropic effects, i.e. instruments with `alpha` $\ne$ 0. With random error terms set to 0 and no causal effect present (i.e. `rand_error = FALSE` and `causal_effect = FALSE`), the estimated gene-outcome coefficient estimated using any given instrument will equal the pleiotropic effects of that instrument (i.e. `coeff_G_Y = alpha`), and therefore will only be non-zero for invalid instruments with non-zero pleiotropic effects on the outcome . Plotting `coeff_G_Y` against `alpha` for simulated data with no causal effect or random error should therefore yield a graph where

- For valid instruments: gene-outcome coefficient = alpha = 0
- For invalid instruments:  gene-outcome coefficient = alpha $\ne$  0, with values spread uniformly between `alpha_min` and `alpha_max`




\newpage
Similarly, with random error terms set to 0 (`rand_error = FALSE`) and no causal effect present (`causal_effect = FALSE`), gene-exposure coefficients estimated for each instrument should exactly match the actual values simulated, i.e. `coeff_G_X = gamma` for all instruments:



\newpage
### Gene-Exposure Coefficient Versus Gene-Outcome Coefficient Plots
\leavevmode\newline For the next phase of testing, a function (`plot_GY_GX`) was written to plot the coefficients for gene-exposure versus gene-outcome as estimated using the previously created linear models:



\newpage
With random error terms set to 0 (`rand_error = FALSE`) and no causal effect present, a graph of gene-exposure coefficients versus gene-outcome coefficients should be a straight line through the origin with gradient = 0; causal effect of $\beta$ = 0.1  present (`beta_val = 0.1`, `causal_effect = TRUE`), the slope of a graph of gene-exposure coefficients versus gene-outcome coefficients from the same sample should be a straight line through the origin with gradient = 0.1:




\newpage
### Random Errors
\leavevmode\newline Re-plotting the same graphs with non-zero random error terms (`rand_error = TRUE`) should produce similar graphs with Gaussian spread around lines passing through the origin with gradients of 0 and 0.1 for no causal effect and causal effect, respectively:



\newpage
### One versus Two Sample MR
\leavevmode\newline Where gene-exposure coefficients and gene-outcome coefficients are estimated from two separate samples rather than one (i.e. `two_sample = TRUE`, simulating 2 sample MR), even with random error terms set to zero, error will be introduced into causal effect estimation through random sampling of different combinations of effect alleles. However, where a causal effect is not present, the effect estimated will consistently be zero regardless of the combinations of alleles sampled, so random error should not be introduced:



\newpage
### Invalid Instruments 
\leavevmode\newline Where invalid instruments are present (i.e. `prop_invalid` $\ne$ `0`) and random error terms are set to 0, graphs of gene-exposure coefficients versus gene-outcome coefficients should be straight lines through the origin and all points representing valid instruments; the invalid instruments should appear as outliers to this line:



\newpage
### Balanced Versus Directional Pleiotropy
\leavevmode\newline Replotting the above with unbalanced pleiotropy present (`balanced_pleio = FALSE`), the invalid instruments should all appear as outliers in the positive direction, i.e. steepening the line of best fit and leading to overestimation of the causal effect: 


<!-- https://pmc.ncbi.nlm.nih.gov/articles/PMC4469799/ -->
<!-- https://mr-dictionary.mrcieu.ac.uk/term/inside/ -->

\newpage
### InSIDE Assumption and Phi
\leavevmode\newline The variable phi represents additional pleiotropic effects of each invalid instrument when the \acr{InSIDE} assumption is not satisfied. The \acr{InSIDE} assumption states that the gene-exposure association is not correlated with the pleiotropic path gene-outcome path of any invalid genetic instruments. This assumption can be violated if e.g.:

- several invalid genetic instruments influence the outcome via the same pleiotropic path

- several invalid genetic instruments are related to the same (unmeasured) confounders of the exposure:outcome relationship, aka correlated pleiotropy. 

As such, when the \acr{InSIDE} assumption is violated, even "strong" instruments (i.e. those with a strong gene-exposure relationship) may not allow accurate estimation of the true causal effect, as pleiotropic effects may scale with instrument strength. If pleiotropic effects are balanced, InSIDE assumption violation may lead to greater imprecision in causal effect estimation; if pleiotropic effects are directional, \acr{InSIDE} assumption violation may lead to bias.

Bowden et al [@bowden_consistent_2016] modeled phi as the pleiotropic effects of unmeasured genetic confounders of the exposure:outcome relationship.  Phi adds additional error to causal effect estimation in scenarios with directional pleiotropic effects (`0 < alpha < 0.2`) and \acr{InSIDE} assumption violation. As such, switching `InSIDE_satisfied` from `TRUE` to `FALSE` should add scatter to the linear association expected when plotting alpha versus gene-outcome coefficients with random error terms set to zero:




\newpage
Setting `InSIDE_satisfied = TRUE` should mean `phi = 0`; `InSIDE_satisfied=FALSE` should result in `phi` $\propto$ gene-outcome coefficient, with scatter only in the positive direction of gene-outcome coefficients given the model also requires directional pleiotropy before `phi` is used:




\newpage
## Summary Table {#appendix-sim-summ}

A function (`get_summary_MR_tib_row`) was written to take models generated from each simulated dataset, estimate causal effect using both weighted median and MR-Hevo methodologies, then output a summary formatted as per Tables 2 & 3 in Bowden et al [@bowden_consistent_2016]:











\newpage

# Appendix: R Packages Used {#appendix-pkg}

## Package Citations

This work was completed using R version 4.4.3 [@base] with the following R packages: acronymsdown v. 0.11.1 [@acronymsdown], bookdown v. 0.43 [@bookdown2016; @bookdown2025], car v. 3.1.3 [@car], cowplot v. 1.1.3 [@cowplot], crayon v. 1.5.3 [@crayon], devtools v. 2.4.5 [@devtools], ggdag v. 0.2.13 [@ggdag], gghighlight v. 0.4.1 [@gghighlight], grateful v. 0.2.12 [@grateful], grid v. 4.4.3 [@grid], here v. 1.0.1 [@here], infer v. 1.0.8 [@infer], kableExtra v. 1.4.0 [@kableExtra], knitr v. 1.50 [@knitr2014; @knitr2015; @knitr2025], matrixStats v. 1.5.0 [@matrixStats], medicaldata v. 0.2.0 [@medicaldata], parallel v. 4.4.3 [@parallel], rmarkdown v. 2.29 [@rmarkdown2018; @rmarkdown2020; @rmarkdown2024], rstan v. 2.32.7 [@rstan], tidyverse v. 2.0.0 [@tidyverse], TwoSampleMR v. 0.6.16 [@TwoSampleMR2017; @TwoSampleMR2018], wordcountaddin v. 0.3.0.9000 [@wordcountaddin].

## Session Information


```
##  setting  value
##  version  R version 4.4.3 (2025-02-28 ucrt)
##  os       Windows 11 x64 (build 26100)
##  system   x86_64, mingw32
##  ui       RTerm
##  language (EN)
##  collate  English_United Kingdom.utf8
##  ctype    English_United Kingdom.utf8
##  tz       Europe/London
##  date     2025-06-07
##  pandoc   3.4 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
##  quarto   NA @ C:\\PROGRA~1\\RStudio\\resources\\app\\bin\\quarto\\bin\\quarto.exe
```

```
## # A tibble: 22 x 5
##    package      ondiskversion loadedversion date       source                   
##    <chr>        <chr>         <chr>         <chr>      <chr>                    
##  1 acronymsdown 0.11.1        0.11.1        2025-06-07 Github (rchaput/acronyms~
##  2 bookdown     0.43          0.43          2025-04-15 CRAN (R 4.4.3)           
##  3 cowplot      1.1.3         1.1.3         2024-01-22 CRAN (R 4.4.3)           
##  4 dplyr        1.1.4         1.1.4         2023-11-17 CRAN (R 4.4.3)           
##  5 forcats      1.0.0         1.0.0         2023-01-29 CRAN (R 4.4.3)           
##  6 gghighlight  0.4.1         0.4.1         2023-12-16 CRAN (R 4.4.3)           
##  7 ggplot2      3.5.2         3.5.2         2025-04-09 CRAN (R 4.4.3)           
##  8 grateful     0.2.12        0.2.12        2025-04-30 CRAN (R 4.4.3)           
##  9 here         1.0.1         1.0.1         2020-12-13 CRAN (R 4.4.3)           
## 10 infer        1.0.8         1.0.8         2025-04-14 CRAN (R 4.4.3)           
## 11 kableExtra   1.4.0         1.4.0         2024-01-24 CRAN (R 4.4.3)           
## 12 lubridate    1.9.4         1.9.4         2024-12-08 CRAN (R 4.4.3)           
## 13 medicaldata  0.2.0         0.2.0         2021-08-16 CRAN (R 4.4.3)           
## 14 purrr        1.0.4         1.0.4         2025-02-05 CRAN (R 4.4.3)           
## 15 readr        2.1.5         2.1.5         2024-01-10 CRAN (R 4.4.3)           
## 16 rstan        2.32.7        2.32.7        2025-03-10 CRAN (R 4.4.3)           
## 17 StanHeaders  2.32.10       2.32.10       2024-07-15 CRAN (R 4.4.3)           
## 18 stringr      1.5.1         1.5.1         2023-11-14 CRAN (R 4.4.3)           
## 19 tibble       3.2.1         3.2.1         2023-03-20 CRAN (R 4.4.3)           
## 20 tidyr        1.3.1         1.3.1         2024-01-24 CRAN (R 4.4.3)           
## 21 tidyverse    2.0.0         2.0.0         2023-02-22 CRAN (R 4.4.3)           
## 22 TwoSampleMR  0.6.16        0.6.16        2025-06-05 https://mrcieu.r-univers~
```


\newpage

<!--chapter:end:9_Appendices.Rmd-->

