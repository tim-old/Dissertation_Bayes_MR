# References

<div id="refs"></div>

<!-- export Zotero as Western encoding -->

# (APPENDIX) Appendix {-} 

\newpage
 
# Appendix A: List of Abbreviations
 
\printacronyms
 
\newpage
 
# Appendix B: Simulation Code


 
### Generating Data and Models

 The data generating model used was from Appendix 3 of Bowden et al (ref); the relevant section describing their model is reproduced below:

 _"..._

 \begin{equation}
 U_i = \sum^J_{j=1} \phi_jG_{ij} + \epsilon_i^U
 \end{equation}


 \begin{equation}
 X_i = \sum^J_{j=1} \gamma_jG_{ij} + U_i + \epsilon_i^X
 \end{equation}

 \begin{equation}
 Y_i = \sum^J_{j=1} \alpha_jG_{ij} + \beta X_i + U_i + \epsilon_i^Y
 \end{equation}

 _for participants indexed by $i = 1, . . . , N$, and genetic instruments indexed by $j = 1, . . . , J$._

 _The error terms $\epsilon_i^U , \epsilon_i^X$ and $\epsilon_i^Y$ were each drawn independently from standard normal distributions. The genetic effects on the xposure γj are drawn from a uniform distribution between 0.03 and 0.1. Pleiotropic effects $\alpha_j$ and $\phi_j$ were set to zero if the genetic instrument was a alid instrumental variable. Otherwise (with probability 0.1, 0.2, or 0.3):_

 _1. In Scenario 1 (balanced pleiotropy, InSIDE satisfied), the $\alpha_j$ parameter was drawn from a uniform distribution between −0.2 and 0.2._

 _2. In Scenario 2 (directional pleiotropy, InSIDE satisfied), the $\alpha_j$ parameter was drawn from a uniform distribution between 0 and 0.2._

 _3. In Scenario 3 (directional pleiotropy, InSIDE not satisfied), the $\phi_j$ parameter was drawn from a uniform distribution between −0.2 and 0.2._


 _The causal effect of the exposure on the outcome was either $\beta X = 0$ (null causal effect) or $\beta X = 0.1$ (positive causal effect). A total of 10 000 simulated datasets were generated for sample sizes of N = 10 000 and 20 [sic] participants. Only the summary data, that is genetic associations with the exposure and with the outcome and their standard errors as estimated by univariate regression on the genetic instruments in turn, were used by the analysis methods. In the two-sample setting, data were generated on 2N participants, and genetic associations with the exposure were estimated in the first N participants, and genetic associations with the outcome in the second N participants."_ @bowden_consistent_2016

 To reproduce this model, code was written in R to generate the relevant participant level data. First, a function (`simulate_MR_data`) was written which included arameters specified by Bowden et al, and also to allow testing of data simulation:





 This initial simulation function generated data in the following format:




 A function (`extract_models`) was then written to create linear models from each dataset generated as per Bowden et al:





 These models generated estimates of the coefficient of gene:exposure association (`coeff_G_X`), coefficient of gene:outcome association (`coeff_G_Y`), and the elevant standard errors of these estimates. The values of parameters inputted were also returned to aid in further testing of data/model generation, i.e. actual ene:exposure associations (`gamma`), pleiotropic effects of invalid instruments (`alpha`), additional pleiotropic effects when InSIDE assumption not satified `phi`), causal effect of exposure on outcome (`beta`) and the proportion of invalid genetic instruments with pleiotropic effects on the outcome (`prop_invalid`).



### Testing Generation of Data and Models

 A series of test plots were used to verify that data were simulated as intended under the various conditions specified by input parameters. Test plots were not reated for the parameters `n_participants`, `n_instruments` or `n_datasets`, as the functioning of these parameters could be readily inferred from the structure of he  datasets outputted, as above.

 The `prop_invalid` parameter specifies the proportion of invalid genetic instruments simulated, i.e. the proportion of genetic instruments affecting the outcome via irect/pleiotropic effects, and thus not solely via the exposure of interest. If simulated correctly, increasing the value of `prop_invalid` should increase the umber of instruments with pleiotropic effects, i.e. instruments with `alpha` =/= 0. With random error terms set to 0 and no causal effect present (i.e. `rand_error = ALSE` and `causal_effect = FALSE`), the estimated gene:outcome coefficient estimated using any given instrument will equal the pleiotropic effects of that instrument i.e. `coeff_G_Y = alpha`), and therefore will only be non-zero for invalid instruments with non-zero pleiotropic effects on the outcome . Plotting `coeff_G_Y` gainst `alpha` for simulated data with no causal effect or random error should therefore yield a graph where

   - For valid instruments: gene:outcome coefficient = alpha = 0
   - For invalid instruments:  gene:outcome coefficient = alpha =/=  0, with values spread uniformly between `alpha_min` and `alpha_max`




 Similarly, with random error terms set to 0 and no causal effect present, gene:exposure coefficients estimated for each instrument should exactly match the actual alues simulated, i.e. `coeff_G_X = gamma` for all instruments:



 For the next phase of testing, a function (`plot_GY_GX`) was written to plot the coefficients for gene:exposure versus gene:outcome as estimated using the reviously created linear models:








































\newpage

# Appendix C: Citation Search Strategy



