---
title: manuscript
author: parris t humphrey
---

## Indirect effects of antibiotic production can neutralize competitive hierarchies among phyllosphere bacteria

#### Parris T. Humphrey<sup>1,2</sup>, Trang N. Nguyen<sup>2</sup>, Noah K. Whiteman<sup>1,2</sup>

<sup>1</sup>Organismic & Evolutionary Biology, Harvard University, Cambridge, MA 02138, USA<br>
<sup>2</sup>Integrative Biology, University of California, Berkeley, CA 94720, USA<br>
<sup>3</sup>Department of Plant Sciences, University of Arizona, Tucson, AZ 85721, USA

Corresponding author: PTH ([phumphrey@g.harvard.edu](phumphrey@g.harvard.edu))

<hr>

### Abstract

Competition is rife within bacterial communities, but we poorly understand the traits that underlie various types of competitive fitness, the correlations among them, and their distributions in natural bacterial populations. In this study, we characterize the types and strengths of competitive interactions that occur among a suite of *Pseudomonas* spp. strains isolated from a native forb at a single field site. Outcomes from pairwise competitions in spatial microcosms revealed strong competitive dominance of *Pseudomonas fluorescens* strains over *P. syringae* strains. *P. fluorescens* strains were better resources competitors and more often produced antibiotics to which few *P. syringae* strains expressed resistance. Within each bacterial clade, growth rates of individual strains were not correlated with overall competitive ability or saturation density *in vitro*; however, *P. fluorescens* strains with shorter lag times *in vitro* were on average better competitors, suggesting that spatial preemption of a shared resource promoted competitive fitness in the context of our experiments. We then quantified how patterns of resistance among *P. syringae* to *P. fluorescens* secretions might influence relative competitive ranks *within* the *P. syrinage* clade. This analysis showed that a high proportion of competitive outcomes would be reversed in the presence of antibiotics produced by *P. fluorescens*, and these indirect benefits accrued most to strains of lower average competitiveness. Thus, *P. fluorescens* are not only potent direct competitors of putatively phytopathogenic *P. syringae* strains but can also indirectly change the outcome of a large number of competitive interactions among them when present within the same spatial context. Testing the relevance of these interaction types in the context of living plant tissues will help reveal the role of *P. fluorescens* in suppressing the growth and/or upending the fitness ranks of potentially phytopathogenic *P. syringae*.


**Keywords:** interference competition; species interactions; *Pseudomonas*; co-existence; facilitation


**Author contributions:** PTH and NKW designed the study; PTH and TNN collected the data; PTH wrote the manuscript with help from NKW.

<hr>

### Introduction



### Materials & Methods

#### Bacterial strains

Of the 51 *Pseudomonas* spp. strains isolated from bittercress and described by Humphrey et al. (2014), we selected a set of 39 (26 *P. syringae*, 14 *P. fluorescens*) that represent the extent of observed diversity. The laboratory strain *P. syringae* pv. maculicola str. ES4326 (hereafter Psm4326) was used as a reference owing to its phylogenetic similarity to strains isolated from bittercress and its extensive characterization in the laboratory as a pathogen of *Arabidopsis thaliana* (Cui:2005dn, Cui:2002gp, Groen:2013bt, Groen:2015bv).

#### Competition assays

Pairwise competition assays were conducted on 1% agar MM plates (100 mm diameter) onto which we overlaid 4 ml of 0.5% soft agar MM containing a bacterial suspension of each resident strain inoculated at $5 \times 10^{5} \text{ml}^{-1}$ while soft agar was still molten ($\sim 42 \text{ C}$). Suspensions of each of the 40 invader strains were then spotted at the same concentration in 4 $\mu\text{L}$ aliquots spaced every 0.5 cm in parallel rows using an 8-channel pipettor. Plates were incubated face up for 12 h, followed by face down incubation at $28\text{ C}$ for 10 days. Megacolony spots were scored by hand for growth on days 1, 3, 5, 8, and 10. Data used for the following analyses are from day 10, by which time all interactions dynamics had reached equilibrium.

We scored growth of each invader as $0$ for no visible growth of the invader above a negative control spot containing MM alone, $0.5$ for a largely translucent 'megacolony', which reflected a definite presence of
growth but which was relatively suppressed and confined to the megacolony margin, and $1$ for obvious and robust megacolony growth. Examples of each can be seen in Fig. S3. We scored inhibition interactions as the presence of a zone of clearance (halo) $\geq 1 \text{mm}$ surrounding the extent of the invader megacolony (Fig. S3). Inhibition interactions were ultimately scored as 0 or 1 regardless of the spatial extent of the halo, although variation in halo width was recorded. We also scored any morphological variation among megacolonies for particular strains, and later we relate such variation for particular strains to the phylogenetic position of their competitors (see *Analyzing the distribution of competitive outcomes*, below).

#### Calculating indexes of competitiveness

Each strain was assayed under 40 different conditions both as resident strain and as invader, comprising an interaction network with 1600 entries (including self vs. self). One version of the interaction network represents the outcome of resource competition and details the extent of growth of each invader, while the other captures the presence or absence of inhibitory interactions indicated by zones of clearance in the resident population. For resource competitions, we calculate the invasiveness ($C_o$) and defense capacity (i.e. territoriality; $C_d$) of each strain. $C_o$ for each strain $i$ was calculated as $C_{o,i} = \frac{1}{n_{\text{ij}}}\sum_{i \neq j}^{n}x_{ij}$, where $x_{ij} \in \{ 0,0.5,1\}$ and $n_{ij}$ is the total number of scored interactions for each strain as the invader with all non-self resident strains. $C_o$ is thus the expected value of growth attained by each strain as the invader across the population of residents.

$C_d$ was calculated similarly except the focal strain $j$ is in the resident state, $x_{ji} \in \{ 0,0.5,1\}$ is as before but has a subscript reversal and indicates the degree to which the resident prevented the growth of each invader $i$, and $n_{ji}$ is the number of interactions occurring between each focal resident and its non-self invaders. $C_d$ can thus be interpreted as the expected amount of growth each resident strain can prevent among the population of invaders assayed. We calculated an overall exploitative competition index $C_w$ for each strain as $C_{o,i} - {(1 - C}_{d,i})$, where $-1 \leq C_{w} \leq 1$. These extremes represent absolute competitive inferiority ($-1$), where a strain failed to prevent any growth of any invader and similarly failed to invade any other strain, to absolute competitive dominance ($1$), where a strain fully invaded all residents and fully prevented growth of all invaders.

We also calculated $C_t$ and $C_r$ based on the interaction matrix for interference competition, where $C_{t}$ is the proportion of successful invasions (i.e. given growth of 0.5 or above) that also resulted in halo formation, indicative of inhibition of the resident. $C_r$ for a strain is the proportion of contests with all invading inhibitor strains (i.e., all strains with $C_t > 0$) that failed to result in halo formation, which we took as evidence of resistance. An overall interference competition index $I_w$ was calculated for each strain as $I_{w,i} = C_{t,i} - {(1 - C}_{r,i})$, where $-1 \leq I_w \leq 1$, which is equal to the aggressiveness index (AI) of Vestigian et al. (2011).

#### *In vitro* growth assays

Strains were re-streaked from $-80 \text{ C}$ stocks in 50% glycerol onto King's B (KB) agar plates +10 mM MgSO$^4$ and incubated at $28 C$ for 3 days. Bittercress isolates had undergone only one prior cycle of isolation--growth--freezing since initial isolate on KB plates from surface-sterilized homogenates from bittercress leaf samples (Humphrey et al 2014). Single colonies were picked and inoculated into 1 mL minimal media at pH 5.6 ('MM'; 10 mM fructose, 10 mM mannitol, 50 mM KPO$^{4}$, 7.6 mM (NH$_{4}$)$^{2}$SO$^4$, and 1.7 mM MgCl$^{2}$; Mudgett:1999tu, Barrett:2011fo) and grown overnight in a shaking incubator (250 rpm) at $28 \text{ C}$. MM at pH 5.6 has been shown to induce the expression of the type-III secretion system (T3SS) in a diversity of *Pseudomonas* spp. (Huynh:1989ux), in contrast to KB, which results in negligible T3SS expression. T3SS expression was important for maximizing the potential relevance of our in vitro assay environments to those of plants, in which T3SS expression is expected. Each 1 mL overnight culture was spun down for 3 m at $3000 \times \text{g}$ and the supernatant was replaced with 500 $\mu\text{L}$ fresh MM. The density of each culture was adjusted to $\text{OD}_{600} = 0.2$ prior to 1:100 dilution into a total of 180 $\mu\text{L}$ MM inside the wells of sterile polystyrene 96-well plates (Falcon). Each 96-well plate was covered with optically clear, gas-permeable plastic tape (Sigma \#Z380059) and incubated for 60 hr in a BioTek 600 plate reader in which $\text{OD}_{600}$ measurements were taken every 5 min with continuous orbital shaking. Identical growth assays were performed on separate days in duplicate.

#### Estimating resource competition traits

We used R package `grofit` (Kahm:2010vv) to fit smoothed functions to the bacterial growth data. Curve fits generated using logistic, Richards, Gompertz, or modified Gompertz equations failed to produce estimates with $r \geq 0.5$ and we therefore used a non-parametric locally-weighted smoothing function to estimate the following growth curve parameters: maximum growth rate $r_m$, lag phase $L$, and maximum yield $K$. Lag phase represents the length of time (min) prior to initiation of exponential growth, while maximum yield is the maximum $\text{OD}_{600}$ attained during 60 h of growth. Curves for long lag-phased strains never leveled off (Fig. S1, e.g. strain 17A); in these cases, $K$ was set as the final $\text{OD}_{600}$. When growth trajectories exhibited multiple exponential phases (diauxic shifts), $r_m$ was estimated during the initial exponential phase (e.g. strain 20A; Fig. S1).

To examine fundamental axes of trait co-variances, we conducted principle components analysis (PCA) using the matrix of mean-centered and scaled competitive indexes and growth parameters for all strains (40 x 9 matrix) as input. We also constructed linear multiple regression models to estimate the contribution of $r_m$, $L$, and $K$ to variation among *P. syringae* and *P. fluorescens* strains in each of the overall competitive indexes $C_w$ and $I_w$.

#### Analyzing the distribution of competitive outcomes

We determined when the outcomes of all pairwise interactions between strains $i$ and $j$ ($i \neq j$) took following forms: reciprocal invasibility (RI), where strains $i$ and $j$ each invade one another; reciprocal non-invasibility (RNI), where strains $i$ and $j$ cannot invade each other; and asymmetric (AS), where strain $i$ invades strain $j$ but $j$ cannot invade $i$. We constructed binomial generalized linear models (GLMs) in \texttt{R} with the canonical logit link function to estimate the probability of RI, RNI, and AS as a function of genetic distance as well as trait distance between strains. Genetic distance ($D_g$) was calculated as the pairwise raw nucleotide distance between 2690 bp of sequence comprised of four partial housekeeping gene sequences previously generated for each strain from Humphrey et al. (2014). Orthologous sequences from the genome of Psm4326 were derived from its published genome sequence (Baltrus et al 2011; RefSeq ID `NZ_AEAK00000000.1`).

Euclidean distances between each growth trait $r_m$, $L$, and $K$ for all pairs of strains were measured as $D_{ij} = \sqrt{{(x_{i} - x_{j})}^{2}}$. We first examined a binomial model for each outcome type using $D_g$ as the only fixed effect, and then computed models including each growth trait, which took the form

$$
\text{logit}\left(P(y_{\text{ij}}|x_{\text{ij}}) \right)\ \sim\ \beta_{0} + \beta_{d}x_{d} + \beta_{r_{m}}x_{r_{m}} + \beta_{L}x_{L} + \beta_{K}x_{K}
$$

To test for genetic correlations in trait values, we ran Mantel tests between pairs of trait and genetic distance matrixes in R using package `vegan`. Test statistics were compared with those generated from 1000 matrix permutations (veganCommunityEco:2012uw). We noted instances where megacolony morphology differed between strain pairings for particular isolates (e.g., *P. fluorescens* str. RM43A), and we compared the incidence of each discrete phenotype to the phylogenetic position of the competitor strains using the genetic data described above (data from Humphrey:2014ga). To test the significance of a phylogenetic correlation between the phylogenetic position of competitor strains and the induced megacolony morphology of the focal strains, we conducted a permutation analysis of variance (perMANOVA) using `vegan` package with 10000 permutations.

#### Inferring indirect interactions from the pair-wise network

We assembled all possible combinations of strain trios and evaluated whether their patterns of interactions fulfilled the criteria for facilitation described in the Introduction. Specifically, we calculated the net first order linkage effects on each strain serving as the focal strain in the presence of each other as the 'competitor' strain, where the interaction between the two is mediated by a nearby third strain. Briefly, facilitation can occur by strain A releasing strain C from inhibition from B (where A also has to be resistant to B's inhibitors), or from resource competition from the superior competitor B. This analysis remains agnostic to mechanism, but merely calculates the proportion of conditions under which facilitation of an otherwise less-fit competitor is expected to arise. We also determined the prevalence of strain trios resulting in R--P--S intransitivity as well as synergistic inhibition, where strain C is both out-competed by B and inhibited by A, to which the B strain is resistant. For each strain, we calculated the net effect of antagonistic vs. facilitative indirect interactions across all possible trios and compare this to underlying fitness metrics derived from the pair-wise interaction network.

We calculated the magnitude of such opportunity costs of being non-resistant to quantify the predicted dependence of the indirect cost of susceptibility and the indirect benefit of resistance. We compared these two variables to their underlying dependence on relative competitiveness in terms of resource use ($C_w$). Strong competitors lose more by being sensitive, because they would have already won most resource contests. In contrast, weaker resource competitors have much to gain by being resistant, but little to lose: the relative improvement in fitness increases dramatically as more contests are won owing to their increased resistance. We modeled how the costs and benefits of susceptibility and resistance depend on underlying resource competitiveness by simulating so-called first-order linkage indirect effects, where the focal interaction is impacted by a third associate strain.


### Results


### Discussion
