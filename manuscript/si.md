---
title: Online Supplementary Information
author:
- name: Parris T Humphrey
  affiliation: Department of Organismic and Evolutionary Biology, Harvard University, Cambridge, MA. USA \newline Department of Ecology and Evolutionary Biology, University of Arizona, Tucson, AZ, USA
  email: phumphrey@g.harvard.edu
- name: Trang T Satterlee
  affiliation: Department of Ecology and Evolutionary Biology, University of Arizona, Tucson, AZ, USA
  email: ttnguyen@arizona.edu
- name: Noah K Whiteman
  affiliation: Department of Integrative Biology, University of California, Berkeley, CA. USA \newline Department of Ecology and Evolutionary Biology, University of Arizona, Tucson, AZ, USA
  email: whiteman@berkeley.edu
thanks: "Code and data available at \\texttt{https://github.com/phumph/competitive\\_hierarchies}."
keywords: Pseudomonas, indirect interactions, phyllosphere, microbiome, phytopathogen
csl: manuscript/ecology-letters.csl
bibliography: manuscript/competition.bib

header-left: "\\textbf{\\textit{Competition in the phyllosphere}}"
header-center: "\\hspace{1cm}"
header-right: "Humphrey *et al.* (2020) \\textbf{OSM}"
footer-left: "v.2020-09-27"
footer-right: "\\thepage"
---

# Supplemental Methods

## *In vitro* growth conditions

We characterized canonical growth traits of each bacterial strain separately by tracking cell density through a complete lag--exponential--saturation cycle *in vitro*. We re-streaked each strain from -80 C stocks onto King's B (KB) agar plates +10 mM MgSO<sup>4</sup> and incubated at 28 C for 3 days. We picked and inoculated single colonies into 1 mL minimal media at pH 5.6 ('MM'; 10 mM fructose, 10 mM mannitol, 50 mM KPO<sub>4</sub>, 7.6 mM (NH<sub>4</sub>)<sub>2</sub>SO<sub>4</sub>, and 1.7 mM MgCl<sub>2</sub>; [@Mudgett99a; @Barrett11a]) and grew cultures overnight in a shaking incubator (250 rpm) at 28 C. MM at pH 5.6 induces expression of the type-III secretion system (T3SS) in a diversity of *Pseudomonas* spp. [@Huynh89a], in contrast to KB, which results in negligible T3SS expression. T3SS expression was important for maximizing the potential relevance of our *in vitro* assay environments to those of plants, in which T3SS expression is expected. Thus, MM was used for this and all other culture experiments in this study.

Overnight cultures of each strain in MM were spun down for 3 m at $3000 \times \text{g}$ and the supernatant replaced with 500 $\mu$L fresh MM. The density of each culture was adjusted to $\text{OD}_{600} = 0.2$ prior to 1:100 dilution into a total of 180 $\mu\text{L}$ MM inside the wells of sterile polystyrene 96-well plates (Falcon). Each 96-well plate was covered with optically clear, gas-permeable plastic tape (Sigma `Z380059`) and incubated for 60 hr in a BioTek 600 plate reader in which $\text{OD}_{600}$ measurements were taken every 5 min with continuous random orbital shaking. Identical growth assays were performed on separate days in duplicate.

## Estimating traits related to growth

We used **R** package `grofit` [@Kahm10a] to fit smoothed functions to the bacterial growth data. Curve fits generated using logistic, Richards, Gompertz, or modified Gompertz equations failed to produce estimates with $r \geq 0.5$ and we therefore used a non-parametric locally-weighted smoothing function to estimate the following growth curve parameters: maximum growth rate $r_m$, lag phase $l$, and maximum yield $K$. Lag phase represents the length of time (min) prior to initiation of exponential growth, while maximum yield is the maximum $\text{OD}_{600}$ attained during 60 h of growth. Curves for long lag-phased strains never leveled off (Fig. S1, e.g., strain `17A`); in these cases, $K$ was set as the final $\text{OD}_{600}$. When growth trajectories exhibited multiple exponential phases (diauxic shifts), $r_m$ was estimated during the initial exponential phase (e.g., strain `20A`; Fig. S1).

## Pairwise competition assays

In 100 mm diameter Petri dishes, we overlaid 4 mL of 0.5% (w/v) soft agar inoculated with a resident strain overtop a sterile base layer of growth medium (minimal medium + 1% w/v agar). The top agar overlay for each plate contained a suspension of a single resident strain inoculated at $5 \times 10^{5} \text{ml}^{-1}$ while soft agar was still molten ~42 C. Once cooled, suspensions of each of the 40 invader strains were then spotted at the same concentration in 4 $\mu$L aliquots spaced every 0.5 cm in parallel rows using an 8-channel pipettor. Resident and invader strain suspensions were made from exponential phase cultures in MM 3 mL with shaking at 28 C. Plates were incubated at 28 C face up for 12 h and then face down incubation for an addition 10 days. Megacolony spots were scored by hand for growth after days 1, 3, 5, 8, and 10. Data used for the following analyses are from day 10, by which time all interactions dynamics had leveled off.

\newpage
# Supplemental Figures

<!-- ********* -->
<!-- FIGURE S1 -->
<!-- ********* -->

![**Growth cycle profiles of 40 *Pseudomonas* spp. strains measured individually in minimal media *in vitro*.** Overlain individual growth curves are shown grouped by nearest phylogenetic neighboring strains (see Fig. 1, main text). Strain `08B` is sister to *P. syringae* clade. `4326` = laboratory strain *P. syringae* pv. maculicola str. `ES4326`. \label{figS1}](manuscript/figures/figureS1.png){ width=450px }

<!-- ********* -->
<!-- FIGURE S2 -->
<!-- ********* -->

\newpage
![**Distribution of interaction outcomes within and between *Pseudomonas* clades**. Accompanying outcome counts and statistical results are displayed in Tables S1 and S2, respectively. Total strain pairings = 772 between 40 strains, with the 20 self-interactions and 8 no-data interactions removed. `RI` = reciprocal invasibility; `RNI` = reciprocal non-invasibility; `Asymm.` = asymmetric dominance. \label{figS2}](manuscript/figures/figureS2.png){ width=250px }

<!-- ********* -->
<!-- FIGURE S3 -->
<!-- ********* -->

\newpage
![**Morphological variation observed during pairwise inhibition assays reflects strain-specific induction of motility and inhibition phenotypes in *P. fluorescens*.** **a-c**. Photos of interaction plates against *P. syringae* resident strains against which *P. fluorescens* strain `43A` adopted distinct megacolony morphologies (white arrows). Black arrows denote invading strains producing tox+ halos against the resident strains. **d-f**. close-ups of smooth morph (SM), wrinkly spreader-like (WS-like), and the smooth spreader (SS) morphs from above. **g**. Close-up of putative growth facilitation of non-focal strain via secretion of effectors from 43A, arising from either cell cycle modification by secreted regulators, or competitive release via killing of resident strain. H. Phylogenetic specificity of morphological induction and its relationship with toxicity for strain `43A`. Phylogenetic tree from Fig. 1 (main text). Scale bars in **a–c** = 1 cm; **d–f** = 0.50 cm; **g** = 0.25 cm. \label{figS3}](manuscript/figures/figureS3.png){ width=450px }

<!-- ********* -->
<!-- FIGURE S4 -->
<!-- ********* -->

\newpage
![**Life history trait correlations exhibited by *P. syringae* and *P. fluorescens* strains.** Slopes of clade-wise linear regression are over-plotted in each panel for illustration purposes only regardless of $p$ value of correlation coefficient (see Fig. 1c, main text). \label{figS4}](manuscript/figures/figureS4.png){ width=450px }

<!-- ********* -->
<!-- FIGURE S5 -->
<!-- ********* -->

\newpage
![**Correlations among competitive traits exhibited by *P. syringae* and *P. fluorescens* strains.** Note that $C_{t}$ for *P. syringae* are all 0 because none produced detectable inhibitors scored in our pairwise soft-agar competition assays (see Extended Data Fig. 1, main text). Slopes of clade-wise linear regression are over-plotted in each panel for illustration purposes only regardless of $p$ value of correlation coefficient (see Fig. 1c, main text). \label{figS5}](manuscript/figures/figureS5.png){ width=400px }

<!-- ********* -->
<!-- FIGURE S6 -->
<!-- ********* -->

\newpage
![**Correlations between competitive traits and life-history traits exhibited by *P. syringae* and *P. fluorescens* strains.** Slopes of clade-wise linear regression are over-plotted in each panel for illustration purposes only regardless of $p$ value of correlation coefficient (see Fig. 1c, main text). \label{figS6}](manuscript/figures/figureS6.png){ width=325px }

<!-- ********* -->
<!-- FIGURE S7 -->
<!-- ********* -->

\newpage
![**Multivariate trait dispersion within *P. syringae* and *P. fluorescens* clades**. Pairwise Euclidean distances along trait PCs 1-3 are greater for *P. syringae* strains than among *P. fluorescence* strains (unequal variance two-sample Welch's $t$-test, $t = 8.71; p < 10^{-5}$). \label{figS7}](manuscript/figures/figureS7.png){ width=150px }

<!-- ********* -->
<!-- FIGURE S8 -->
<!-- ********* -->

\newpage
![**Strain-wise and overall rank differences with and without indirect interactions from secretion-producing *P. fluorescens* strains**. **a-b.** Indirect effects on *P. syringae* reduce fitness most for strains with high base rank and facilitate weaker strains. Community-wide effects of indirect facilitation is to weaken fitness hierarchies among *P. syringae* strains (base rank versus final rank rank correlation = 0.52). **c-d.** In *P. fluorescens*, fitness ranks are far less affected by intra-clade indirect effects than in *P. syringae* (base rank versus final rank rank correlation $= 0.94$). Strains along the $y$ axis are plotted from high (top) to low (bottom) baseline fitness ($C_{w}$), while *P. fluorescens* strains are plotted from most facilitating (left) to least facilitating (right) along the $x$ axis. \label{figS8}](manuscript/figures/figureS8.png){ width=450px }

\newpage

# Supplemental Tables

\input{manuscript/tables}

\clearpage
\newpage

# References
