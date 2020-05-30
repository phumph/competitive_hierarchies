---
title: "Online Supplemental Material"
author: "Parris T Humphrey"
---

# Online Supplemental Material:


## Competitive hierarchies, antibiosis, and the distribution of bacterial life history traits in a microbiome

#### Parris T. Humphrey<sup>1,2,3*</sup>, Trang N. Nguyen<sup>4</sup>, Noah K. Whiteman<sup>2</sup>

<sup>1</sup>Department of Organismic & Evolutionary Biology, Harvard University, Cambridge, MA, USA<br>
<sup>2</sup>Department of Integrative Biology, University of California, Berkeley, CA, USA<br>
<sup>3</sup>Department of Ecology and Evolutionary Biology, University of Arizona, Tucson, USA<br>
<sup>4</sup>Department of Plant Sciences, University of Arizona, Tucson, AZ, USA

**\*Correspondence and requests for materials** should be addressed to N.K.W. ([whiteman@berkeley.edu](whiteman@berkeley.edu))

<hr>

## Supplemental Methods

### *In vitro* growth conditions

We characterized canonical growth traits of each bacterial strain separately by tracking cell density through a complete lag--exponential--saturation cycle *in vitro*. We re-streaked each strain from -80°C stocks onto King's B (KB) agar plates +10 mM MgSO<sup>4</sup> and incubated at 28°C for 3 days. We picked and inoculated single colonies into 1 mL minimal media at pH 5.6 ('MM'; 10 mM fructose, 10 mM mannitol, 50 mM KPO<sub>4</sub>, 7.6 mM (NH<sub>4</sub>)<sub>2</sub>SO<sub>4</sub>, and 1.7 mM MgCl<sub>2</sub>; Mudgett:1999tu, Barrett:2011fo) and grew cultures overnight in a shaking incubator (250 rpm) at 28°C. MM at pH 5.6 induces expression of the type-III secretion system (T3SS) in a diversity of *Pseudomonas* spp. (Huynh:1989ux), in contrast to KB, which results in negligible T3SS expression. T3SS expression was important for maximizing the potential relevance of our *in vitro* assay environments to those of plants, in which T3SS expression is expected. Thus, MM was used for this and all other culture experiments in this study.

Overnight cultures of each strain in MM were spun down for 3 m at $3000 \times \text{g}$ and the supernatant replaced with 500 µL fresh MM. The density of each culture was adjusted to $\text{OD}_{600} = 0.2$ prior to 1:100 dilution into a total of 180 $\mu\text{L}$ MM inside the wells of sterile polystyrene 96-well plates (Falcon). Each 96-well plate was covered with optically clear, gas-permeable plastic tape (Sigma \#Z380059) and incubated for 60 hr in a BioTek 600 plate reader in which $\text{OD}_{600}$ measurements were taken every 5 min with continuous random orbital shaking. Identical growth assays were performed on separate days in duplicate.

### Estimating traits related to growth

We used R package `grofit` (Kahm:2010vv) to fit smoothed functions to the bacterial growth data. Curve fits generated using logistic, Richards, Gompertz, or modified Gompertz equations failed to produce estimates with $r \geq 0.5$ and we therefore used a non-parametric locally-weighted smoothing function to estimate the following growth curve parameters: maximum growth rate $r_m$, lag phase $l$, and maximum yield $K$. Lag phase represents the length of time (min) prior to initiation of exponential growth, while maximum yield is the maximum OD<sub>600</sub> attained during 60 h of growth. Curves for long lag-phased strains never leveled off (Fig. S1, e.g. strain 17A); in these cases, $K$ was set as the final OD<sub>600</sub>. When growth trajectories exhibited multiple exponential phases (diauxic shifts), $r_m$ was estimated during the initial exponential phase (e.g. strain 20A; Fig. S1).

### Pairwise competition assays

In 100 mm diameter Petri dishes, we overlaid 4 mL of 0.5% (w/v) soft agar inoculated with a resident strain overtop a sterile base layer of growth medium (minimal medium + 1% w/v agar; see ). The top agar overlay for each plate contained a suspension of a single resident strain inoculated at $5 \times 10^{5} \text{ml}^{-1}$ while soft agar was still molten ~42°C. Once cooled, suspensions of each of the 40 invader strains were then spotted at the same concentration in 4 µL aliquots spaced every 0.5 cm in parallel rows using an 8-channel pipettor. Resident and invader strain suspensions were made from exponential phase cultures in MM 3 mL with shaking at 28°C. Plates were incubated at 28°C face up for 12 h and then face down incubation for an addition 10 days. Megacolony spots were scored by hand for growth after days 1, 3, 5, 8, and 10. Data used for the following analyses are from day 10, by which time all interactions dynamics had leveled off.


<br>
<hr>

## Supplemental Figures

<!-- ********* -->
<!-- FIGURE S1 -->
<!-- ********* -->

<br>

<center><img src="https://i.imgur.com/nx1AZMl.png" width = "600"></center><br>

><strong>Fig. S1. Growth cycle profiles of 40 *Pseudomonas* spp. strains measured individually in minimal media <em>in vitro</em></strong>. Overlain individual growth curves are shown grouped by nearest phylogenetic neighboring strains (see Fig. 1, main text). Strain 08B is sister to *P. syringae* clade. 4326 = laboratory strain *P. syringae* pv. maculicola str. ES4326.

<br>
<br>


<!-- ********* -->
<!-- FIGURE S2 -->
<!-- ********* -->

<center><img src="https://i.imgur.com/lHXFOwm.png" width = "300"></center><br>

><strong>Fig S2.</strong> <strong>Distribution of interaction outcomes within and between <i>Pseudomonas</i> clades</strong>. Accompanying outcome counts and statistical results are displayed in Tables S1 and S2, respectively. Total strain pairings = 772 between 40 strains, with the 20 self-interactions and 8 no-data interactions removed. RI = reciprocal invasibility; RNI = reciprocal non-invasibility; Asymm. = asymmetric dominance.

<br>
<br>

<!-- ********* -->
<!-- FIGURE S3 -->
<!-- ********* -->

<center>
<img src="https://i.imgur.com/u3ywT7h.png" width = "450"></center><br>

><strong>Fig. S3.</strong> <strong>Morphological variation observed during pairwise inhibition assays reflects strain-specific induction of motility and inhibition phenotypes in *P. fluorescens*.</strong> **a-c**. Photos of interaction plates against *P. syringae* resident strains against which *P. fluorescens* strain 43A adopted distinct megacolony morphologies (white arrows). Black arrows denote invading strains producing tox+ halos against the resident strains. **d-f**. close-ups of smooth morph (SM), wrinkly spreader-like (WS-like), and the smooth spreader (SS) morphs from above. **g**. Close-up of putative growth facilitation of non-focal strain via secretion of effectors from 43A, arising from either cell cycle modification by secreted regulators, or competitive release via killing of resident strain. H. Phylogenetic specificity of morphological induction and its relationship with toxicity for strain 43A. Phylogenetic tree from Fig. 1 (main text). Scale bars in **a–c** = 1 cm; **d–f** = 0.50 cm; **g** = 0.25 cm.


<br>
<br>

<!-- ********* -->
<!-- FIGURE S4 -->
<!-- ********* -->

<center><img src="https://i.imgur.com/tQIzT4e.png" width = "450"></center><br>

><strong>Fig. S4.</strong> <strong>Life history trait correlations exhibited by *P. syringae* and *P. fluorescens* strains.</strong> Slopes of clade-wise linear regression are over-plotted in each panel for illustration purposes only regardless of $p$ value of correlation coefficient (see Fig. 1c, main text).

<br>
<br>


<!-- ********* -->
<!-- FIGURE S5 -->
<!-- ********* -->

<center><img src="https://i.imgur.com/vwxokxu.png" width = "600"></center><br>

><strong>Fig. S5.</strong> <strong>Correlations among competitive traits exhibited by *P. syringae* and *P. fluorescens* strains.</strong> Note that $C_{t}$ for *P. syringae* are all 0 because none produced detectable inhibitors scored in our pairwise soft-agar competition assays (see Extended Data Fig. 1, main text). Slopes of clade-wise linear regression are over-plotted in each panel for illustration purposes only regardless of $p$ value of correlation coefficient (see Fig. 1c, main text).
<br>
<br>


<!-- ********* -->
<!-- FIGURE S6 -->
<!-- ********* -->

<center><img src="https://i.imgur.com/8OZeMFy.png" width = "400"></center><br>

><strong>Fig. S6.</strong> <strong>Correlations between competitive traits and life-history traits exhibited by *P. syringae* and *P. fluorescens* strains.</strong> Slopes of clade-wise linear regression are over-plotted in each panel for illustration purposes only regardless of $p$ value of correlation coefficient (see Fig. 1c, main text).

<br>
<br>

<!-- ********* -->
<!-- FIGURE S7 -->
<!-- ********* -->

<center><img src="https://i.imgur.com/iu8Mulc.png" width = "200"></center><br>

><strong>Fig. S7.</strong><strong> Multivariate trait dispersion within <i>P. syringae</i> and <i>P. fluorescens</i> clades</strong>. Pairwise Euclidean distances along trait PCs 1-3 are greater for <i>P. syringae</i> strains than among <i>P. fluorescence</i> strains (unequal variance two-sample [Welch's] $t$-test, $t = 8.71; p < 10^{-5}$).

<br>

<!-- ********* -->
<!-- FIGURE S8 -->
<!-- ********* -->
![]()

<center><img src="https://i.imgur.com/Ot3lihb.png" width = "600"></center><br>

><strong>Fig. S8. </strong><strong>Strain-wise and overall rank differences with and without indirect interactions from secretion-producing <i>P. fluorescens</i> strains</strong>. <strong>a-b.</strong> Indirect effects on <i>P. syringae</i> reduce fitness most for strains with high base rank and facilitate weaker strains. Community-wide effects of indirect facilitation is to weaken fitness hierarchies among <i>P. syringae</i> strains (base rank versus final rank rank correlation = 0.52). <strong>c-d.</strong> In <i>P. fluorescens</i>, fitness ranks are far less affected by intra-clade indirect effects than in <i>P. syringae</i> (base rank versus final rank rank correlation = 0.94). Strains along the <i>y</i> axis are plotted from high (top) to low (bottom) baseline fitness (<i>C<sub>w</sub></i>), while <i>P. fluorescens</i> strains are plotted from most facilitating (left) to least facilitating (right) along the <i>x</i> axis.

<br>
<hr>
<br>

## Supplemental Tables

### Table S1. Summary statistics of life history and competitive trait distributions between *Pseudomonas* spp. clades.

<small>
<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> trait </th>
   <th style="text-align:left;"> <i>P. fluorescens</i> </th>
   <th style="text-align:left;"> <i>P. syringae</i> </th>
  </tr>
 </thead>
<tbody>
  <tr grouplength="9"><td colspan="3" style="border-bottom: 1px solid;"><strong>mean (stdev)</strong></td></tr>
<tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>r</i> </td>
   <td style="text-align:left;"> 0.57 (0.15) </td>
   <td style="text-align:left;"> 0.34 (0.22) </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>L</i> </td>
   <td style="text-align:left;"> 672.13 (216.43) </td>
   <td style="text-align:left;"> 1549.86 (656.52) </td>
  </tr>
  <tr>
  <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>K</i> </td>
   <td style="text-align:left;"> 0.49 (0.06) </td>
   <td style="text-align:left;"> 0.42 (0.18) </td>
  </tr>
  <tr>
  <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>i<sub>w</sub></i> </td>
   <td style="text-align:left;"> 0.15 (0.19) </td>
   <td style="text-align:left;"> -0.28 (0.22) </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>c<sub>w</sub></i> </td>
   <td style="text-align:left;"> 0.65 (0.2) </td>
   <td style="text-align:left;"> -0.35 (0.3) </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>c<sub>t</sub></i> </td>
   <td style="text-align:left;"> 0.18 (0.16) </td>
   <td style="text-align:left;"> 0 (0) </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>c<sub>r</sub></i> </td>
   <td style="text-align:left;"> 0.97 (0.05) </td>
   <td style="text-align:left;"> 0.72 (0.22) </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>c<sub>o</sub></i> </td>
   <td style="text-align:left;"> 0.82 (0.1) </td>
   <td style="text-align:left;"> 0.26 (0.13) </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>c<sub>d</sub></i> </td>
   <td style="text-align:left;"> 0.83 (0.11) </td>
   <td style="text-align:left;"> 0.39 (0.2) </td>
  </tr>
  <tr grouplength="9"><td colspan="3" style="border-bottom: 1px solid;"><strong>[min; max]</strong></td></tr>
<tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>r</i> </td>
   <td style="text-align:left;"> [0.41; 0.91] </td>
   <td style="text-align:left;"> [0.03; 1.09] </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>L</i> </td>
   <td style="text-align:left;"> [446.1; 1126.92] </td>
   <td style="text-align:left;"> [470.52; 2716.58] </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>K</i> </td>
   <td style="text-align:left;"> [0.36; 0.57] </td>
   <td style="text-align:left;"> [0.03; 0.71] </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>i<sub>w</sub></i> </td>
   <td style="text-align:left;"> [-0.11; 0.56] </td>
   <td style="text-align:left;"> [-0.77; 0] </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>c<sub>w</sub></i> </td>
   <td style="text-align:left;"> [0.3; 0.88] </td>
   <td style="text-align:left;"> [-0.85; 0.31] </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>c<sub>t</sub></i> </td>
   <td style="text-align:left;"> [0; 0.56] </td>
   <td style="text-align:left;"> [0; 0] </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>c<sub>r</sub></i> </td>
   <td style="text-align:left;"> [0.83; 1] </td>
   <td style="text-align:left;"> [0.23; 1] </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>c<sub>o</sub></i> </td>
   <td style="text-align:left;"> [0.64; 0.94] </td>
   <td style="text-align:left;"> [0.04; 0.67] </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> <i>c<sub>d</sub></i> </td>
   <td style="text-align:left;"> [0.65; 1] </td>
   <td style="text-align:left;"> [0.03; 0.64] </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>a</sup> <i>r</i> displayed as x 1000</td>
</tr>
</tfoot>
</table>
</small>

<br>

### Table S2. Distribution of outcomes by pairing type.

<small>
<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> predictor_col </th>
   <th style="text-align:right;"> ASYM </th>
   <th style="text-align:right;"> RI </th>
   <th style="text-align:right;"> RNI </th>
  </tr>
 </thead>
<tbody>
  <tr grouplength="3"><td colspan="4" style="border-bottom: 1px solid; text-align:left"><strong>Counts</strong></td></tr>
<tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> Between </td>
   <td style="text-align:right;"> 353 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> Within (<i>P. fluorescens</i>) </td>
   <td style="text-align:right;"> 72 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
  <td style="text-align:left; padding-left: 2em;" indentlevel="1"> Within (<i>P. syringae</i>) </td>
   <td style="text-align:right;"> 259 </td>
   <td style="text-align:right;"> 20</td>
   <td style="text-align:right;"> 45 </td>
  </tr>
  <tr grouplength="3"><td colspan="4" style="border-bottom: 1px solid; text-align:left"><strong>Frequencies</strong></td></tr>
<tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> Between </td>
   <td style="text-align:right;"> 0.989 </td>
   <td style="text-align:right;"> 0.006 </td>
   <td style="text-align:right;"> 0.006 </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> Within (<i>P. fluorescens</i>) </td>
   <td style="text-align:right;"> 0.79 </td>
   <td style="text-align:right;"> 0.18 </td>
   <td style="text-align:right;"> 0.03 </td>
  </tr>
  <tr>
   <td style="padding-left: 2em;" indentlevel="1"> Within (<i>P. syringae</i>) </td>
   <td style="text-align:right;"> 0.80 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:right;"> 0.14 </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = "text-align:right; padding: 0; border:0;" colspan='100%'><sup>a</sup> RNI = Reciprocal non-invasion</td>
</tr>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>b</sup> RI = Reciprocal invasion</td>
</tr>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>c</sup> ASYM = Asymmetric dominance</td>
</tr>
</tfoot>
</table>
</small>
<br>


### Table S3. Multinomial model coefficients.

<small>
<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> outcome </th>
   <th style="text-align:left;"> term </th>
   <th style="text-align:right;"> coefficient </th>
   <th style="text-align:right;"> std.err </th>
   <th style="text-align:right;"> <i>z</i> </th>
   <th style="text-align:right;"> <i>p</i> </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> RI </td>
   <td style="text-align:left;"> intercept </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.71 </td>
   <td style="text-align:right;"> -7.30 </td>
   <td style="text-align:right;"> 0.0000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RI </td>
   <td style="text-align:left;"> <i>P. fluorescens</i> </td>
   <td style="text-align:right;"> 39.22 </td>
   <td style="text-align:right;"> 0.76 </td>
   <td style="text-align:right;"> 4.82 </td>
   <td style="text-align:right;"> 0.0000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RI </td>
   <td style="text-align:left;"> <i>P. syringae</i> </td>
   <td style="text-align:right;"> 13.63 </td>
   <td style="text-align:right;"> 0.75 </td>
   <td style="text-align:right;"> 3.50 </td>
   <td style="text-align:right;"> 0.0005 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNI </td>
   <td style="text-align:left;"> intercept </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.71 </td>
   <td style="text-align:right;"> -7.30 </td>
   <td style="text-align:right;"> 0.0000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNI </td>
   <td style="text-align:left;"> <i>P. fluorescens</i> </td>
   <td style="text-align:right;"> 7.36 </td>
   <td style="text-align:right;"> 0.92 </td>
   <td style="text-align:right;"> 2.16 </td>
   <td style="text-align:right;"> 0.0305 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNI </td>
   <td style="text-align:left;"> <i>P. syringae</i> </td>
   <td style="text-align:right;"> 30.67 </td>
   <td style="text-align:right;"> 0.73 </td>
   <td style="text-align:right;"> 4.71 </td>
   <td style="text-align:right;"> 0.0000 </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>a</sup> RNI = Reciprocal non-invasion</td>
</tr>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>b</sup> RI = Reciprocal invasion</td>
</tr>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>c</sup> coefficients = odds</td>
</tr>
</tfoot>
</table>
</small>

<br>

### Table S4. Linear regression of multivariate trait distance versus phylogenetic distance by clade.
<small>
<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> term </th>
   <th style="text-align:right;"> coefficient </th>
   <th style="text-align:right;"> std.err </th>
   <th style="text-align:right;"> <i>t</i> </th>
   <th style="text-align:right;"> <i>p</i> </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> intercept </td>
   <td style="text-align:right;"> 1.03 </td>
   <td style="text-align:right;"> 0.17 </td>
   <td style="text-align:right;"> 6.22 </td>
   <td style="text-align:right;"> <10<sup>-5</sup> </td>
  </tr>
  <tr>
  <td style="text-align:left;"> <i>D<sub>g</sub></i> </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 6.46 </td>
   <td style="text-align:right;"> <10<sup>-5</sup> </td>
  </tr>
  <tr>
  <td style="text-align:left;"> clade (<i>P. syringae</i>) </td>
   <td style="text-align:right;"> 0.88 </td>
   <td style="text-align:right;"> 0.15 </td>
   <td style="text-align:right;"> 5.76 </td>
   <td style="text-align:right;"> <10<sup>-5</sup> </td>
  </tr>
</tbody>
</table>
</small>