# R_IC user guide

Implementation of the Index of Connectivity in R language

Author details: Tommaso Baggio, Lorenzo Martini, Loris Torresani.
Script and data info: This script performs the Index of Connectivity (Cavalli et al., 2013).
Copyright statement: This script is the product of the work of Tommaso Baggio, Lorenzo Martini and Loris Torresani.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

The release performs the Index of Connectivity (IC) starting from the version developed by Cavalli et al. (2013). We implemented IC in the R programming language to make it more accessible and versatile for final users. Furthermore, beyond the traditional implementation of computing IC, we developed two algorithms to compute IC recursively (batch process) and a more flexible version, allowing user-defined weighting factor and the possibility to perform an analysis of the outputs along a given longitudinal profile.

### Dependencies
In order to run the three scripts provided in the release, TauDEM has to be installed in order to successively call the executables. You can download and install the latest version of the complete windows installer from https://hydrology.usu.edu/taudem/taudem5/downloads.html To develop the three scripts, version number 5.3.7 was used, however, R_IC is expected to work also with earlier versions.
Since the script is written in R language, R software is required https://cran.r-project.org/bin/windows/base/ The adopted version is the 4.0.3. To compute operations, we further used built-in functions from different R packages (before installing and loading them remember to install Rtools 4 from https://cran.r-project.org/bin/windows/Rtools/rtools40.html). The necessary packages are raster (v 3.4.5), shapefiles (v 0.7), rgdal (v 1.5.19), cowplot (v 1.1.1) and ggplot2 (v 3.3.2). The reported versions are those used to develop the three scripts. However, we expect that earlier or updated versions of R and packages will be suitable as well.

### Usage 
Three scripts are available in the repository of R_IC. The scripts are structured in the same way: the first lines are reserved for the input files, followed by operations to check the input data. If the data are correct, the script continues with the core algorithm to compute IC and to save the outputs. The common inputs of the three scripts are:
-	Definition of the working directory: directory where input files are located and where output files will be saved. Remember to not use the symbol “/” as last character.
-	Names of the input files: the raster file representing the depitted DEM (the scripts do not perform any type of control on the input raster in terms of hydrological consistency) and the shapefile representing the target of the IC analysis. Note that these files need to be in the same reference coordinate system. The formats of the input files are those accepted by the functions raster and readOGR of the package raster and rgdal, respectively. In the script, batch_R_IC remember to not specify the extension of the shapefile (as reported in the training data).
-	Parameters representing (i) the size of the moving window to derive the roughness map and then the weighting factor expressed as cell number (only odd numbers greater than 1 are acceptable) and (ii) a flag to save intermediate outputs (flag = 1) or to keep only essential outputs as IC, Ddown, Dup and roughness (flag=0).



#################################################################
                      Details on statistical tests
#################################################################

Having continuous annual data from 1991 to 2021 across ~220 catchments elevates your study from a simple static comparison to a 31-year Spatio-Temporal Panel (Longitudinal) Analysis (≈6,820 catchment-year observations).
Because you have repeated measurements on the same catchments over 31 years, standard textbook ANOVA and simple regressions cannot be used directly without violating the assumption of independent observations (consecutive years within a catchment are temporally correlated).
Below is a structured breakdown of the most scientifically rigorous statistical tests, models, and visualizations suited for your specific objectives.

Module 1: How IC Varies Over Time (1991–2021 Trends)
1. Mann-Kendall Trend Test with Sen’s Slope Estimator (Catchment & Regional Scale)
What it tests: Whether there is a statistically significant monotonic upward (degradation / increasing connectivity) or downward (greening / trapping) trend over the 31-year period.
Why it fits: Non-parametric, robust to missing years, non-normal distributions, and extreme climate anomalies (e.g., severe flood or drought years).
Two levels of application:
Catchment-by-Catchment: Computes a Sen’s slope (ΔIC/year) and p-value for every single basin.
Regional Kendall Test (Hirsch-Slack): Tests whether Peninsular India as a whole experienced a uniform subcontinent-wide trend.
Recommended Plot:
A spatial map of Peninsular India where catchment outlets/boundaries are colored by trend direction and significance (Green = significant decrease, Red = significant increase, Gray = non-significant).
A continuous 31-year regional time-series line plot showing the mean annual IC curve flanked by a 95% confidence envelope.
2. Non-Parametric Breakpoint / Change-Point Detection (Pettitt’s Test or PELT)
What it tests: Did connectivity change gradually, or was there an abrupt shift (regime shift) around a specific year (e.g., post-2000 afforestation policies, major dam construction, or severe drought cycles)?
Recommended Plot: Cumulative deviation curve showing the critical year where the slope of connectivity inverted or stepped.

Module 2: Variation Across Geological and Climate Zones (ANOVA Equivalents)
1. Linear Mixed-Effects Model (LMM / Hierarchical ANOVA) — Primary Recommendation
The Formulation:
IC_it = β_0 + β_1(Geology_i) + β_2 (Climate_i) + β_3 (Regime_i) + β_4 (Year_t) + u_i + ε_it
where u_i is a catchment-specific random intercept.

What it tests (Type-III ANOVA table):
Main Effects: Does IC differ significantly between Geological classes (Igneous, Metamorphic, Sedimentary, Unconsolidated) and Climate zones after controlling for temporal trends?
Interaction Effects (Geology × Year): Does connectivity change faster over time on Sedimentary vs. Igneous bedrock?
Why it fits: It functions like a classic multi-factor ANOVA, but properly accounts for repeated measures within catchments without deflating p-values.

2. Repeated-Measures Permutation ANOVA (PERMANOVA)
What it tests: A distribution-free alternative that evaluates differences between grouped centroids (Geology, Climate) through 9,999 randomized permutations of the data matrix.
Recommended Plot:
A multi-line annual time series (1991–2021) with 4 colored ribbons (one for each geology class) showing how the gap between lithologies evolved over three decades.
Faceted boxplots comparing annual variability across Climate zones.
Module 3: Relation with the Annual Binary Dataset
(Whether this binary variable represents Drought vs. Normal years, El Niño vs. La Niña, or Pre- vs. Post-dam installation):
1. Paired Wilcoxon Signed-Rank Test (Within-Basin Comparison)
What it tests: For each catchment, compute its average IC during "State 0" years versus "State 1" years. Tests whether switching to State 1 induces a consistent, significant upward or downward shift in connectivity across all basins.
Recommended Plot:
Paired Slope / Dumbbell Plot: Each catchment is represented by two connected points (State 0 → State 1) to visualize whether basins uniformly shift toward higher or lower connectivity.
2. Event-Study / Difference-in-Differences (DiD) Framework
What it tests: If the binary event occurred at different times in different catchments (e.g., dam commissioning or local policy changes), this model estimates the causal impact of the event while removing regional climate background noise.
Module 4: Relation with Continuous Annual Datasets (Rainfall, Temperature, NDVI, etc.)
1. Panel Fixed-Effects ("Within-Catchment") Regression
What it tests:
(IC_it - IC_i) = β_1 (X_(1, it) - X_(1, i)) + β_2 (X_(2, it) - X_(2, i)) + ... + ε_it
 
Why it fits: By subtracting each basin's 31-year mean (IC_i), this method completely eliminates all static confounders (elevation, basin area, slope, lithology). The resulting coefficients (β_1, β_2) tell you strictly: "When annual rainfall or temperature deviates from normal by +10%, by how much does connectivity shift?"
2. Repeated-Measures Correlation (rmcorr)
What it tests: Evaluates the common within-individual linear association between continuous annual climate drivers and IC without committing the "ecological fallacy" of pooling longitudinal data into a simple scatter plot.
3. Lagged Cross-Correlation Analysis (t, t-1, t-2)
What it tests: Does catchment connectivity react instantaneously to annual climate drivers, or is there a 1- to 2-year lag (e.g., rainfall in year t−1 driving vegetation growth in year t, which controls connectivity)?
Recommended Plot: Cross-correlation correlogram (r vs. lag years −3 to +3).
Panel scatter plots with catchment-level regression slopes.