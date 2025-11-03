# Wildfire Occurrence Probability Modeling in Villaflores Municipality, Chiapas
Augusto Pablo Salonio Carbó
2025-09-16

The municipality of Villaflores, Chiapas, located within a high-density
wildfire cluster in Mexico, has faced an increase in both the frequency
and burned area of wildfires between 2019 and 2024, surpassing the
averages for 2012–2024, particularly in coniferous and montane cloud
forests. This study modeled the probability of wildfire occurrence
(2012–2024) using logistic regression, incorporating topo-climatic
predictors (e.g., wind speed), landscape predictors (e.g., green biomass
density, forest–agriculture interface), and demographic predictors
(e.g., population density). A wildfire probability map was generated,
and the relative importance of predictors was assessed.

The model, evaluated across 100 iterations, showed a useful
discriminatory capacity (mean AUC-ROC 0.791, standard deviation 0.018)
but a low positive predictive value (34%), suggesting a stronger
identification of areas potentially susceptible to wildfires than of
actual ignition points. The most influential predictors were
satellite-derived proxies of green biomass density (WDVI_Q3, 15% drop in
AUC-ROC when excluded) and herbaceous stratum density (MSAVI_Q1, 10%),
percentage of forest affected per pixel (Perc_Des), and wind speed
(Wind). Interactions between anthropogenic and landscape variables were
key to analyzing the two zones of highest probability: “West” (85%,
conifers, low accessibility) and “East” (95%, montane cloud forest,
atypical fires in February) near population centers. The 1 km²
resolution, which smooths local variations, and the exclusion of fire
records lacking burned-area polygons may potentially bias spatial
patterns.
