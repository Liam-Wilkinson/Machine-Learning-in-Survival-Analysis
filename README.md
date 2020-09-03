# Simulation Study for Machine Learning Methods in Survival Analysis

Simulate survival times under a Weibull proportional hazards model with either linear or linear and non-linear effects of the covariates on the hazard and varying levels of censoring.

Splits simulated datasets into training and test datasets.

Fits a random forest; support vector machine; Cox model with elastic-net penalty and and penalized splines (flexible parametric) model to the training datasets.

Records how long the various methods take to fit.

Uses the fitted models from the training datasets to predict risk scores in the test datasets. 

Computes the concordance index comparing predicted risk scores in the test datasets with the observed censoring status and survival times.

