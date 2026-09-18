# Bayesian-Inference-and-Data-Assimilation-for-Bivariate-Log-Price-Models
This project applies Bayesian inference to weekly average log-prices for two agricultural commodities observed over 30 weeks. A bivariate normal model with a fixed correlation structure is used to estimate the unknown mean log-price of one commodity and a shared precision parameter, using a conditional normal prior for the mean and a Gamma prior for the precision.

The analysis derives marginal and conditional distributions, constructs the likelihood and joint posterior, and obtains the full conditional posterior distributions needed for Gibbs sampling. A MATLAB implementation alternates between normal and Gamma draws, discarding 5,000 burn-in iterations and retaining 10,000 samples for posterior analysis.

Posterior uncertainty is summarised through estimated means, standard deviations, and 95% credible intervals. Trace plots and posterior histograms are used to visually assess chain mixing and explore the sampled distributions.

## Project Report
[View full report](DataAssimilation_1.pdf)

## License
This project is provided for viewing and evaluation purposes only.
Reuse, modification, or distribution is not permitted without explicit permission.
