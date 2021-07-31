# GDFM-CHF_VaRES

Codes used in Hallin and Trucíos (2021)' paper (Econometrics and Statistics). The main codes are in this repository, however, additional matlab codes from the [MFE](https://www.kevinsheppard.com/code/matlab/mfe-toolbox/) Toolbox of Kevin Sheppard are also needed.

> Some Matlab codes are modifications from Mateo Barigozzi codes available on his [website](http://www.barigozzi.eu/Codes.html) and/or Mario Forni's codes kindly provided for him.

## Instructions

- The folder `cov_estimation` provides routines for computing the one-step-ahead portfolio returns based on different conditional covariance estimation procedures.
- The `VaRES_Forecast_Backtesting.R` function performs the bootstrap procedure to forecast both risk measures (using the result from `cov_estimation` folder) and also performs the backtesting exercise.


## Paper's abstract 
Beyond their importance from the regulatory policy point of view, Value-at-Risk (VaR) and Expected Shortfall (ES) play an important role in risk management, portfolio allocation, capital level requirements, trading systems, and hedging strategies. However, due to the curse of dimensionality, their accurate estimation and forecast in large portfolios is quite a challenge. To tackle this problem, two procedures are proposed. The first one is based on a filtered historical simulation method in which high-dimensional conditional covariance matrices are estimated via a general dynamic factor model with infinite-dimensional factor space and conditionally heteroscedastic factors; the other one is based on a residual-based bootstrap scheme. The two procedures are applied to a panel with concentration ratio close to one. Backtesting and scoring results indicate that both VaR and ES are accurately estimated under both methods, which both outperform the existing alternatives.



## References
- Trucíos, C., Mazzeu, J. H. G., Hallin, M., Hotta, L. K., Valls Pereira, P. L., & Zevallos, M. (2020). Forecasting conditional covariance matrices in high-dimensional time series: a general dynamic factor approach. Available at SSRN 3399782.
- Hallin, M. and Trucíos C. (2021). [Forecasting Value-at-Risk and Expected Shortfall in Large Portfolios: A General Dynamic Factor Model Approach.](https://www.sciencedirect.com/science/article/abs/pii/S2452306221000563?via%3Dihub)
