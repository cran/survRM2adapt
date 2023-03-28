
# survRM2adapt: Flexible and Coherent Test/Estimation Procedure Based on Restricted Mean Survival Times

The survRM2adapt package includes the follwoing two main functions – 1)
rmst2adapt and 2)ltrmst2adapt.

1)  The rmst2adapt function estimates the restricted mean survival time
    (RMST) with a time window \[0, tau\] (ie, the standard RMST), where
    tau is selected during the procedure from among the time points
    provided by the user. For details, see Horiguchi et al. (2018) \[1\]

2)  The ltrmst2adapt function estimates the RMST with a time window
    \[tau1, tau2\] (ie, the long-term RMST), where tau1 is selected
    during the procedure, same as above. For details, see Horiguchi et
    al. (2023) \[2\]

## Installation

To install the package from the CRAN:

``` r
install.packages("survRM2adapt")
```

To install the most updated version:

``` r
devtools::install_github("uno1lab/survRM2adapt")
```

## Usage

Please see the package vignette on the website
<https://www.uno1lab.com/survRM2adapt/>

## Reference

\[1\] Horiguchi M, Cronin A, Takeuchi M, Uno H. A flexible and coherent
test/estimation procedure based on restricted mean survival times for
censored time-to-event data in randomized clinical trials. Statistics in
Medicine 2018;37(15):2307-2320. <https://doi.org/doi:10.1002/sim.7661>.

\[2\] Horiguchi M, Tian L, Uno H. On assessing survival benefit of
immunotherapy using long-term restricted mean survival time. Statistics
in Medicine 2023;42(8):1139-1155.
<https://doi.org/doi:10.1002/sim.9662>.
