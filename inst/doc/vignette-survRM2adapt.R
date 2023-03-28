## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  install.packages("survRM2adapt")

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  install.packages("devtools") #-- if the devtools package has not been installed
#  devtools::install_github("uno1lab/survRM2adapt")

## ---- echo=TRUE, eval=TRUE, message=FALSE-------------------------------------
library(survRM2adapt)

nrow(cm214_pfs)

cm214_pfs[1:10,]

## ---- echo=FALSE, eval=TRUE, fig.height=6, fig.width=6------------------------
library(survival)
plot(survfit(Surv(time, status)~arm, data=cm214_pfs), col=c("blue","red"), lwd=2, mark.time=F, xlab="Time (month)", ylab="Probability")
legend("bottomleft", c("Placebo (arm=0)","Treatment (arm=1)"), col=c("blue","red"), lwd=2)

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
ltrmst2adapt(indata=cm214_pfs, tau1=3, tau2=10, seed=123)

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
ltrmst2adapt(indata=cm214_pfs, tau1=c(0,1,2,3), tau2=10, seed=123)

