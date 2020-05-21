## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("gnomonicM")

## ---- eval=FALSE--------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("ejosymart/gnomonicM")

## ---- echo = FALSE, message = FALSE-------------------------------------------
library(gnomonicM)

## ----echo=TRUE----------------------------------------------------------------
model <- gnomonic(nInterval   = 7, 
                  eggDuration = 2, 
                  longevity   = 365, 
                  fecundity   = 200000, 
                  a_init      = 2)


## ----echo=TRUE----------------------------------------------------------------
modelAddInfo <- gnomonic(nInterval   = 7, 
                         eggDuration = 2, 
                         addInfo     = c(4, NA, NA, 40, NA, NA),
                         longevity   = 365, 
                         fecundity   = 200000, 
                         a_init      = 2)


## ----echo = TRUE--------------------------------------------------------------
#No additional information in the duration of gnomonic intervals
print(model)

plot(model)

#With additional information in the duration of gnomonic intervals
print(modelAddInfo)

plot(modelAddInfo, xlab = "My X label", ylab = "My Y label", cex = 3, bg = "blue")

## ----echo=TRUE----------------------------------------------------------------
modelUniform <- gnomonicStochastic(nInterval     = 7, 
                                   eggDuration   = 2,
                                   longevity     = 365,
                                   distr         = "uniform", 
                                   min_fecundity = 100000, 
                                   max_fecundity = 300000, 
                                   niter         = 1000, 
                                   a_init        = 2)


modelUniformAddInfo <- gnomonicStochastic(nInterval     = 7, 
                                          eggDuration   = 2,
                                          addInfo       = c(4, NA, NA, 40, NA, NA),
                                          longevity     = 365,
                                          distr         = "uniform", 
                                          min_fecundity = 100000, 
                                          max_fecundity = 300000, 
                                          niter         = 1000, 
                                           a_init        = 2)


modelNormal <- gnomonicStochastic(nInterval     = 7, 
                                  eggDuration   = 2,
                                  longevity     = 365,
                                  distr         = "normal", 
                                  fecundity     = 200000, 
                                  sd_fecundity  = 50000, 
                                  niter         = 1000, 
                                  a_init        = 2)


modelTriangle <- gnomonicStochastic(nInterval     = 7, 
                                    eggDuration   = 2,
                                    longevity     = 365,
                                    distr         = "triangle", 
                                    fecundity     = 200000,
                                    min_fecundity = 100000,
                                    max_fecundity = 300000,
                                    niter         = 1000, 
                                    a_init        = 2)

## ----echo = TRUE, results = 'hide'--------------------------------------------
#The results are not shown.
print(modelUniform)


## ----echo = TRUE--------------------------------------------------------------
plot(modelUniform, main = "Uniform distribution in MLF")
plot(modelUniformAddInfo, main = "Uniform distribution in MLF \nwith additional information in \nsome gnomonic intervals")
plot(modelNormal, main = "Normal distribution in MLF")
plot(modelTriangle, main = "Triangular distribution in MLF")

