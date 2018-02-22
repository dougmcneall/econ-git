# tol2018.R
# Reproducing the statistical analysis of Tol (2018) "The economic impacts
# of climate change"
# Doug McNeall dougmcneall@gmail.com

dat = read.table('https://raw.githubusercontent.com/dougmcneall/econ-git/master/tol2018.txt',
                 header = TRUE)


tol_fitted = read.table('tol2018_fitted_values.txt', header = TRUE)

fitorder = c('PWL', 'linear', 'parabolic','quadratic', 'exp',
             'double_exp', 'T^6', 'T^7')
attach(dat)
warming2 = warming^2
warming6 = warming^6
warming7 = warming^7
warmingexp = exp(warming)
warmingexpexp = exp(exp(warming))

# fit a linear model
fit0 = lm(impact ~ 0 + warming) 
# fit a parabola
fit1 = lm(impact ~ 0 + warming2) 
# fit a simple polynomial
fit2 = lm(impact ~ 0 + warming + warming2)
# fit exponential
fit3 = lm(impact ~ 0 + warmingexp)
# fit double exponential
fit4 = lm(impact ~ 0 + warmingexpexp)
# miss fit5 for labelling purposes :)
fit6 = lm(impact ~ 0 + warming2 + warming6)
fit7 = lm(impact ~ 0 + warming2 + warming7)
# piecewise linear fit
library(segmented)
fitseg = segmented(fit0, seg.Z = ~ 0 + warming)

warmingseq = seq(from = -1, to = 6.5, by = 0.1)

newdata = data.frame(warming = warmingseq,
                     warming2 = warmingseq^2,
                     warmingexp = exp(warmingseq),
                     warmingexpexp = exp(exp(warmingseq)),
                     warming6 = warmingseq^6,
                     warming7 = warmingseq^7
                     )

pred0= predict(fit0, newdata = newdata, se.fit = TRUE)
pred1 = predict(fit1, newdata = newdata, se.fit = TRUE)
pred2= predict(fit2, newdata = newdata, se.fit = TRUE)
pred3= predict(fit3, newdata = newdata, se.fit = TRUE)
pred4= predict(fit4, newdata = newdata, se.fit = TRUE)
pred6= predict(fit6, newdata = newdata, se.fit = TRUE)
pred7= predict(fit7, newdata = newdata, se.fit = TRUE)


plot(dat)
lines(newdata$warming, pred1$fit)
lines(newdata$warming, pred1$fit+(2*pred1$se.fit), lty = 'dashed')
lines(newdata$warming, pred1$fit-(2*pred1$se.fit), lty = 'dashed')

plot(dat, pch = 19, xlim = c(-1, 6), ylim = c(-9, 3))
points(warming, fit1$fitted.values, col = 'red', pch = 19)
points(warming, c(tsq, recursive = TRUE), col = 'pink', pch = 19)

plot(dat, pch = 19, xlim = c(-1, 6), ylim = c(-9, 3))
for(i in 1:ncol(tol_fitted)){
  points(warming, tol_fitted[,i], type = 'o')
  
}

predseg = predict(fitseg, newdata = newdata, se.fit = TRUE )

pdf('segmented_vs_parabolic.pdf')
plot(dat)
lines(newdata$warming, predseg$fit, lwd = 2)
lines(newdata$warming, predseg$fit+(2*predseg$se.fit), lty = 'dashed')
lines(newdata$warming, predseg$fit-(2*predseg$se.fit), lty = 'dashed')

lines(newdata$warming, pred1$fit, col = 'red', lwd = 2)
lines(newdata$warming, pred1$fit+(2*pred1$se.fit), lty = 'dashed', col = 'red')
lines(newdata$warming, pred1$fit-(2*pred1$se.fit), lty = 'dashed', col = 'red')
legend('topright', legend = c('piecewise linear','parabolic'), lwd = 2, col = c('black', 'red'))
abline(h = 0)
abline(v = 0)
dev.off()



# Plot all the predictions
plot(dat)
lines(newdata$warming, predseg$fit, col = 'black', lwd = 2)
lines(newdata$warming, pred0$fit, col = 'tomato2', lwd = 2)
lines(newdata$warming, pred1$fit, col = 'red', lwd = 2)
lines(newdata$warming, pred2$fit, col = 'blue', lwd = 2)
lines(newdata$warming, pred3$fit, col = 'green', lwd = 2)
lines(newdata$warming, pred4$fit, col = 'grey', lwd = 2)
lines(newdata$warming, pred6$fit, col = 'pink', lwd = 2)
lines(newdata$warming, pred7$fit, col = 'orange', lwd = 2)

# find the relative likelihood of the various fits
# https://en.wikipedia.org/wiki/Akaike_information_criterion#How_to_apply_AIC_in_practice
rl = function(x){
  out = exp( (min(x) - x)/2)
  out
}
# piecewise linear actually less likely using AIC?
aicvec = rl(c(AIC(fitseg), AIC(fit0), AIC(fit1), AIC(fit2), AIC(fit3), AIC(fit4), AIC(fit6), AIC(fit7)))
names(aicvec) = fitorder

bicvec = rl(c(BIC(fitseg),BIC(fit0), BIC(fit1), BIC(fit2), BIC(fit3),
              BIC(fit4), BIC(fit6), BIC(fit7)))
names(bicvec) = fitorder

# Relative likelihood ratio, calculated by aic, 
# expressed compared to total
raic = (aicvec/sum(aicvec)) * 100
# Relative likelihood ratio, calculated by bic,
# expressed compared to total
rbic = (bicvec/sum(bicvec)) * 100

# I think that this is how Tol is, in effect, doing his likelihood
# ratio test. This isn't appropriate for non-nested models.
# The key piece of information that I need is: are the PWL model
# and the parabolic model nested?
# What are nested models?
# https://www.theanalysisfactor.com/what-are-nested-models/
# https://en.wikipedia.org/wiki/Likelihood-ratio_test

# These give (kind of) similar answers to Tol in terms of
# likelihood ratio
expl = c(exp(logLik(fitseg)), exp(logLik(fit0)), exp(logLik(fit1)), exp(logLik(fit2)),
            exp(logLik(fit3)), exp(logLik(fit4)), exp(logLik(fit6)), exp(logLik(fit7)))

# relative likelihood ratio, expressed in a similar
# way to the spreadsheet
rexpl = (expl/sum(expl)) *100
names(rexpl) = fitorder

round(cbind(rexpl, raic, rbic))

# Check that AIC and loglikelihood are being used in the same way
# npar = length(coef(fit)) + 1, as I think AIC counts the estimated 
# uncertainty as a parameter. This shouldn't matter, as it's the 
# *difference* between number of parameters that counts in the AIC
# comparison.

npar = 2
AIC(fit0)
(-2*logLik(fit0)) + 2*npar

npar = 2
AIC(fit1)
(-2*logLik(fit1)) + 2*npar

npar = 3
AIC(fit2)
(-2*logLik(fit2)) + 2*npar

AIC(fitseg)
npar = 4
(-2*logLik(fitseg)) + 2*npar

# hand calulate AIC using only residuals and number
# of parameters in the model
calcAICdiff = function(fit){
  
  rss = sum(fit$residuals^2)
  # I think this should be longer (i.e. include variance), but shouldn't
  # matter if we are just comparing AIC in similar models
  # i.e., it's the difference in parameters that matters.
  k = length(coef(fit)) # + 1
  n = length(fit$residuals)
  
  aic = 2*k + n*(log(rss))
  aic
}

testAIC = rl(c(calcAICdiff(fitseg),
            calcAICdiff(fit0),
            calcAICdiff(fit1),
            calcAICdiff(fit2), 
            calcAICdiff(fit3), 
            calcAICdiff(fit4),
            calcAICdiff(fit6),
            calcAICdiff(fit7))
          )

# These should be the same
cbind(testAIC, aicvec)


# Can we reproduce the likelihood (and relative likelihood)
# Values from the spreadsheet?
toll = function(fit, rez = NULL){
  # this is my attempt at backing out the likelihood formula
  # from the spreadsheet, but it doesn't seem to produce the 
  # same results.
  if(is.null(rez)){
    rez = fit$residuals
  }
  
  rss = sum(rez^2)
  k = length(coef(fit)) 
  n = length(rez)
  ll = -(rss/2) /  (  (sqrt(rss/(n-k)))^2 - n*(log(sqrt(rss/(n-k)))))
  
  ll
}
tollvec = c(exp(toll(fitseg)), exp(toll(fit0)), exp(toll(fit1)), exp(toll(fit2)), 
            exp(toll(fit3)), exp(toll(fit4)), exp(toll(fit6)), exp(toll(fit7))) 

# relative likelihood
rtol = (tollvec/sum(tollvec))*100
names(rtol) = fitorder

#table of the relative likelihood measures
out = round(cbind(rtol, rexpl, raic, rbic))

colnames(out) = c('Tol RL', 'RL', 'AIC RL', 'BIC RL')
print(out)

# Received email from R. Tol 16/02/2018
#Likelihood for one observation = 1/sqrt(2pi)/sigma e^(-0.5((x-mu)/sigma)^2)

#Take log
#(x-mu)^2 is the squared residual
#Loglikelihood for all observations is thus proportional to the sum of squared residuals
#The maximum likelihood estimator of sigma^2 is the sum of squared residuals over the number of degrees of freedom
#With an uniformative prior, the likelihood of the data given the model (calculated above) is proportional to the likelihood of the model given the data. For a countable set of candidate models, you rescale to sum to one.

# Calculate sigma by hand, with info from Tol's email
testfun = function(fit, rez=NULL){
  
  if(is.null(rez)){
    rez = fit$residuals
  }
  
 # if(is.null(rez)){
  p = length(coef(fit))
#  }

  n = length(rez)
  rss = sum(rez^2)
  sigma = sqrt(rss/(n))
  #ll = log(1/(sqrt(2*pi)/sigma) * exp(-0.5*(sum(fit$residuals)/sigma)^2))
  #ll = log(1/sqrt(2*pi*sigma) * exp(-0.5*(sum(rez^2)/sigma^2)))
  #ll = log
  truell = -(n/2)*log(2*pi) - (n/2)*log(sigma^2) - (1/(2*sigma^2))*rss
  truell
}

rtestfun = (testfunvec/sum(testfunvec))*100
names(rtestfun) = fitorder

# Table of relative likelihoods.
RLtab = round(cbind(rtestfun, rtol, rexpl, raic, rbic))
rownames(RLtab) = fitorder
print(RLtab)

# Here's Tol's function working on the fits from Tol's own spreadsheet.
# This should reproduce Tol's results, apart from mistakenly swapped
# Parabolic and Quadratic labels.
# It doesn't, meaning I must have the calculation slighly wrong
tolfit = c(exp(toll(fitseg, rez = tol_fitted[,2] - impact)),
           exp(toll(fit0, rez = tol_fitted[,8] - impact)),
           exp(toll(fit1, rez = tol_fitted[,3] - impact)),
           exp(toll(fit2, rez = tol_fitted[,1] - impact)),
           exp(toll(fit3, rez = tol_fitted[,6] - impact)),
           exp(toll(fit4, rez = tol_fitted[,7] - impact)),
           exp(toll(fit6, rez = tol_fitted[,5] - impact)),
           exp(toll(fit7, rez = tol_fitted[,4] - impact))
)
rtolfit = (tolfit/sum(tolfit))*100

# Here's my first-principles likelihood calculation, working on
# Tol's fitted values from the spreadsheet. This DOES reproduce the
# results from the spreadsheet.
testfunvec_tol = c(exp(testfun(fitseg, rez = tol_fitted[,2] - impact)),
                   exp(testfun(fit0, rez = tol_fitted[,8] - impact)),
                   exp(testfun(fit1, rez = tol_fitted[,3] - impact)),
                   exp(testfun(fit2, rez = tol_fitted[,1] - impact)),
                   exp(testfun(fit3, rez = tol_fitted[,6] - impact)),
                   exp(testfun(fit4, rez = tol_fitted[,7] - impact)),
                   exp(testfun(fit6, rez = tol_fitted[,5] - impact)),
                   exp(testfun(fit7, rez = tol_fitted[,4] - impact))
)

rtestfun_tol = (testfunvec_tol/sum(testfunvec_tol))*100

#Same answers as Tol, but with parabolic/quadratic reversed?
# HA! Finally found the mistake. Checked on the spreadsheet - Tol has mistakenly reversed
# "parabolic" and "quadratic" labels in the final table. 
# He calls cell AJ9 = W4 Quadratic (actually T^2) and AJ3 = S4 Parabolic (actually T +T^2)
RL_tolfits = round(cbind(rtolfit, rtestfun_tol ))
rownames(RL_tolfits) = fitorder
print(RL_tolfits)


# AIC with from-first-principles likelihood
fp_AIC = function(fit, rez=NULL){
  
  if(is.null(rez)){
    rez = fit$residuals
  }
  
  L = exp(testfun(fit, rez = rez))
  k = length(coef(fit)) +1

  aic = (2*k)  - 2*log(L) 
  aic
  
}

# This is what happens when you apply AIC with from-first-principles
# likelihood to Tol's fits
rl_AIC_fp_tolfits = rl(c(fp_AIC(fitseg, rez = tol_fitted[,2] - impact),
                   fp_AIC(fit0, rez = tol_fitted[,8] - impact),
                   fp_AIC(fit1, rez = tol_fitted[,3] - impact),
                   fp_AIC(fit2, rez = tol_fitted[,1] - impact),
                   fp_AIC(fit3, rez = tol_fitted[,6] - impact),
                   fp_AIC(fit4, rez = tol_fitted[,7] - impact),
                   fp_AIC(fit6, rez = tol_fitted[,5] - impact),
                   fp_AIC(fit7, rez = tol_fitted[,4] - impact)
)
)
round((rl_AIC_fp_tolfits/sum(rl_AIC_fp_tolfits))*100)

# This is what happens when you apply AIC with from-first-principles
# likelihood to the correct fits
rl_AIC_fp = rl(c(fp_AIC(fitseg),
                 fp_AIC(fit0),
                 fp_AIC(fit1),
                 fp_AIC(fit2),
                 fp_AIC(fit3),
                 fp_AIC(fit4),
                 fp_AIC(fit6),
                 fp_AIC(fit7)
)
)

round((rl_AIC_fp/sum(rl_AIC_fp))*100)

# First principles likelihood using Tol's RL method 
  

# --------------------------------------------------------------------
# What do the fits look like when removing Tol2002?
# --------------------------------------------------------------------

dat.trunc = dat[-15, ]
attach(dat.trunc)
plot(dat.trunc, ylim = c(-7,5))
warming2 = warming^2

# Initial linear fit
fit0.trunc = lm(impact~warming, data = dat.trunc)

fit1.trunc = lm(impact~warming2, data = dat.trunc)
pred1.trunc = predict(fit1.trunc, newdata = newdata, se.fit = TRUE)

fitseg.trunc =segmented(fit0.trunc, seg.Z = ~warming)
predseg.trunc = predict(fitseg.trunc, newdata = newdata, se.fit = TRUE )

pdf(file = 'excludeTol2002.pdf')
par(las = 1)
plot(dat.trunc, ylim = c(-7,3), pch = 19)
# without Tol2002
lines(newdata$warming, pred1.trunc$fit, lwd = 2)
lines(newdata$warming, predseg.trunc$fit, lwd = 2)
# with Tol2002
lines(newdata$warming, pred1$fit, col = 'red', lwd = 2)
lines(newdata$warming, predseg$fit, col = 'red', lwd = 2)
points(dat[15,], col = 'red', pch = 19)
abline(h = 0)
abline(v = 0)
legend('topright', legend = c('with Tol (2002)', 'without Tol (2002)'), col = c('red', 'black'), lwd = 2)
dev.off()


