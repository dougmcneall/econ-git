# tol2018.R

# NOTE need to refactor to remove intercepts!
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
fit2 = lm(impact ~ 0+ warming + warming2)
# fit exponential
fit3 = lm(impact ~ 0 + warmingexp)
# fit double exponential
fit4 = lm(impact ~ 0 + warmingexpexp)
# miss fit5 for labelling purposes :)
fit6 = lm(impact ~ 0 + warming2 + warming6)
fit7 = lm(impact ~ 0 + warming2 + warming7)

warmingseq = seq(from = -1, to = 6.5, by = 0.1)

newdata = data.frame(warming = warmingseq,
                     warming2 = warmingseq^2,
                     warmingexp = exp(warmingseq),
                     warmingexpexp = exp(exp(warmingseq)),
                     warming6 = warmingseq^6,
                     warming7 = warmingseq^7
                     )

pred1 = predict(fit1, newdata = newdata, se.fit = TRUE)
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

# piecewise linear
library(segmented)
fitseg = segmented(fit0, seg.Z = ~0 + warming)
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

pred0= predict(fit0, newdata = newdata, se.fit = TRUE )
pred2= predict(fit2, newdata = newdata, se.fit = TRUE )
pred3= predict(fit3, newdata = newdata, se.fit = TRUE )
pred4= predict(fit4, newdata = newdata, se.fit = TRUE )
pred6= predict(fit6, newdata = newdata, se.fit = TRUE )
pred7= predict(fit7, newdata = newdata, se.fit = TRUE )

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

# Relative likelihood ratio, calculated by aic, expressed in the same way
raic = (aicvec/sum(aicvec)) * 100

cbind(rexpl, raic)

# Check that AIC and loglikelihood are being used in the same way
npar = 3
AIC(fit0)
(-2*logLik(fit0)) + 2*npar

npar = 3
AIC(fit1)
(-2*logLik(fit1)) + 2*npar

npar = 4
AIC(fit2)
(-2*logLik(fit2)) + 2*npar

AIC(fitseg)
npar = 5
(-2*logLik(fitseg)) + 2*npar

# Looks like the PWL fit has 5 parameters, which is why it is being more heavily
# penalised by AIC measure.

# hand calulate AIC using only residuals and number
# of parameters in the model
calcAICdiff = function(fit){
  
  rss = sum(fit$residuals^2)
  # I think this should be longer (i.e. include variance), but shouldn't
  # matter if we are just comparing AIC in similar models
  # i.e., it's the difference in parameters that matters.
  k = length(coef(fit)) 
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
toll = function(fit){
  # this is my attempt at backing out the likelihood formula
  # from the spreadsheet, but it doesn't seem to produce the 
  # same results.
  
  rss = sum(fit$residuals^2)
  k = length(coef(fit)) -1 # don't count the intercept, apparently
  n = length(fit$residuals)
  ll = -(rss/2) /  (  (sqrt(rss/(n-k)))^2 - n*(log(sqrt(rss/(n-k)))))
  
  ll
}
tollvec = c(exp(toll(fitseg)), exp(toll(fit0)), exp(toll(fit1)), exp(toll(fit2)), 
            exp(toll(fit3)), exp(toll(fit4)), exp(toll(fit6)), exp(toll(fit7))) 

# relative likelihood
rtol = (tollvec/sum(tollvec))*100
names(rtol) = fitorder

#table of the relative likelihood measures
out = cbind(round(rtol, 3), round(raic, 3), round(rexpl, 3))

colnames(out) = c('Tol_RL', 'AIC', 'R_RL')
print(out)



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


