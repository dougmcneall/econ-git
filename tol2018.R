# tol2018.R

dat <- read.table('tol2018.txt', header = TRUE)

attach(dat)
warming2 = warming^2

# fit a linear model (useful for later)
fit0 = lm(impact ~ warming) 

# fit a parabola
fit1 = lm(impact ~ warming2)
warmingseq = seq(from = -1, to = 6, by = 0.1)
newdata = data.frame(warming = warmingseq,
                     warming2 = warmingseq^2)

pred1 = predict(fit1, newdata = newdata, se.fit = TRUE)
plot(dat)
lines(newdata$warming, pred1$fit)
lines(newdata$warming, pred1$fit+(2*pred1$se.fit), lty = 'dashed')
lines(newdata$warming, pred1$fit-(2*pred1$se.fit), lty = 'dashed')

# fit a simple polynomial
fit2 = lm(impact ~ warming + warming2)

# piecewise linear
library(segmented)
fitseg = segmented(fit0, seg.Z = ~warming)
predseg = predict(fitseg, newdata = newdata, se.fit = TRUE )

plot(dat)
lines(newdata$warming, predseg$fit)
lines(newdata$warming, predseg$fit+(2*predseg$se.fit), lty = 'dashed')
lines(newdata$warming, predseg$fit-(2*predseg$se.fit), lty = 'dashed')
abline(h = 0)
abline(v = 0)


# find the relative likelihood of the various fits
# https://en.wikipedia.org/wiki/Akaike_information_criterion#How_to_apply_AIC_in_practice
rl = function(x){
  out = exp( (min(x) - x)/2)
  out
}
# piecewise linear actually less likely using AIC?
rl(c(AIC(fitseg), AIC(fit0), AIC(fit1), AIC(fit2)))

# a look at the residuals
sqrt(mean(fit0$residuals^2))
sqrt(mean(fit1$residuals^2))
sqrt(mean(fit2$residuals^2))
sqrt(mean(fitseg$residuals^2))

plot(resid(fit1))
abline(h = 0)
points(resid(fitseg), col = 'red')

# I think that this is how Tol is, in effect, doing his likelihood
# ratio test. This isn't appropriate for non-nested models.
# The key piece of information that I need is: are the PWL model
# and the parabolic model nested?
# What are nested models?
# https://www.theanalysisfactor.com/what-are-nested-models/
# https://en.wikipedia.org/wiki/Likelihood-ratio_test

# These give similar answers to Tol in terms of
# likelihood ratio
exp(logLik(fit1))/exp(logLik(fitseg))
exp(logLik(fit0))/exp(logLik(fitseg))

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



# What do the fits look like when removing Tol2002?
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


