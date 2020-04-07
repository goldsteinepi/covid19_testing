#################
# Performance of testing in COVID-19
# Citation: Goldstein ND, Burstyn I. On the importance of early testing even when imperfect in a pandemic such as COVID-19. Manuscript in preparation.
# 3/17/20 -- Igor Burstyn and Neal Goldstein
#################

#define parameter ranges
sn_pcr = seq(0.6, 0.9, by=0.01)
sp_pcr = seq(0.9, 0.99, by=0.01)

pop = 1000 #population size
prev = seq(0, 0.5, by=0.05)  #COVID-19 prevalance

fn_hi = NA
fn_lo = NA
fp_hi = NA
fp_lo = NA

#tp = NA 
#tn = NA

for (i in 1:length(prev))
{
  true_pos = prev[i] * pop
  true_neg = (1-prev[i]) * pop
  test_pos_true_pos = round(sn_pcr * true_pos)
  test_neg_true_pos = true_pos - test_pos_true_pos
  test_neg_true_neg = round(sp_pcr * true_neg)
  test_pos_true_neg = true_neg - test_neg_true_neg
  
  fn_hi = c(na.omit(fn_hi), max(test_neg_true_pos))
  fn_lo = c(na.omit(fn_lo), min(test_neg_true_pos))
  fp_hi = c(na.omit(fp_hi), max(test_pos_true_neg))
  fp_lo = c(na.omit(fp_lo), min(test_pos_true_neg))
}
rm(i)

true_pos = prev * pop
true_neg = (1-prev) * pop
test_pos_true_pos = round(sn_pcr * true_pos)
#test_pos_true_pos = round(sn_pcr[4] * true_pos)
test_neg_true_pos = true_pos - test_pos_true_pos
test_neg_true_neg = round(sp_pcr * true_neg)
test_pos_true_neg = true_neg - test_neg_true_neg

plot(x=prev, y=fn_hi, type="n", ylim=c(0,pop/2), xlab="Hypothetical COVID-19 prevalence", ylab="Count of tests")
polygon(x=c(min(prev),max(prev),max(prev),min(prev)), y=c(min(fn_hi),max(fn_hi),max(fn_lo),min(fn_lo)), col=rgb(red=0.8,blue=0.8,green=0.8,alpha=1), border="NA")
polygon(x=c(min(prev),max(prev),max(prev),min(prev)), y=c(max(fp_hi),min(fp_hi),min(fp_lo),max(fp_lo)), col=rgb(red=0.2,blue=0.2,green=0.2,alpha=.75), border="NA")
lines(x=prev, y=prev * pop, lwd=3)
legend("topleft", legend=c("False positives","False negatives","True prevalence"), fill=c(rgb(red=0.2,blue=0.2,green=0.2,alpha=.75),rgb(red=0.8,blue=0.8,green=0.8,alpha=1),NA), border=c("black","black","white"), lty=c(NA,NA,1), cex=0.8)

#5% prevalence
i=2
true_pos = prev[i] * pop
true_neg = (1-prev[i]) * pop
test_pos_true_pos = round(sn_pcr * true_pos)
test_neg_true_pos = true_pos - test_pos_true_pos
test_neg_true_neg = round(sp_pcr * true_neg)
test_pos_true_neg = true_neg - test_neg_true_neg

range(test_neg_true_pos) #fn
range(test_pos_true_neg) #fp

#25% prevalence
i=6
true_pos = prev[i] * pop
true_neg = (1-prev[i]) * pop
test_pos_true_pos = round(sn_pcr * true_pos)
test_neg_true_pos = true_pos - test_pos_true_pos
test_neg_true_neg = round(sp_pcr * true_neg)
test_pos_true_neg = true_neg - test_neg_true_neg

range(test_neg_true_pos) #fn
range(test_pos_true_neg) #fp

