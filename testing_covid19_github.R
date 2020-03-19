#################
# Performance of testing in COVID-19
# Citation: Goldstein ND, Burstyn I. On the importance of early testing even when imperfect in a pandemic such as COVID-19. Manuscript in preparation.
# 3/17/20 -- Igor Burstyn and Neal Goldstein
#################


### SIMULATION PARAMETERS ###

#reproducibility
set.seed(777)

#simulation size
n=5000


### DATA FOR  TEST ACCURACY ###

#range of PCR values
sp_pcr = 0.99
sn_pcr = c(0.60,0.70,0.80,0.90)

#sn_sero/sp_sero from gold standard serologic data https://www.ncbi.nlm.nih.gov/pubmed/32104917
sn_sero=rbeta(n, 352, (397-352))
sp_sero=rbeta(n, 116, 12) #128-116

sn_sero_mean=mean(sn_sero)
sp_sero_mean=mean(sp_sero)


### DATA FOR PREV and TESTING DEPLOYMENT ###

#as of March 19, 2020

#https://www.cdc.gov/coronavirus/2019-ncov/cases-in-us.html
cases=10442

#https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/testing-in-us.html
tested=4484+33340
cases/tested


### PLOT ###

pr=seq(0.01, .5, length.out=20)
pospredval=sn_pcr[1]*pr/(sn_pcr[1]*pr+(1-sp_pcr)*(1-pr))
negpredval=sp_pcr*(1-pr)/((1-sn_pcr[1])*pr+sp_pcr*(1-pr))

plot(100*pr, 100*(1-pospredval), type="b", col=gray.colors(6)[2], lwd=2, ylim=c(0,100), xlim=c(0,50), 
     ylab="Predictive value of the test (%)", xlab="Hypothetical COVID-19 prevalence (%)",  pch=16) 
lines(100*pr, 100*(1-negpredval), col=gray.colors(6)[2], lwd=2, type="b", pch=1)

pospredval=sn_pcr[2]*pr/(sn_pcr[2]*pr+(1-sp_pcr)*(1-pr))
negpredval=sp_pcr*(1-pr)/((1-sn_pcr[2])*pr+sp_pcr*(1-pr))
lines(100*pr, 100*(1-pospredval), col=gray.colors(6)[3], lwd=2, type="b", pch=16)
lines(100*pr, 100*(1-negpredval), col=gray.colors(6)[3], lwd=2, type="b", pch=1)

pospredval=sn_pcr[3]*pr/(sn_pcr[3]*pr+(1-sp_pcr)*(1-pr))
negpredval=sp_pcr*(1-pr)/((1-sn_pcr[3])*pr+sp_pcr*(1-pr))
lines(100*pr, 100*(1-pospredval), col=gray.colors(6)[4], lwd=2, type="b", pch=16)
lines(100*pr, 100*(1-negpredval), col=gray.colors(6)[4], lwd=2, type="b", pch=1)

pospredval=sn_pcr[4]*pr/(sn_pcr[4]*pr+(1-sp_pcr)*(1-pr))
negpredval=sp_pcr*(1-pr)/((1-sn_pcr[4])*pr+sp_pcr*(1-pr))
lines(100*pr, 100*(1-pospredval), col=gray.colors(6)[5], lwd=2, type="b", pch=16)
lines(100*pr, 100*(1-negpredval), col=gray.colors(6)[5], lwd=2, type="b", pch=1)

pospredval=sn_sero_mean*pr/(sn_sero_mean*pr+(1-sp_sero_mean)*(1-pr))
negpredval=sp_sero_mean*(1-pr)/((1-sn_sero_mean)*pr+sp_sero_mean*(1-pr))
lines(100*pr, 100*(1-pospredval), col="black", lwd=2, type="b", pch=16)
lines(100*pr, 100*(1-negpredval), col="black", lwd=2, type="b", pch=1)

legend(x=20, y=95, legend=c("Not infected when test positive", "Infected when test negative"), col="black", pch=c(16,1), cex=0.8, title="True infection status")
legend(x=42, y=95, legend=c("60%", "70%", "80%", "90%"), col=c(gray.colors(6)[2],gray.colors(6)[3],gray.colors(6)[4],gray.colors(6)[5]), lty=1, lwd=2, cex=0.8, title="Sensitivity")


### HYPOTHETHICAL ACCURACY ###

pop = 100000 #population size
prev = 0.28  #COVID-19 prevalance

true_pos = prev * pop
true_neg = (1-prev) * pop
test_pos_true_pos = round(sn_pcr[1] * true_pos)
test_pos_true_pos = round(sn_pcr[4] * true_pos)
test_neg_true_pos = true_pos - test_pos_true_pos
test_neg_true_neg = round(sp_pcr * true_neg)
test_pos_true_neg = true_neg - test_neg_true_neg


# ### DISTRIBUTIONS of PREDICTED COUNTS UNDER PCR TESTING ###
# 
# #3x3 plotting window
# par(mfrow=c(3,3))
# 
# #high prevalence: at equillibrium of 50% (China, Italy) for R=2
# #pr=rbeta(n, 9, 1) #90+/-10%%
# pr=rbeta(n, 100, 100) #50+/-3%%
# mean(pr)
# sd(pr)
# ppv=sn_pcr*pr/(sn_pcr*pr+(1-sp_pcr)*(1-pr))
# npv=sp_pcr*(1-pr)/((1-sn_pcr)*pr+sp_pcr*(1-pr))
# plot(ppv, npv, main="Prevalence at equilibirum 50%+/3% for R=2; medians indicated by solid lines", xlab="Pr(infected if test positive)", ylab="Pr(not infected if test negative)", col="grey")
# abline(v=median(ppv))
# abline(h=median(npv))
# TP=testpos*ppv
# hist(TP, main=bquote("actual infected among "~ .(testpos)~ "tested positive in USA;"~"median"~.(round(median(TP,2)))), xlab="true positive")
# abline(v=median(TP), col="red")
# #	TP=testposp_pcrA*ppv
# #	hist(TP, main="actual infected vs among tested positive in PA")
# #	abline(v=median(TP), col="red")
# FN=(testedUSA-testpos)*(1-npv)
# hist(FN, main=bquote("infected among"~.(testedUSA-testpos)~ "who tested negative in USA;"~"median"~.(round(median(FN,2)))), xlab="false negative")
# abline(v=median(FN), col="red")
# 
# #low prevalence (USA)
# pr=rbeta(n, 1, 9) #10+/-10%%
# mean(pr)
# sd(pr)
# ppv=sn_pcr*pr/(sn_pcr*pr+(1-sp_pcr)*(1-pr))
# npv=sp_pcr*(1-pr)/((1-sn_pcr)*pr+sp_pcr*(1-pr))
# plot(ppv, npv, main="Low 10%+/10% prevalence (e.g. USA); medians indicated by solid lines", xlab="Pr(infected if test positive)", ylab="Pr(not infected if test negative)", col="grey")
# abline(v=median(ppv))
# abline(h=median(npv))
# TP=testpos*ppv
# hist(TP, main=bquote("actual infected among "~ .(testpos)~ "tested positive in USA;"~"median"~.(round(median(TP,2)))), xlab="true positive")
# abline(v=median(TP), col="red")
# #	TP=testposp_pcrA*ppv
# #	hist(TP, main="actual infected among 22 tested positive in PA")
# #	abline(v=median(TP), col="red")
# FN=(testedUSA-testpos)*(1-npv)
# hist(FN, main=bquote("infected among"~.(testedUSA-testpos)~ "who tested negative in USA;"~"median"~.(round(median(FN,2)))), xlab="false negative")
# abline(v=median(FN), col="red")
# 
# #very low prevalence (USA)
# pr=rbeta(n, 1, 99) #1+/-1%%
# mean(pr)
# sd(pr)
# ppv=sn_pcr*pr/(sn_pcr*pr+(1-sp_pcr)*(1-pr))
# npv=sp_pcr*(1-pr)/((1-sn_pcr)*pr+sp_pcr*(1-pr))
# plot(ppv, npv, main="Very low 1%+/1% prevalence (e.g. USA); medians indicated by solid lines", xlab="Pr(infected if test positive)", ylab="Pr(not infected if test negative)", col="grey")
# abline(v=median(ppv))
# abline(h=median(npv))
# TP=testpos*ppv
# hist(TP, main=bquote("actual infected among "~ .(testpos)~ "tested positive in USA;"~"median"~.(round(median(TP,2)))), xlab="true positive")
# abline(v=median(TP), col="red")
# #	TP=testposp_pcrA*ppv
# #	hist(TP, main="actual infected among 22 tested positive in PA")
# #	abline(v=median(TP), col="red")
# FN=(testedUSA-testpos)*(1-npv)
# hist(FN, main=bquote("infected among"~.(testedUSA-testpos)~ "who tested negative in USA;"~"median"~.(round(median(FN,2)))), xlab="false negative")
# abline(v=median(FN), col="red")
# 
# 
# ### DISTRIBUTIONS of PREDICTED COUNTS UNDER SEROLOGICAL TESTING ###
# 
# #3x3 plotting window
# par(mfrow=c(3,3))
# 
# #high prevalence: at equillibrium of 50% (China, Italy) for R=2
# #pr=rbeta(n, 9, 1) #90+/-10%%
# pr=rbeta(n, 100, 100) #50+/-3%%
# mean(pr)
# sd(pr)
# ppv=sn_sero*pr/(sn_sero*pr+(1-sp_sero)*(1-pr))
# npv=sp_sero*(1-pr)/((1-sn_sero)*pr+sp_sero*(1-pr))
# plot(ppv, npv, main="Prevalence at equilibirum 50%+/3% for R=2; medians indicated by solid lines", xlab="Pr(infected if test positive)", ylab="Pr(not infected if test negative)", col="grey")
# abline(v=median(ppv))
# abline(h=median(npv))
# TP=testpos*ppv
# hist(TP, main=bquote("actual infected among "~ .(testpos)~ "tested positive in USA;"~"median"~.(round(median(TP,2)))), xlab="true positive")
# abline(v=median(TP), col="red")
# #	TP=testposp_seroA*ppv
# #	hist(TP, main="actual infected vs among tested positive in PA")
# #	abline(v=median(TP), col="red")
# FN=(testedUSA-testpos)*(1-npv)
# hist(FN, main=bquote("infected among"~.(testedUSA-testpos)~ "who tested negative in USA;"~"median"~.(round(median(FN,2)))), xlab="false negative")
# abline(v=median(FN), col="red")
# 
# #low prevalence (USA)
# pr=rbeta(n, 1, 9) #10+/-10%%
# mean(pr)
# sd(pr)
# ppv=sn_sero*pr/(sn_sero*pr+(1-sp_sero)*(1-pr))
# npv=sp_sero*(1-pr)/((1-sn_sero)*pr+sp_sero*(1-pr))
# plot(ppv, npv, main="Low 10%+/10% prevalence (e.g. USA); medians indicated by solid lines", xlab="Pr(infected if test positive)", ylab="Pr(not infected if test negative)", col="grey")
# abline(v=median(ppv))
# abline(h=median(npv))
# TP=testpos*ppv
# hist(TP, main=bquote("actual infected among "~ .(testpos)~ "tested positive in USA;"~"median"~.(round(median(TP,2)))), xlab="true positive")
# abline(v=median(TP), col="red")
# #	TP=testposp_seroA*ppv
# #	hist(TP, main="actual infected among 22 tested positive in PA")
# #	abline(v=median(TP), col="red")
# FN=(testedUSA-testpos)*(1-npv)
# hist(FN, main=bquote("infected among"~.(testedUSA-testpos)~ "who tested negative in USA;"~"median"~.(round(median(FN,2)))), xlab="false negative")
# abline(v=median(FN), col="red")
# 
# #very low prevalence (USA)
# pr=rbeta(n, 1, 99) #1+/-1%%
# mean(pr)
# sd(pr)
# ppv=sn_sero*pr/(sn_sero*pr+(1-sp_sero)*(1-pr))
# npv=sp_sero*(1-pr)/((1-sn_sero)*pr+sp_sero*(1-pr))
# plot(ppv, npv, main="Very low 1%+/1% prevalence (e.g. USA); medians indicated by solid lines", xlab="Pr(infected if test positive)", ylab="Pr(not infected if test negative)", col="grey")
# abline(v=median(ppv))
# abline(h=median(npv))
# TP=testpos*ppv
# hist(TP, main=bquote("actual infected among "~ .(testpos)~ "tested positive in USA;"~"median"~.(round(median(TP,2)))), xlab="true positive")
# abline(v=median(TP), col="red")
# #	TP=testposp_seroA*ppv
# #	hist(TP, main="actual infected among 22 tested positive in PA")
# #	abline(v=median(TP), col="red")
# FN=(testedUSA-testpos)*(1-npv)
# hist(FN, main=bquote("infected among"~.(testedUSA-testpos)~ "who tested negative in USA;"~"median"~.(round(median(FN,2)))), xlab="false negative")
# abline(v=median(FN), col="red")

