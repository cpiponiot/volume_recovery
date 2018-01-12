setwd("C:/Users/camille.piponiot/Google Drive/volume")


source("gfbi/codes/volume_equations.R")

### provide aP[nrun, S], vmax[nrun,s], theta[nrun], bP[nun], bM[nrun] ###
###              and sdV[nrun, S], sdG[nrun,s], sdM[nrun],            ### 
###              and maxl; and logged and ns and np; and 'title'      ###
###               and ti[nrun,S], and tlog[nrun, logged]              ###    

##############################################################################
##########       predictions on volume, growth and mortality,       ########## 
##########       with confidence interval and data, per site        ########## 
##############################################################################


samp = sample(nrow(aP), 1000)
tvec = 1:500

ti_pred_maxl = ti[maxl, ns];  
for (i in 1:length(logged)) { ti_pred_maxl[np==logged[i]] = tlog[maxl,i]}

  
pdf(paste("interm results/pred",title,".pdf"), width=4, height=8)
par(mfrow=c(3,1), oma = c(5,0,5,0), mar=c(0,4,0,1))

sapply(1:12, function(s){
  predV = sapply(samp, function(i){
    rlnorm(length(tvec), mean=log(volume(tvec, ag = aP[i,s], 
                                         am = aP[i,s] - vmax[i,s]*theta[i], 
                                         th = theta[i], bg=bP[i], bm=bM[i])), sd=sdV[i])
  })
  predG = sapply(samp, function(i){
    rlnorm(length(tvec), mean=log(dVG(tvec, ag = aP[i,s], 
                                         am = aP[i,s] - vmax[i,s]*theta[i], 
                                         th = theta[i], bg=bP[i], bm=bM[i])), sd=sdG[i])
  })
  predM = sapply(samp, function(i){
    rlnorm(length(tvec), 
           mean=log(dVM(tvec, bm=bM[i], am=aP[i,s]-vmax[i,s]*theta[i])), 
           sd=sdM[i])
  })
  
  ICV = apply(predV, 1, function(X) quantile(X,probs = c(0.025,0.975)))
  ICG = apply(predG, 1, function(X) quantile(X,probs = c(0.025,0.975)))
  ICM = apply(predM, 1, function(X) quantile(X,probs = c(0.025,0.975)))
  
  # volume plot
  plot(y = V[ns==s], x = t[ns==s] + ti_pred_maxl[ns==s], 
       ylab="V (m3/ha)", xaxt="n",
       xlim=c(0,500),ylim=c(0,max(ICV)),
       pch=16, col=md_inf$np[md_inf$ns==s])
  polygon(x = c(tvec,rev(tvec)), y = c(ICV[1,], rev(ICV[2,])), col = "#00000020", border=NA)
  mtext(outer = TRUE, unique(md_plot$site_name)[s], side = 3, line=2)
  # growth
  plot(y = md_inf$dVGdt[ns==s], x = t[ns==s] + ti_pred_maxl[ns==s], 
       ylab="Growth (m3/ha/yr)", xaxt="n",
       xlim=c(0,500),ylim=c(0,max(max(ICG), md_inf$dVGdt[ns==s])),
       pch=16, col=md_inf$np[md_inf$ns==s])
  polygon(x = c(tvec,rev(tvec)), y = c(ICG[1,], rev(ICG[2,])), col = "#00000020", border=NA)
  # mortaltiy
  plot(y = md_inf$dVMdt[ns==s], x = t[ns==s] + ti_pred_maxl[ns==s], 
       ylab="Mortality (m3/ha/yr)", xlab="",
       xlim=c(0,500),ylim=c(0,max(max(ICM), md_inf$dVMdt[ns==s])),
       pch=16, col=md_inf$np[md_inf$ns==s])
  polygon(x = c(tvec,rev(tvec)), y = c(ICM[1,], rev(ICM[2,])), col = "#00000020", border=NA)
  mtext(outer=TRUE, "Maturity (yr)", side=1, line=3)
})
dev.off()


##############################################################################
##########     Goodness of fit, with 95 credibility intervals       ########## 
##############################################################################

samp = c(maxl, samp)

ti_pred = ti[, ns];  
for (i in 1:length(logged)) { ti_pred[,np==logged[i]] = tlog[,i]}
# see differences in ti_pred among sites and treatments
# res = sapply(1:S, function(site){
#   sapply(unique(md_inf$treat[md_inf$ns==site]), function(tr){
#     ti_trs = ti_pred[,md_inf$ns==site&md_inf$treat==tr]
#     return(quantile(c(ti_trs), probs=c(0.025,0.5,0.975)))
#   })
# })
# res1 = data.table(t(do.call(cbind,res)))
# sites_names = unique(md_inf$site[order(ns)])
# res1$site = rep(sites_names, unlist(lapply(res,ncol)))
# res1$treat = unlist(lapply(res,colnames))
# plot(res1$`50%`, col=as.numeric(res1$treat)+1, pch=16, ylab="ti")
# segments(y0=res1$`2.5%`, y1=res1$`97.5%`,x0=1:nrow(res1), col=as.numeric(res1$treat)+1)

pred = sapply(1:N, function(n){
  ind_pred = sapply(samp, function(i){
    predV = rlnorm(1, mean=log(volume(ti_pred[i,n], ag = aP[i,ns[n]], 
                                      am = aP[i,ns[n]] - vmax[i,ns[n]]*theta[i], 
                                      th = theta[i], bg=bP[i], bm=bM[i])), sd=sdV[i])
    predG = rlnorm(1, mean=log(cum_dVG(md_inf$t[n]+ti_pred[i,np[n]], ag = aP[i,ns[n]], t0 = ti_pred[i,np[n]],
                                       am = aP[i,ns[n]] - vmax[i,ns[n]]*theta[i], 
                                       th = theta[i], bg=bP[i], bm=bM[i])), sd=sdG[i])
    predM =  rlnorm(1, mean=log(cum_dVM(md_inf$t[n]+ti_pred[i,np[n]], bm=bM[i],  t0 = ti_pred[i,np[n]],
                                        am = aP[i,ns[n]] - vmax[i,ns[n]]*theta[i])), sd=sdM[i])
    return(c(predV,predG,predM))
  })
  CI = rbind(ind_pred[,1], apply(ind_pred,1,function(x) quantile(x, probs=c(0.025,0.975))))
  return(CI)
})

pred = data.table(t(pred))
colnames(pred) = c("V_maxl","V2_5","V97_5","G_maxl","G2_5","G97_5","M_maxl","M2_5","M97_5")

pdf(paste("interm results/gof_",title,".pdf", sep=""), height=4, width=9)
par(mfrow=c(1,3), mar = c(4,2,4,1), oma = c(3,3,3,0))
sapply(1:S, function(s){
  plot(pred$V_maxl[ns==s]~md_inf$V[ns==s], pch=16, col=np, cex=1.5, xlab = "Volume",ylab="", xlim=c(0,350), ylim=c(0,350)); abline(0,1)
  segments(y0=pred$V2_5[ns==s], y1=pred$V97_5[ns==s], x0=md_inf$V[ns==s], col=np)
  plot(pred$G_maxl[ns==s]~md_inf$cumVG[ns==s], pch=16, col=np, cex=1.5, xlab="Cumulative growth", ylab="",xlim=c(0,35), ylim=c(0,35)); abline(0,1)
  segments(y0=pred$G2_5[ns==s], y1=pred$G97_5[ns==s], x0=md_inf$cumVG[ns==s], col=np)
  plot(pred$M_maxl[ns==s]~md_inf$cumVM[ns==s], pch=16, col=np, cex=1.5, xlab="Cumulative mortality", ylab="",xlim=c(0,35), ylim=c(0,35)); abline(0,1)
  segments(y0=pred$M2_5[ns==s], y1=pred$M97_5[ns==s], x0=md_inf$cumVM[ns==s], col=np)
  mtext(outer = T, unique(md_inf$site_name[ns==s]), side=3)
  mtext(outer = TRUE, "Observed", side=1, line=1)
  mtext(outer = TRUE, "Predicted", side=2, line=1)
})
dev.off()



