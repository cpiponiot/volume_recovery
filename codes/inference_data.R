######################################################################################
##############                       ANALYSIS                       ##################
######################################################################################

setwd("C:/Users/camille.piponiot/Google Drive/volume")
library(rstan)
library(data.table)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("C:/Users/camille.piponiot/Google Drive/radam/site_info.Rdata")
# source("codes/write_metadata.R")
load("data/metadata.Rdata")

### change t0 so that for logged plots it corresponds to the minimum volume during 
### the 5 years following logging events; if the minimum volume occurs during a different period, discard the plot
minV = md_plot[order(t),.(t[which.min(V)], first(V) - min(V[t%in%0:5])),.(idplot,site,treat)]
colnames(minV)[4:5] <- c("tmin","deltaV")
minV_keep = subset(minV, treat == 0 | ( tmin <= 5 & tmin >= 0) |  (tmin<0 & treat==3 & site=="prc") )  ## in paracou, t=0 in treatment 3 in 1990, 3 years after logging
plot_pb = minV$idplot[!(minV$idplot %in% minV_keep$idplot)]
# par(mfrow=c(2,3))
# sapply(minV_keep$idplot, function(i) plot(md_plot$V[md_plot$idplot==i]~md_plot$t[md_plot$idplot==i], type="l", main=paste(i, unique(md_plot$treat[md_plot$idplot==i]))))

md_plot = merge(md_plot, site_info, by="site")
md_plot = merge(md_plot, minV, by=c("idplot","site","treat"))
md_plot$t0 = md_plot$tmin
md_plot$t0[md_plot$treat==0] = 0
md_plot$t_before = md_plot$t
md_plot$t = md_plot$t - md_plot$t0
md_plot$ns = as.numeric(as.factor(md_plot$site))
md_plot$np = as.numeric(as.factor(md_plot$idplot))
md_inf = subset(md_plot, t>0)
columns = c("site","prec","seas","rad","bkd","cec","cfr","psd","dep")
site_md = unique(md_inf[,columns,with=FALSE])
site_md = cbind(site_md, scale(site_md[,-1]))
colnames(site_md) = c(columns, paste(columns[-1], "sd", sep="_"))

## cumulative fluxes
md_inf = md_inf[order(t)]
cumvol = md_inf[,.(t,cumsum(dVGdt),cumsum(dVMdt)),.(np)]
colnames(cumvol) = c("np","t","cumVG","cumVM")
md_inf = merge(md_inf, cumvol)
## plot info
per_plot = md_plot[order(np,t),.(s=unique(ns), treat=unique(treat), deltaV=unique(deltaV)),.(np)]
ps = per_plot$s
logged = which(per_plot$treat>0)
md_inf$logged = md_inf$np %in% logged
deltaV = per_plot$deltaV[per_plot$treat>0]

## choice of priors ##
# radambrasil data
load("C:/Users/camille.piponiot/Google Drive/radam/vol_data.Rdata")
top95 = subset(vol_data, vol_tot>=quantile(vol_data$vol_tot, 0.95))$vol_tot
mean(top95); sd(top95)

# save(md_inf, file="data/md_inference.Rdata")

N=nrow(md_inf);  S=max(md_inf$ns); P=max(md_inf$np) 
L = length(logged); ns = md_inf$ns; np = md_inf$np
t = md_inf$t; cVG = md_inf$cumVG; cVM = md_inf$cumVM; V = md_inf$V

## pre-logging volumes: 95 and 99.9 percentiles to bound vmax
# Vpre = md_plot[order(t),first(V),.(idplot,site)]
Vpre = md_plot[,.(V0=V[which.min(t)]),.(site, np)]
Vlog = Vpre[,.(mean(log(V)), sd(log(V)), mean(V), sd(V), quantile(V, 0.95)), .(site)]
colnames(Vlog) = c("site","meanVlog","sdVlog","meanV","sdV","q95")
Vlog = Vlog[order(Vlog$site)]
meanVlog = Vlog$meanVlog; dim(meanVlog) = S
sdVlog =  Vlog$sdVlog; dim(sdVlog) = S
v90 = qlnorm(p=0.90, meanlog = meanVlog, sdlog = sdVlog)
v95 = qlnorm(p=0.95, meanlog = meanVlog, sdlog = sdVlog)
v999 = qlnorm(p=0.999, meanlog = meanVlog, sdlog = sdVlog)

weight = tapply(md_plot$plot.size, md_plot$np, unique)
weight = weight/sum(weight)

meanV = md_plot[order(ns), .(mean(V)),.(ns)]$V1
sdV = md_plot[order(ns), .(sd(V)),.(ns)]$V1

dtV0 = md_plot[, .(V0 = V[which.min(year)]),.(site,idplot)]
meanV0 = dtV0[, .(mean(V0)),.(site)]$V1
sdV0 = dtV0[, .(sd(V0)),.(site)]$V1

boxplot(V0~site, data=dtV0)

plot(meanV0, ylim=c(0, max(meanV0+2*sdV0)), xlim=c(0,15))
segments(x0 = 1:14, y0 = meanV0 + 2*sdV0, y1= meanV0 - 2*sdV0)
text(1:14, meanV0, labels = sort(unique(dtV0$site)), pos = 2)
