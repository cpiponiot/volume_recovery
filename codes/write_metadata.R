
setwd("C:/Users/camille.piponiot/Google Drive/volume")

library(data.table)
library(stringdist)
library(scales)
library(tidyr)

###################################################################################
###################                  OPEN DATA                  ################### 
###################################################################################

# source("C:/Users/camille.piponiot/Google Drive/Data TmFO cleaned/opening_and_correcting_data.R")
load("C:/Users/camille.piponiot/Google Drive/Data TmFO cleaned/new_data/dataDBH.Rdata")
data <- data.table(data)

data = subset(data, dbh_c_50>=50 & !is.na(dbh_c_50))
data$dbh = data$dbh_c_50
data = data[,c("site","plot","idtree","year","dbh"),with=FALSE]
data = data[order(site,plot,idtree,year)]

# repeat last measurement on year of death
source("C:/Users/camille.piponiot/Google Drive/Data TmFO cleaned/functions/add_missing.R")
source("C:/Users/camille.piponiot/Google Drive/Data TmFO cleaned/functions/tree_status.R")
# (1) add all censuses per site/plot
missing = data[,.(add_missing(yr=year,id=idtree)),.(site,plot)]
missing$idtree = tstrsplit(missing$V1,split = "\\s+")[[1]]
missing$year = as.numeric(tstrsplit(missing$V1,split = "\\s+")[[2]])
missing$V1 = NULL
data = merge(data, missing, by=c("site","plot","idtree","year"),all=TRUE)
# (2) calculate status and keep status = 1 (alive) and status = 0 (dead)
status = data[,.(year,tree_status(dbh)),.(idtree)]
colnames(status)[3] = "status"
data = merge(data, status, by=c("idtree","year"))
data = subset(data, !is.na(status))
data = data[order(site,idtree,year)]
# keep only first census dead and replicate last DBH
data$diffstatus = c(0,diff(data$status))
data = subset(data, !(diffstatus==0 & status==0))
data$dbh[data$status==0] = data$dbh[which(data$status==0)-1]
data$diffstatus = NULL

data = subset(data, !(site=="prc" & plot %in% c("13","14","15","16","17","18")))
data$subplot = tstrsplit(data$idtree, split = "_")[[3]]
data = subset(data, subplot!=0)

## get plot size and treatment ##
load("C:/Users/camille.piponiot/Google Drive/Data TmFO cleaned/new_data/plot_data.Rdata")
plot_data = subset(plot_data, (site %in% c("lch","inp") & plot.size > 16)|!(site%in%c("lch","inp")))[,c("site","plot","plot.size","treat")]
plot_data = unique(plot_data)
plot_data$plot.size[plot_data$site=="prc"] = 6.25
data = merge(data, plot_data, by=c("site","plot"))
## add treatment info
### 0 = control, 1 = conventional logging, 2 = RIL, 3 = sylvicultural treatments
data$treat[data$treat == "ctrl"] = 0
data$treat[data$treat %in% c("CL","int")] = 1
data$treat[data$treat == "RIL"] = 2
data$treat[data$treat == "sylv"] = 3

# calculate volumes
# source("codes/vol_equations.R")
load("C:/Users/camille.piponiot/Google Drive/radam/site_info.Rdata")
data = merge(data, site_info[,c("site","a","b","site_name")], by="site")
data$vol = data$a*data$dbh^data$b
data$a = NULL
data$b = NULL

### starting year
data$idplot = paste(data$site, data$plot, sep="_")
ti = sapply(unique(data$idplot), function(x){
  tm = sort(unique(data$year[data$idplot==x]))
  site = unique(data$site[data$idplot==x])
  ttm = as.numeric(unique(data$treat[data$idplot==x]))
  if (site!="prc" & ttm>0) return(tm[2]) else if (ttm==0) return(tm[1]) else if (site=="prc" & ttm==1) return(tm[4])  else return(tm[7])
})
ti = data.table(idplot= names(ti), ti=ti)
data = merge(data, ti, by="idplot")
data$t = data$year - data$ti

## selecting plots in TmFO ##
data = subset(data, !(site=="ita" & plot=="B4") )  ### a trail was built in the middle


###################################################################################
###################                 ANALYZE DATA                ################### 
###################################################################################

# ## for paracou & lachonta: because plots are very large, do it by  subplot
# data$idplot[data$site=="prc"] = paste(data$idplot[data$site=="prc"], data$subplot[data$site=="prc"], sep="_")
# data$plot.size[data$site=="prc"] = 6.25/4
# data$idplot[data$site=='lch'] = paste(data$idplot[data$site=="lch"], data$subplot[data$site=="lch"], sep="_")
 
# ### calculate volume ###
# source("codes/vol_equations.R")
# volumes = data[,.(vol_equations(dbh, site)),.(year,idtree)]
# colnames(volumes)[3] = "vol"
# data = merge(data, volumes, by=c("year","idtree"))

# divide paracou into subplots ("carrés")
# data$idplot[data$site=='prc'] = paste(data$idplot[data$site=="prc"], data$subplot[data$site=="prc"],sep="_")

## rq: V of commercial species, C of all species
md_plot = data[,.(Vmort = sum(vol * (1-status)/plot.size),            # Vmort /ha
                  V = sum(vol * status/plot.size)),                   # V /ha
               .(t,year,site,idplot,treat,plot.size)]
md_plot = md_plot[order(t)]

diffv = md_plot[,.(t,c(0,diff(V)/diff(t)), c(0,Vmort[-1]/diff(t))),.(site,idplot)]
colnames(diffv)=c("site","idplot","t","dVdt","dVMdt")
md_plot = merge(diffv,md_plot,by=c("site","idplot","t"))
md_plot$dVGdt = md_plot$dVdt + md_plot$dVMdt

vol_data = subset(data, status==1)[,.(sum(vol/plot.size)),.(site,idplot,year,treat,plot.size, site_name)]

rm(data)

# ### volume per site/plot
# vol_data = subset(data, status==1)[,.(sum(vol/plot.size)),.(site,idplot,year,treat)]
# par(mfrow=c(2,3))
# sapply(unique(vol_data$site), function(i){ 
#   dt = subset(vol_data, site==i)
#   dt = dt[order(idplot,year)]
#   # mat = spread(data = dt, key = year, value = V1)
#   plot(dt$V1~dt$year, col=as.numeric(dt$treat)+1, pch=16,main=i,xlab="t",ylab="Volume (m3/ha)")
#   sapply(unique(dt$idplot), function(j) lines(dt$V1[dt$idplot==j]~dt$year[dt$idplot==j], col=as.numeric(dt$treat[dt$idplot==j])+1))
# })
#  
# sapply(unique(md_plot$site), function(i){ 
#   dt = subset(md_plot, site==i)
#   dt = dt[order(idplot,t)]
#   plot(dt$dVdt~dt$t, col=as.numeric(dt$treat)+1, pch=16,main=i,xlab="t",ylab="Volume (m3/ha)")
#   sapply(unique(dt$idplot), function(j) lines(dt$dVdt[dt$idplot==j]~dt$t[dt$idplot==j], col=as.numeric(dt$treat[dt$idplot==j])+1))
# })
# 
# ## lachonta
# dt = subset(vol_data, site=="lch")
# dt = dt[order(year),]
# par(mfrow=c(1,1), mar=c(4,4,2,1))
# plot(dt$V1~dt$year, col=as.numeric(dt$treat)+1, pch=16, ylim=c(0,max(dt$V1)))
# sapply(unique(dt$plot), function(i){
#   lines(dt$V1[dt$plot==i]~dt$year[dt$plot==i], col=as.numeric(dt$treat[dt$plot==i])+1)
# })
# 
