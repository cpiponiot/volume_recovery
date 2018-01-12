
smallplots = subset(md_plot, plot.size<5)
smallplots$newplot = paste(smallplots$site, smallplots$treat, sep="_")
smallplots$newplot[smallplots$site=="ita"] = paste(smallplots$newplot[smallplots$site=="ita"], substring(smallplots$idplot[smallplots$site=="ita"], 5,5), sep="")

newplots = smallplots[,.(dVdt=mean(dVdt), dVMdt=mean(dVMdt), plot.size = sum(plot.size), Vmort = mean(Vmort), V=mean(V), dVGdt=mean(dVGdt)),.(site, newplot,year,treat)]
newplots$idplot = newplots$newplot
tmin = newplots[,.(tmin=year[which.min(V[year<min(year)+5])], t0 = min(year)),.(idplot)]
newplots = merge(newplots, tmin, by="idplot")
newplots$tmin[newplots$treat==0] = newplots$t0[newplots$treat==0]
newplots$t = newplots$year-newplots$tmin

## 1 missing inventory for treatments 2,3,4 in tapajos
newplots$plot.size[newplots$site=="tpj"] = 3

## merge big and small plots
cols = c("site","idplot","t","year","treat","plot.size","V","Vmort","dVdt","dVGdt","dVMdt")
md_grouped = rbind.data.frame(subset(md_plot[,cols,with=FALSE], plot.size>=5), newplots[,cols,with=FALSE])

