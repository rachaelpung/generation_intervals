load('data/el.10m.RData')


# plot el
plot.el = el[ct.grp == 1 & ct == 1]

plot.el = plot.el %>% 
  group_by(id.x, id.y, day) %>% 
  mutate(diff.interval = ifelse(row_number()==1, NA, interval - lag(interval, default = interval[1])))
plot.el=data.table(plot.el)

index.start = plot.el[, which(diff.interval == 1)]
diff.index.start = c(999, diff(index.start))
index.start = index.start[diff.index.start != 1]-1

index.end = plot.el[, which(diff.interval == 1)]
diff.index.end = c(diff(index.end),999)
index.end = index.end[diff.index.end != 1]

plot.el[index.start, interval.end:=plot.el[index.end, interval]]
plot.el=plot.el[is.na(diff.interval) | (!is.na(diff.interval) & !is.na(interval.end))]
plot.el[,diff.interval:=NULL]
setnames(plot.el, c('interval'), c('interval.start'))
plot.el[is.na(interval.end), interval.end:=interval.start]
plot.el[,interval.end:=interval.end+1]

plot.el[,duration:=interval.end-interval.start]

# el.pairs=unique(plot.el[,.(id.x, id.y)])
el.pairs=plot.el[,.(sum(duration), sum(ct)), by=.(id.x, id.y)]
setorder(el.pairs, V1, V2)


el.pairs[,.N]/2
el.pairs=rbindlist(list(el.pairs[1:41],
                        data.table(id.x = c(9999,9999,9999), 
                                   id.y = c(9999,9999,9999),
                                   V1 = c(9999,9999,9999),
                                   V2 = c(9999,9999,9999)),
                        el.pairs[42:82],
                        data.table(id.x = c(9999,9999,9999), 
                                   id.y = c(9999,9999,9999),
                                   V1 = c(9999,9999,9999),
                                   V2 = c(9999,9999,9999))))


paneller=function(row = 1,column=1, p)
{
  xlm = c(0, 3*24*60/5) # number of 5-min intervals over 3 days
  ylm = c(-1,1)
  
  innermargins = c(1,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # set data
  data=plot.el[id.x == el.pairs[p, id.x] & id.y == el.pairs[p, id.y]]
  
  colTl = '#1B191999'
  colCt = '#00468BFF'
  
  if(data[,.N] != 0 ){
    grid.lines(c(xlm[1],xlm[2]),c(0, 0),default.units = 'native',gp=gpar(col=colTl, lwd=0.25))
    
    for(i in 1:data[,.N]){
      
      grid.polygon(c(data$interval.start[i],data$interval.start[i],data$interval.end[i],data$interval.end[i]),
                   c(-0.075,0.075,0.075,-0.075), default.units = 'native', gp=gpar(col=NA,fill=colCt))
      
    }
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  if(row==43) {
    grid.xaxis(at=c(0,288,576,864),label = c('','','',''), gp=gpar(fontsize=unit(8,'pt')))
    grid.text('Thu',x=unit(2.25,'lines'), y=unit(-1,'lines'))
    grid.text('Fri',x=unit(6.7,'lines'), y=unit(-1,'lines'))
    grid.text('Sat',x=unit(11.2,'lines'), y=unit(-1,'lines'))
    grid.text('Time',y=unit(-2,'lines'))
  }
  
  popViewport()
  popViewport()
  
}


png('figure/hh_edgelist_supp.png',height=16,width=16,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=46, ncol=2)))


# paneller(1,1,1)

for(p in 1:el.pairs[,.N]){ #el.pairs[,.N]
  
  col=ifelse(p<=44, 1, 2)
  row=p-44*(col-1)
  paneller(row,col,p)
  
  
}


popViewport()
popViewport()
dev.off()

rm(paneller)

