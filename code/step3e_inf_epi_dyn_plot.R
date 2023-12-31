source('code/step0_load_library.R')

lightup = function(c, alpha)
{
  z=col2rgb(c)/255
  return(rgb(z[1],z[2],z[3],alpha))  # 0.125 for col3
}

# distributions unadjusted
paneller=function(row = 1,column=1)
{
  xlm=c(0,25); ylm=c(0,0.2)
  
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # data
  if(column==1){
    # is.na(mean.iso) & is.na(var.iso)
    row.index.ref = param.unadj[inf.func.name=='spline.covid.delta' & mean.incub==4 & growth_rate=='0.2', set.no]
    row.index.alt = param.unadj[inf.func.name=='spline.covid.wild' & mean.incub==5 & growth_rate=='-0.2', set.no]
    row.index = c(row.index.ref, row.index.alt); space = c(0,0.008*2) #reference, alternate
  }
  if(column==2){
    row.index.ref = param.unadj[inf.func.name=='spline.covid.delta' & mean.incub==4 & growth_rate==0, set.no]
    row.index.alt = param.unadj[inf.func.name=='spline.covid.wild' & mean.incub==5 & growth_rate==0, set.no]
    row.index = c(row.index.ref, row.index.alt); space = c(0,0.008*2) #reference, alternate
  }
  if(column==3){
    row.index.ref = param.unadj[inf.func.name=='spline.covid.delta' & mean.incub==4 & growth_rate==0, set.no]
    row.index.alt = param.unadj[inf.func.name=='spline.covid.wild' & mean.incub==5 & growth_rate=='0.2', set.no]
    row.index = c(row.index.ref, row.index.alt); space = c(0,0.008*2) #reference, alternate
  }
  
  
  # set colour and symbol
  if(column==1){ col.lines = c('#a50f15', lightup('#fb6a4a',0.5)) }
  if(column==2){ col.lines = c('#08519c', lightup('#6baed6',0.5)) }
  if(column==3){ col.lines = c('#caa502', lightup('#fddc49',0.5)) }
  
  # gridlines
  for(y.gl in seq(0,0.2,0.008)){
    grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
    
  for(x.gl in seq(0,25,1)){
    grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }

  
  
  # hist plots 
  for(i in 1:2){
    data=hist(gen.dist[row.index[i], ], breaks=c(0:30), plot = F)
    
    for(b in 1:(length(data$breaks)-1)){
      if(data$breaks[b+1]<=25){
        grid.polygon(c(data$breaks[b], data$breaks[b+1], data$breaks[b+1], data$breaks[b]),
                     c(0,0,data$counts[b]/sum(data$counts),data$counts[b]/sum(data$counts)), default.units = 'native', gp=gpar(fill=col.lines[i], col=col.lines[i]))
        
      }
    }
    
    mean=round(mean(gen.dist[row.index[i], ]), 1)
    sd=round(sd(gen.dist[row.index[i], ]), 1)
    if(i==1) {
      grid.text(bquote('Ref GI' == .(mean) ~ ' (' ~ .(sd) ~')'), x=13.5, y=0.188, just='left', default.units = 'native', gp = gpar(fontsize = 8))
      grid.polygon(c(12,13,13,12,12), c(0.184,0.184,0.192,0.192,0.184), default.units = 'native', gp=gpar(fill=col.lines[i], col=col.lines[i]))
    }
    if(i==2) {
      grid.text(bquote('Alt GI' == .(mean) ~ ' (' ~ .(sd) ~')'), x=13.5, y=0.172, just='left', default.units = 'native', gp = gpar(fontsize = 8))
      grid.polygon(c(12,13,13,12,12), c(0.168,0.168,0.176,0.176,0.168), default.units = 'native', gp=gpar(fill=col.lines[i], col=col.lines[i]))
    }
    
    
  }
  
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  grid.xaxis(at=seq(0,25,5),label=seq(0,25,5),gp=gpar(fontsize=unit(7,'pt')))
  grid.yaxis(at=seq(0,0.2,0.05),label=seq(0,20,5),gp=gpar(fontsize=unit(7,'pt')))


  
  # labels
  if(row == 1 & column == 1) grid.text('A',x=unit(-2.5,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-2.5,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 3) grid.text('C',x=unit(-2.5,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  
  grid.text('Probability (%)',x=unit(-2.5,'lines'),rot=90)
  grid.text('GI (day)',y=unit(-2,'lines'))
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
    
  popViewport()
  popViewport()
  
}

png('figure/epi_dynamics_unadj_supp.png',height=8,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3)))

paneller(1,1)
paneller(1,2)
paneller(1,3)


popViewport()
popViewport()
dev.off()




# distributions adjusted
paneller=function(row = 1,column=1)
{
  xlm=c(0,25); ylm=c(0,0.35)
  
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # generate data
  if(column==1){ growth_rate = log(2)/3.5 }
  if(column==2){ growth_rate = 0 }
  if(column==3){ growth_rate = -log(2)/3.5 }
  
  set.seed(555)
  t=seq(0,30,1) # 0.000001
  inc.f = dlnorm(t, meanlog = 1.51827713, sdlog = 0.4269913) # incub mean=5, var=5
  inc.f = inc.f/sum(inc.f)
  inc.b = exp(-growth_rate*t)*inc.f/sum(exp(-growth_rate*t)*inc.f)
  inc.b = sample(t, 25, replace=TRUE, prob=inc.b)
  inc.f = sample(t, 1000, replace=TRUE, prob=inc.f)
  # inc = sample(t, nsim, replace=TRUE, prob=inc.b)
  
  # sensitivity analysis
  # adjust incubation period
  inc.b.table = data.table(table(inc.b))
  setnames(inc.b.table, c('t','N'))
  inc.b.table[, P:=N/sum(N)]
  inc.b.table[, t:=as.numeric(t)]
  inc.b.f = exp(growth_rate*inc.b.table$t)*inc.b.table$P/sum(exp(growth_rate*inc.b.table$t)*inc.b.table$P)
  inc = sample(inc.b.table$t, 25, replace=T, prob = inc.b.f)
  
  
  
  # set colour and symbol
  if(column==1){ col.lines = c('#a50f15', lightup('#fb6a4a',0.5)) } #lightup('gray70',0.5)
  if(column==2){ col.lines = c('#08519c', lightup('#6baed6',0.5)) }
  if(column==3){ col.lines = c('#caa502', lightup('#fddc49',0.5)) }
  
  # gridlines
  for(y.gl in seq(0,0.35,0.014)){
    grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
  
  for(x.gl in seq(0,25,1)){
    grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
  
  
  
  # hist plots 
  for(i in c(1,2)){
    if(i==1) {data = hist(inc.f, breaks=c(0:30), plot = F); data_inc = inc.f}
    if(i==3) {data = hist(inc.b, breaks=c(0:30), plot = F); data_inc = inc.b}
    if(i==2) {data = hist(inc, breaks=c(0:30), plot = F); data_inc = inc}
    
    # plot histogram
    for(b in 1:(length(data$breaks)-1)){
      if(data$breaks[b+1]<=25){
        grid.polygon(c(data$breaks[b], data$breaks[b+1], data$breaks[b+1], data$breaks[b]),
                     c(0,0,data$counts[b]/sum(data$counts),data$counts[b]/sum(data$counts)), default.units = 'native', gp=gpar(fill=col.lines[i], col=col.lines[i]))
        
      }
    }
    
    mean=round(mean(data_inc), 1)
    sd=round(sd(data_inc), 1)
    
    if(i==1) {
      grid.text(bquote('True GI' == .(mean) ~ ' (' ~ .(sd) ~')'), x=13.5, y=0.329, just='left', default.units = 'native', gp = gpar(fontsize = 8))
      grid.polygon(c(12,13,13,12,12), c(0.322,0.322,0.336,0.336,0.322), default.units = 'native', gp=gpar(fill=col.lines[i], col=col.lines[i]))
    }
    if(i==2) {
      grid.text(bquote('Adj GI' == .(mean) ~ ' (' ~ .(sd) ~')'), x=13.5, y=0.301, just='left', default.units = 'native', gp = gpar(fontsize = 8))
      grid.polygon(c(12,13,13,12,12), c(0.294,0.294,0.308,0.308,0.294), default.units = 'native', gp=gpar(fill=col.lines[i], col=col.lines[i]))
    }
    
    
    
  }
  
  
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  grid.xaxis(at=seq(0,25,5),label=seq(0,25,5),gp=gpar(fontsize=unit(7,'pt')))
  grid.yaxis(at=seq(0,0.35,0.05),label=seq(0,35,5),gp=gpar(fontsize=unit(7,'pt')))
  
  
  
  # labels
  if(row == 1 & column == 1) grid.text('A',x=unit(-2.5,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-2.5,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 3) grid.text('C',x=unit(-2.5,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  
  grid.text('Probability (%)',x=unit(-2.5,'lines'),rot=90)
  grid.text('GI (day)',y=unit(-2,'lines'))
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}

png('figure/epi_dynamics_adj_supp.png',height=8,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3)))

paneller(1,1)
paneller(1,2)
paneller(1,3)

popViewport()
popViewport()
dev.off()


paneller=function(row = 1,column=1)
{
  
  if(column==1){xlm=c(-0.3,0.3); ylm=c(-3,3)}
  if(column==2){ xlm=c(-0.2,2); ylm=c(0,1) }
  if(column==3){ xlm=c(0,2); ylm=c(0,1) }
 
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # data
  if(column==1) {
    data = row.pair.unadj[fig=='a']
    space = c(0,0.3,0.6)
  }
  
  if(column==2) {
    row.pair.unadj.pow[,`:=`(diff=diff.mean.gen,
                             pow=pow.gen)]
    data = row.pair.unadj.pow
  }
  
  if(column==3) {
    row.pair.adj.pow[,`:=`(diff=diff.mean.gen,
                           pow=pow.gen)]
    data = row.pair.adj.pow
  }
  
  
  # set colour and symbol
  if(column==1){
    col.lines = c('#fb6a4a', '#6baed6', '#fddc49')
  }

  if(column %in% c(2,3)){
    pch.sample.size = c(16,17)
    col.lines.org = c('#fb6a4a', '#6baed6', '#fddc49')
    col.pts.hue = list(c('#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15', '#67000d'),
                       c('#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b'),
                       c('#feec9a', '#fee472', '#fddc49', '#fdd421', '#f2c602', '#caa502', '#a18402'))
    
    col.lines=sapply(1:length(col.lines.org), function(x) lightup(col.lines.org[x], alpha=0.25))
    val=seq(0,7,length.out = 7)
    
  }
  
  # gridlines
  if(column ==1){
    for(y.gl in seq(-3,3,0.25)){
      grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
    
    for(x.gl in seq(-0.3,0.3,0.025)){
      grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
  }
  
  if(column %in% c(2,3)){
    for(y.gl in seq(0,1,0.05)){
      grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
    
    if(column==2){
      for(x.gl in seq(-0.2,2,0.1)){
        grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
      }
    }
    
    if(column==3){
      for(x.gl in seq(0,2,0.1)){
        grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
      }
    }
  }
  
  # growth plots 
  if(column==1){
    for(g in 1:3){
      
      fit=predict(smooth.spline(data[grp==g]$growth.ref, data[grp==g]$diff.mean.gen), seq(-0.3,0.3,0.01))
      index=which(fit$y>=-3 & fit$y<=3)
      grid.lines(fit$x[index], fit$y[index], default.units = 'native', gp=gpar(col=col.lines[g]))
      grid.lines(c(-0.29,-0.26), 2.15-space[g], default.units = 'native', gp=gpar(col=col.lines[g]))
      
      if(g==1){ index = which(fit$x==0.2)}
      if(g==2){ index = which(fit$x==0)}
      if(g==3){ index = which(fit$x==0)}
      
      grid.points(fit$x[index],fit$y[index],default.units = 'native', pch=15,gp=gpar(col=col.lines[g], cex=0.8))
      
      
    }
  }
  
  # power
  if(column %in% c(2,3)){
    data[,diff.incub:=abs(mean.incub.alt-mean.incub.ref)]
    size=sort(unique(data$sample.size), decreasing = T)
    space = c(0,0.05,0.10)
    
    for(s in 1:length(size)){
      
      data.spline=data[sample.size==size[s]]
      if(column==3) data.points=data[sample.size==size[s] & diff<=1.99 & diff>=0.05 & pow<=0.99 ]
      if(column==2) data.points=data[(sample.size==size[s] & diff<=1.99 & diff>=-0.2 & pow<=0.99) | (mean.incub.ref==4 & mean.incub.alt==5)]
      
      
      for(g in 1:3){
        
        if(s==1) {
          if(column==2) grid.lines(c(0.03,0.17)-0.2, 0.825-space[g], default.units = 'native', gp=gpar(col=col.lines.org[g]))
          if(column==2) grid.points(0.1-0.2,0.825-space[g], default.units = 'native', pch=15, gp=gpar(cex=0.8, col=col.lines.org[g]))
        }
        # fit spline
        fit=predict(smooth.spline(data.spline[grp==g]$diff, data.spline[grp==g]$pow), seq(0,2,0.01))
      
        grid.lines(fit$x, fit$y, default.units = 'native', gp=gpar(col=col.lines[g]))
        
        # assign colour to points
        pal=gradient_n_pal(colours = col.pts.hue[[g]], values = val)
        data.points[grp==g, diff.incub.hex:=pal(diff.incub)]
        
        for(i in 1:data.points[grp==g,.N]){
          
          grid.points(data.points[grp==g]$diff[i],
                      data.points[grp==g]$pow[i],
                      default.units = 'native', pch=pch.sample.size[s],
                      gp=gpar(col=data.points[grp==g]$diff.incub.hex[i], cex=0.3))
          
          if(data.points[grp==g]$mean.incub.ref[i]==4 & data.points[grp==g]$mean.incub.alt[i]==5){
            grid.points(data.points[grp==g]$diff[i],
                        data.points[grp==g]$pow[i],
                        default.units = 'native', pch=15,
                        gp=gpar(col=data.points[grp==g]$diff.incub.hex[i], cex=0.8))
          }
          
        }
      }
    }
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  if(column==1){
    grid.xaxis(at=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3),label=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3),gp=gpar(fontsize=unit(7,'pt')))
    grid.yaxis(at=seq(-3,3,1),label=seq(-3,3,1),gp=gpar(fontsize=unit(7,'pt')))
  }
  if(column %in% c(2,3)){
    grid.xaxis(at=seq(0,2,0.5),label=seq(0,2,0.5),gp=gpar(fontsize=unit(7,'pt')))
    grid.yaxis(at=seq(0,1,0.2),label=seq(0,100,20),gp=gpar(fontsize=unit(7,'pt')))
    
  }
  
  # labels
  if(row == 1 & column == 1) grid.text('A',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 3) grid.text('C',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  
  
  if(column==1){
    grid.text('GI difference',x=unit(-3,'lines'),rot=90)
    grid.text('Growth rate',y=unit(-2,'lines'))
  }
  if(column %in% c(2,3)){
    grid.text('Difference in GI (day)',y=unit(-2,'lines'))
    grid.text('Power (%)',x=unit(-3,'lines'),rot=90)
  }
  
  if(column==1) {
    grid.text('Incub: Ref=4, alt=5', x=-0.29, y=2.75, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('Epi dynamics in alt:', x=-0.29, y=2.45, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('decline', x=-0.25, y=2.15, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('constant', x=-0.25, y=1.85, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('growth', x=-0.25, y=1.55, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  }
  
  if(column==2) {
    grid.text('Diff incub period', x=0.05-0.2, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('& epidemic dynamics (unadj)', x=0.05-0.2, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  }
  
  if(column==3){
    grid.text('Diff incub period', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('& epidemic dynamics (adj)', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  }
  
  # symbols
  if(column==2) {
    
    grid.text('ref grow, alt dec', x=0.2-0.2, y=0.825, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('ref cons, alt cons', x=0.2-0.2, y=0.775, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('ref cons, alt grow', x=0.2-0.2, y=0.725, default.units = 'native', just='left', gp = gpar(fontsize = 8))

    
    grid.text('sample size', x=0.05-0.2, y=0.625, just='left', default.units = 'native', gp = gpar(fontsize = 8))
    grid.points(0.075-0.2,0.575,default.units = 'native', pch=16,gp=gpar(col='grey60', cex=0.5))
    grid.points(0.075-0.2,0.525,default.units = 'native', pch=17,gp=gpar(col='grey60', cex=0.5))
    
    grid.text('100', x=0.15-0.2, y=0.575, just='left', default.units = 'native', gp = gpar(fontsize = 6))
    grid.text('25', x=0.15-0.2, y=0.525, just='left', default.units = 'native', gp = gpar(fontsize = 6))
    
  }
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}

png('figure/epi_dynamics.png',height=8,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3)))

paneller(1,1)
paneller(1,2)
paneller(1,3)


popViewport()
popViewport()
dev.off()




# 
# paneller=function(row = 1,column=1)
# {
#   
#   if(column %in% c(1,2)){xlm=c(-0.3,0.3); ylm=c(-3,3)}
#   if(column==3){ xlm=c(0,2); ylm=c(0,1) }
#   
#   innermargins = c(2,2,2,2)
#   
#   pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
#   pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
#   
#   # data
#   if(column==1) {
#     data = row.pair.unadj[fig=='a']
#     space = c(0,0.3,0.6)
#   }
#   
#   if(column==2) {
#     data = row.pair.adj[fig=='a']
#     space = c(0,0.3,0.6)
#   }
#   
#   if(column==3) {
#     row.pair[,`:=`(diff=diff.mean.gen,
#                    pow=pow.gen)]
#     data = row.pair
#   }
#   
#   
#   # set colour and symbol
#   if(column==1){
#     col.lines = c('#fb6a4a', '#6baed6', '#fddc49')
#   }
#   if(column==2){
#     col.lines = c('#fb6a4a', '#6baed6', '#fddc49')
#     col.ci = c(lightup('#fb6a4a',0.2), lightup('#6baed6',0.2), lightup('#fddc49',0.2))
#   }
#   if(column==3){
#     pch.sample.size = c(16,17)
#     col.lines.org = c('#fb6a4a', '#6baed6', '#fddc49')
#     col.pts.hue = list(c('#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15', '#67000d'),
#                        c('#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b'),
#                        c('#feec9a', '#fee472', '#fddc49', '#fdd421', '#f2c602', '#caa502', '#a18402'))
#     
#     col.lines=sapply(1:length(col.lines.org), function(x) lightup(col.lines.org[x], alpha=0.25))
#     val=seq(0,7,length.out = 7)
#     
#   }
#   
#   # gridlines
#   if(column %in% c(1,2)){
#     for(y.gl in seq(-3,3,0.25)){
#       grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
#     }
#     
#     for(x.gl in seq(-0.3,0.3,0.025)){
#       grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
#     }
#   }
#   
#   if(column==3){
#     for(y.gl in seq(0,1,0.05)){
#       grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
#     }
#     
#     for(x.gl in seq(0,2,0.1)){
#       grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
#     }
#   }
#   
#   # growth plots 
#   if(column==1){
#     for(g in 1:3){
#       
#       fit=predict(smooth.spline(data[grp==g]$growth.ref, data[grp==g]$diff.mean.gen), seq(-0.3,0.3,0.01))
#       index=which(fit$y>=-3 & fit$y<=3)
#       grid.lines(fit$x[index], fit$y[index], default.units = 'native', gp=gpar(col=col.lines[g]))
#       grid.lines(c(-0.29,-0.26), 2.15-space[g], default.units = 'native', gp=gpar(col=col.lines[g]))
#       
#       if(g==1){ index = which(fit$x==0.2)}
#       if(g==2){ index = which(fit$x==0)}
#       if(g==3){ index = which(fit$x==0)}
#       
#       grid.points(fit$x[index],fit$y[index],default.units = 'native', pch=15,gp=gpar(col=col.lines[g], cex=0.8))
#       
#       
#     }
#   }
#   
#   # difference
#   if(column==2){
#     
#     for(g in 3:1){
#       
#       # fit=predict(smooth.spline(data[grp==g]$growth.ref, data[grp==g]$diff.mean.gen), seq(-0.3,0.3,0.01))
#       # index=which(fit$y>=-3 & fit$y<=3)
#       # grid.lines(fit$x[index], fit$y[index], default.units = 'native', gp=gpar(col=col.lines[g]))
#       # grid.lines(c(-0.29,-0.26), 2.15-space[g], default.units = 'native', gp=gpar(col=col.lines[g]))
#       
#       r = c('-0.3','-0.2','-0.1','0','0.1','0.2','0.3')
#       grid.polygon(c(data[grp==g & growth.ref%in%r]$growth.ref, rev(data[grp==g & growth.ref%in%r]$growth.ref)),
#                    c(data[grp==g & growth.ref%in%r]$lwr.diff.mean.gen, rev(data[grp==g & growth.ref%in%r]$upp.diff.mean.gen)), default.units = 'native', gp=gpar(col=col.ci[g], fill=col.ci[g]))
#       grid.lines(data[grp==g & growth.ref%in%r]$growth.ref, data[grp==g & growth.ref%in%r]$med.diff.mean.gen, default.units = 'native', gp=gpar(col=col.lines[g]))
#       
#       grid.lines(c(-0.29,-0.26), 2.15-space[g], default.units = 'native', gp=gpar(col=col.lines[g]))
#       
#     }
#     grid.points(row.pair.adj[growth.ref==0 & growth.alt==0 & nsim==1000]$growth.ref,
#                 row.pair.adj[growth.ref==0 & growth.alt==0 & nsim==1000]$med.diff.mean.gen,
#                 default.units = 'native', pch=15,
#                 gp=gpar(cex=0.3, col='black'))
#     grid.lines(c(row.pair.adj[growth.ref==0 & growth.alt==0 & nsim==1000]$growth.ref,row.pair.adj[growth.ref==0 & growth.alt==0 & nsim==1000]$growth.ref),
#                c(row.pair.adj[growth.ref==0 & growth.alt==0 & nsim==1000]$lwr.diff.mean.gen,row.pair.adj[growth.ref==0 & growth.alt==0 & nsim==1000]$upp.diff.mean.gen),
#                default.units = 'native', gp=gpar(col='black'))
#     
#   }
#   
#   # power
#   if(column==3){
#     data[,diff.incub:=abs(mean.incub.alt-mean.incub.ref)]
#     size=sort(unique(data$sample.size), decreasing = T)
#     space = c(0,0.05,0.10)
#     
#     for(s in 1:length(size)){
#       
#       data.spline=data[sample.size==size[s]]
#       data.points=data[sample.size==size[s] & diff<=1.99 & diff>=0.05 & pow<=0.99]
#       
#       for(g in 1:3){
#         
#         if(s==1) {
#           grid.lines(c(0.03,0.17), 0.825-space[g], default.units = 'native', gp=gpar(col=col.lines.org[g]))
#           grid.points(0.1,0.825-space[g], default.units = 'native', pch=15, gp=gpar(cex=0.8, col=col.lines.org[g]))
#         }
#         # fit spline
#         fit=predict(smooth.spline(data.spline[grp==g]$diff, data.spline[grp==g]$pow), seq(0,2,0.01))
#         
#         grid.lines(fit$x, fit$y, default.units = 'native', gp=gpar(col=col.lines[g]))
#         
#         # assign colour to points
#         pal=gradient_n_pal(colours = col.pts.hue[[g]], values = val)
#         data.points[grp==g, diff.incub.hex:=pal(diff.incub)]
#         
#         for(i in 1:data.points[grp==g,.N]){
#           
#           grid.points(data.points[grp==g]$diff[i],
#                       data.points[grp==g]$pow[i],
#                       default.units = 'native', pch=pch.sample.size[s],
#                       gp=gpar(col=data.points[grp==g]$diff.incub.hex[i], cex=0.3))
#           
#           if(data.points[grp==g]$mean.incub.ref[i]==4 & data.points[grp==g]$mean.incub.alt[i]==5){
#             grid.points(data.points[grp==g]$diff[i],
#                         data.points[grp==g]$pow[i],
#                         default.units = 'native', pch=15,
#                         gp=gpar(col=data.points[grp==g]$diff.incub.hex[i], cex=0.8))
#           }
#           
#         }
#       }
#     }
#   }
#   
#   
#   popViewport()
#   pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
#   
#   
#   if(column %in% c(1,2)){
#     grid.xaxis(at=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3),label=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3),gp=gpar(fontsize=unit(7,'pt')))
#     grid.yaxis(at=seq(-3,3,1),label=seq(-3,3,1),gp=gpar(fontsize=unit(7,'pt')))
#   }
#   if(column==3){
#     grid.xaxis(at=seq(0,2,0.5),label=seq(0,2,0.5),gp=gpar(fontsize=unit(7,'pt')))
#     grid.yaxis(at=seq(0,1,0.2),label=seq(0,100,20),gp=gpar(fontsize=unit(7,'pt')))
#     
#   }
#   
#   # labels
#   if(row == 1 & column == 1) grid.text('A',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))
#   if(row == 1 & column == 2) grid.text('B',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))
#   if(row == 1 & column == 3) grid.text('C',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))
#   
#   
#   if(column %in% c(1,2)){
#     grid.text('GI difference',x=unit(-3,'lines'),rot=90)
#     grid.text('Growth rate',y=unit(-2,'lines'))
#   }
#   if(column==3){
#     grid.text('Difference in GI (day)',y=unit(-2,'lines'))
#     grid.text('Power (%)',x=unit(-3,'lines'),rot=90)
#   }
#   
#   if(column==1) {
#     grid.text('Incub: Ref=4, alt=5', x=-0.29, y=2.75, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     grid.text('Epi dynamics in alt:', x=-0.29, y=2.45, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     grid.text('decline', x=-0.25, y=2.15, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     grid.text('constant', x=-0.25, y=1.85, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     grid.text('growth', x=-0.25, y=1.55, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#   }
#   
#   if(column==2) {
#     grid.text('Incub: Ref=4, alt=5', x=-0.29, y=2.75, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     grid.text('Sample size', x=-0.29, y=2.45, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     grid.text('1000', x=-0.25, y=2.15, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     grid.text('100', x=-0.25, y=1.85, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     grid.text('25', x=-0.25, y=1.55, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#   }
#   
#   if(column==3){
#     grid.text('Diff incub period', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     grid.text('& epidemic dynamics', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#   }
#   
#   # symbols
#   if(column==3) {
#     
#     grid.text('ref grow, alt dec', x=0.2, y=0.825, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     grid.text('ref cons, alt cons', x=0.2, y=0.775, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     grid.text('ref cons, alt grow', x=0.2, y=0.725, default.units = 'native', just='left', gp = gpar(fontsize = 8))
#     
#     
#     grid.text('sample size', x=0.05, y=0.625, just='left', default.units = 'native', gp = gpar(fontsize = 8))
#     grid.points(0.075,0.575,default.units = 'native', pch=16,gp=gpar(col='grey60', cex=0.5))
#     grid.points(0.075,0.525,default.units = 'native', pch=17,gp=gpar(col='grey60', cex=0.5))
#     
#     grid.text('100', x=0.15, y=0.575, just='left', default.units = 'native', gp = gpar(fontsize = 6))
#     grid.text('25', x=0.15, y=0.525, just='left', default.units = 'native', gp = gpar(fontsize = 6))
#     
#   }
#   
#   grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
#   
#   
#   popViewport()
#   popViewport()
#   
# }
# 
# png('figure/epi_dynamics.png',height=8,width=24,units='cm',res=300,pointsize=10)
# pushViewport(plotViewport(c(2,2,1,1)))
# pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3)))
# 
# paneller(1,1)
# paneller(1,2)
# paneller(1,3)
# 
# 
# popViewport()
# popViewport()
# dev.off()
