source('code/step0_load_library.R')
# source('code/step0_load_colour.R')

lightup = function(c, alpha)
{
  z=col2rgb(c)/255
  return(rgb(z[1],z[2],z[3],alpha))  # 0.125 for col3
}


paneller=function(row = 1,column=1)
{
  xlm=c(1,10); ylm=c(-0.25,1.25) 
  
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # set colour and symbol
  pch.shape = c(15,16,17)
  col.lines.gen = c('#fb6a4a', '#6baed6', '#fddc49')
  col.lines.si = c('#FEAC72', '#8c96c6', '#78c679')
  
  # different points for different mean incubation period, different colours for different peak infectiousness
  
  
  # gridlines
  for(y.gl in seq(-0.25,1.25,1/6)){
    grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
  
  for(x.gl in seq(1,10,1)){
    grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
  
  
  inc=c(2,4,6)
  peak.beta=c(0.0001,0.0005,0.001)
  
  # for(i in 1:length(inc)){ 
  i=column
  data = row.cluster[mean.incub==inc[i] & beta==0.0005 & sample.size==100]
  setorder(data, mean.iso)
  
  col.lines=col.lines.gen
  fit.mean=predict(smooth.spline(data$mean.iso, data$diff.mean.gen, nknots=4), seq(1,10,0.01))
  fit.low=predict(smooth.spline(data$mean.iso, data$diff.mean.gen.ci.low, nknots=4), seq(1,10,0.01))
  fit.upp=predict(smooth.spline(data$mean.iso, data$diff.mean.gen.ci.upp, nknots=4), seq(1,10,0.01))
  
  
  grid.polygon(c(fit.low$x,rev(fit.upp$x)), c(fit.low$y,rev(fit.upp$y)), default.units = 'native', gp=gpar(fill=lightup(col.lines[i],0.4), col=lightup(col.lines[i],0.4)))
  
  grid.lines(fit.mean$x, fit.mean$y, default.units = 'native', gp=gpar(col=col.lines[i]))
  
  grid.points(data$mean.iso, data$diff.mean.gen, default.units = 'native', pch=pch.shape[i],
              gp=gpar(col=col.lines[i], cex=0.3))
  
  col.lines=col.lines.si
  fit.mean=predict(smooth.spline(data$mean.iso, data$diff.mean.si, nknots=4), seq(1,10,0.01))
  fit.low=predict(smooth.spline(data$mean.iso, data$diff.mean.si.ci.low, nknots=4), seq(1,10,0.01))
  fit.upp=predict(smooth.spline(data$mean.iso, data$diff.mean.si.ci.upp, nknots=4), seq(1,10,0.01))
  
  if(column!=1) grid.polygon(c(fit.low$x,rev(fit.upp$x)), c(fit.low$y,rev(fit.upp$y)), default.units = 'native', gp=gpar(fill=lightup(col.lines[i],0.4), col=lightup(col.lines[i],0.4)))
  if(column==1) grid.polygon(c(fit.low$x,rev(fit.upp$x)), c(fit.low$y,rev(fit.upp$y)), default.units = 'native', gp=gpar(fill=lightup(col.lines[i],0.4), col=lightup(col.lines[i],0.4)))
  
  grid.lines(fit.mean$x, fit.mean$y, default.units = 'native', gp=gpar(col=col.lines[i]))
  
  grid.points(data$mean.iso, data$diff.mean.si, default.units = 'native', pch=pch.shape[i],
              gp=gpar(col=col.lines[i], cex=0.3))
  
  
  
  
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  grid.xaxis(at=seq(1,10,1),label=seq(1,10,1),gp=gpar(fontsize=unit(7,'pt')))
  grid.yaxis(at=seq(-0.25,1.25,0.25),label=seq(-0.25,1.25,0.25),gp=gpar(fontsize=unit(7,'pt')))
  
  # labels
  if(row == 1 & column == 1) grid.text('A',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 3) grid.text('C',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  
  grid.text('Onset to isolation (days)',y=unit(-2,'lines'))
  grid.text('Difference (days)',x=unit(-3,'lines'),rot=90)
  
  # legend
  
  grid.polygon(c(7/9,7/9+1/36,7/9+1/36,7/9), c(2/9-1/36,2/9-1/36,2/9,2/9), default.units = 'npc', gp=gpar(fill=lightup(col.lines.gen[i],0.4), col=lightup(col.lines.gen[i],0.4)))
  grid.lines(c(7/9,7/9+1/36), c(2/9-1/72,2/9-1/72), default.units = 'npc', gp=gpar(col=col.lines.gen[i]))
  grid.text(bquote(~omega),default.units = 'npc',y=2/9-1/72,x=7/9+2/36, gp=gpar(fontsize=unit(7,'pt')))
  
  if(column!=1) grid.polygon(c(7/9,7/9+1/36,7/9+1/36,7/9), c(2/9-2/36-2/72,2/9-2/36-2/72,2/9-1/36-2/72,2/9-1/36-2/72), default.units = 'npc', gp=gpar(fill=lightup(col.lines.si[i],0.4), col=lightup(col.lines.si[i],0.4)))
  if(column==1) grid.polygon(c(7/9,7/9+1/36,7/9+1/36,7/9), c(2/9-2/36-2/72,2/9-2/36-2/72,2/9-1/36-2/72,2/9-1/36-2/72), default.units = 'npc', gp=gpar(fill=lightup(col.lines.si[i],0.4), col=lightup(col.lines.si[i],0.4)))
  grid.lines(c(7/9,7/9+1/36), c(2/9-2/36-1/72,2/9-2/36-1/72), default.units = 'npc', gp=gpar(col=col.lines.si[i]))
  grid.text(bquote(~sigma),default.units = 'npc',y=2/9-2/36-1/72,x=7/9+2/36, gp=gpar(fontsize=unit(7,'pt')))
  
  grid.text(bquote('Incubation of'~.(inc[i])~'days'),default.units = 'npc',x=0.025,y=0.95, gp=gpar(fontsize=unit(8,'pt')), just='left')
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/hh_gen_si.png',height=8,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3)))

paneller(1,1)
paneller(1,2)
paneller(1,3)

popViewport()
popViewport()
dev.off()

