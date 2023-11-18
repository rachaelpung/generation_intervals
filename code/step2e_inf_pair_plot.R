source('code/step0_load_library.R')

lightup = function(c, alpha)
{
  z=col2rgb(c)/255
  return(rgb(z[1],z[2],z[3],alpha))  # 0.125 for col3
}

paneller=function(row = 1,column=1, type=NULL)
{
  xlm=c(0,2); ylm=c(0,1) 
  
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # data
  if(type=='gen') row.pair[,`:=`(diff=diff.mean.gen,
                                 pow=pow.gen)]
  if(type=='si') row.pair[,`:=`(diff=diff.mean.si,
                                pow=pow.si)]
  if(type=='dgi') row.pair[,`:=`(diff=diff.mean.si,
                                pow=pow.dgi)]
  
  if(row == 1) data = row.pair[fig=='a']
  if(row == 2) data = row.pair[fig=='b']
  if(row == 3) data = row.pair[fig=='c']
  

  # set colour and symbol
  pch.sample.size = c(16,17)
  if(row!=2) col.pts.hue = list(c('#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15', '#67000d'),
                                c('#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b'),
                                c('#feec9a', '#fee472', '#fddc49', '#fdd421', '#f2c602', '#caa502', '#a18402'))
  
  if(row==2) col.pts.hue = list(c('#fdd0a2', '#fdae6b', '#fd8d3c', '#f16913', '#d94801', '#a63603', '#7f2704'),
                                c('#bfd3e6', '#9ebcda', '#8c96c6', '#8c6bb1', '#88419d', '#810f7c', '#4d004b'),
                                c('#d9f0a3', '#addd8e', '#78c679', '#41ab5d', '#238443', '#006837', '#004529'))
  
  if(row!=2) col.lines = c('#fb6a4a', '#6baed6', '#fddc49')
  if(row==2) col.lines = c('#fd8d3c', '#8c96c6', '#78c679')
  
  col.lines=sapply(1:length(col.lines), function(x) lightup(col.lines[x], alpha=0.25))
  
  if(row %in% c(1,3)) val=seq(0,7,length.out = 7)
  if(row == 2) val=seq(0,3.5,length.out = 7)
 
  # gridlines
  for(y.gl in seq(0,1,0.05)){
    grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
    
  for(x.gl in seq(0,2,0.1)){
    grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
  
  
  # 1: power to detect difference in GI due to changes in incubation period
  # 2: power to detect difference in GI due to changes in incubation, infectiousness period and peak infectiousness
  # 3: power to detect difference in GI due to changes in incubation, infectiousness period
  data[,diff.incub:=abs(mean.incub.alt-mean.incub.ref)]
  size=sort(unique(data$sample.size), decreasing = T)
  
  for(s in 1:length(size)){
    
    data.spline=data[sample.size==size[s]]
    if(column != 3) data.points=data[sample.size==size[s] & diff<=1.99 & diff>=0.05 & pow<=0.99]
    if(column == 3) data.points=data[sample.size==size[s] & diff<=1.99 & diff>=0.05] # for aesthetics
    
    for(g in 1:3){
      
      # fit spline
      if(row==1 & column==3 & g==3 & s==1) {
        fit=predict(smooth.spline(data.spline[grp==g]$diff, data.spline[grp==g]$pow), seq(0,0.5,0.01)) # for aesthetics
      } else if(row==3 & column==3 & g==3 & s==1) {
        fit=predict(smooth.spline(data.spline[grp==g]$diff, data.spline[grp==g]$pow), seq(0,0.5,0.01)) # for aesthetics
      } else{
        fit=predict(smooth.spline(data.spline[grp==g]$diff, data.spline[grp==g]$pow), seq(0,2,0.01))
      }
      
      if(!(row==2 & column==3 & g==3))grid.lines(fit$x, fit$y, default.units = 'native', gp=gpar(col=col.lines[g]))
      
      # assign colour to points
      pal=gradient_n_pal(colours = col.pts.hue[[g]], values = val)
      data.points[grp==g, diff.incub.hex:=pal(diff.incub)]
      
      for(i in 1:data.points[grp==g,.N]){
        
        grid.points(data.points[grp==g]$diff[i],
                    data.points[grp==g]$pow[i],
                    default.units = 'native', pch=pch.sample.size[s],
                    gp=gpar(col=data.points[grp==g]$diff.incub.hex[i], cex=0.3))
    
      }
    }
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  grid.xaxis(at=seq(0,2,0.5),label=seq(0,2,0.5),gp=gpar(fontsize=unit(7,'pt')))
  grid.yaxis(at=seq(0,1,0.2),label=seq(0,100,20),gp=gpar(fontsize=unit(7,'pt')))
  
  # fig labels
  if(row == 1 & column == 1) grid.text('A',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 3) grid.text('C',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 1) grid.text('D',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 2) grid.text('E',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 3) grid.text('F',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 3 & column == 1) grid.text('G',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 3 & column == 2) grid.text('H',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 3 & column == 3) grid.text('I',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  
  # top label
  if(row ==1 & column == 1) grid.text('Theoretical GI',y=unit(14.5,'lines'),gp = gpar(fontsize = 10,fontface='bold'))
  if(row ==1 & column == 2) grid.text('Observed SI',y=unit(15.3,'lines'),gp = gpar(fontsize = 10,fontface='bold'))
  if(row ==1 & column == 2) grid.text('(Lower limit for inferred GI)',y=unit(14.5,'lines'),gp = gpar(fontsize = 10,fontface='bold'))
  if(row ==1 & column == 3) grid.text('(Upper limit for inferred GI)',y=unit(14.5,'lines'),gp = gpar(fontsize = 10,fontface='bold'))
  
  # if(column == 1) grid.text(bquote('Difference in '~omega~' (day)'),y=unit(-2,'lines'))
  # if(column == 2) grid.text(bquote('Difference in '~sigma~' (day)'),y=unit(-2,'lines'))
  if(column == 1) grid.text('Difference in GI (day)',y=unit(-2,'lines'))
  if(column == 2) grid.text('Difference in SI (day)',y=unit(-2,'lines'))
  if(column == 3) grid.text('Difference in GI (day)',y=unit(-2,'lines'))
  grid.text('Power (%)',x=unit(-3,'lines'),rot=90)
  
  
  if(row==1) grid.text('Diff incub period', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==1 & column==1) grid.text(bquote('& isolate status, '~omega), x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==1 & column==2) grid.text(bquote('& isolate status, '~sigma), x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  if(row==1) grid.text('& isolate status', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==1 & column==2) grid.text('& isolate status', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  
  if(row==2) grid.text('Diff incub,inf period', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==2 & column==1) grid.text(bquote('& prob of inf, '~omega), x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==2 & column==2) grid.text(bquote('& prob of inf, '~sigma), x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  if(row==2) grid.text('& prob of inf', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==2 & column==2) grid.text('& prob of inf', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  
  if(row==3) grid.text('Diff incub, inf period', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==3 & column==1) grid.text(bquote('& isolate status, '~omega), x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==3 & column==2) grid.text(bquote('& isolate status, '~sigma), x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  if(row==3) grid.text('& isolate status', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==3 & column==2) grid.text('& isolate status', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  
  # symbols
  if(row==1 & column==1) {
    
    grid.text('sample size', x=0.05, y=0.475, just='left', default.units = 'native', gp = gpar(fontsize = 8))
    
    grid.points(0.075,0.425,default.units = 'native', pch=16,gp=gpar(col='grey60', cex=0.5))
    grid.points(0.075,0.375,default.units = 'native', pch=17,gp=gpar(col='grey60', cex=0.5))
    
    grid.text('100', x=0.15, y=0.425, just='left', default.units = 'native', gp = gpar(fontsize = 6))
    grid.text('25', x=0.15, y=0.375, just='left', default.units = 'native', gp = gpar(fontsize = 6))
    
  }
  
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  
  
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # x=1.3, y=0.2
  if(row == 1 & column == 1 & type=='gen'){ colourbar.pathogen(x_bottom_left = 0.05, y_bottom_left = 0.7, y_length = 0.050, x_length = 0.6, palette=col.pts.hue, num=1, row=1, column=1, type='gen')}
  if(row == 2 & column == 1 & type=='gen'){ colourbar.pathogen(x_bottom_left = 0.05, y_bottom_left = 0.7, y_length = 0.050, x_length = 0.6, palette=col.pts.hue, num=2, row=2, column=1, type='gen')}
  if(row == 3 & column == 1 & type=='gen'){ colourbar.pathogen(x_bottom_left = 0.05, y_bottom_left = 0.7, y_length = 0.050, x_length = 0.6, palette=col.pts.hue, num=3, row=3, column=1, type='gen')}
  
  
  popViewport()
  popViewport()
  
}

colourbar.pathogen <- function(x_bottom_left = 1.1, y_bottom_left = 0.4, y_length = 0.2, x_length = 0.25, palette=NULL, num=NULL, row=1, column=1, type='gen'){
  
  y_bottom_left_start = c(y_bottom_left, y_bottom_left-1.5*y_length, y_bottom_left-3*y_length)
  
  for(j in 1:length(y_bottom_left_start)){
    y_bottom_left=y_bottom_left_start[j]
    
    colour_hex = palette[[j]]
    cols <- colorRampPalette(colour_hex)
    
    x_start = seq(x_bottom_left, x_bottom_left+x_length,length.out = 126)[-126]
    x_end = seq(x_bottom_left, x_bottom_left+x_length,length.out = 126)[-1]
    
    for(i in 1:125){
      grid.polygon(c(x_start[i],x_end[i],x_end[i],x_start[i]),
                   c(y_bottom_left,y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length),
                   default.units = 'native',gp=gpar(fill=cols(125)[i], col=cols(125)[i]))
    }
    
    
    if(num %in% c(1,3)){
      x_start <-  seq(x_bottom_left, x_bottom_left+x_length,length.out=3)
      label <- c('0', '3.5', '7')
    } 
    if(num ==2){
      x_start <-  seq(x_bottom_left, x_bottom_left+x_length,length.out=3)
      label <- c('0', '1.75', '3.5')
    } 
    
    for(breaks in 1:length(x_start)){
      
      if(breaks %in% c(2:(length(x_start)-1))){
        grid.lines(c(x_start[breaks], x_start[breaks]),
                   c(y_bottom_left, y_bottom_left+(y_length/5)),default.units = 'native',gp=gpar(col='black'))
        grid.lines(c(x_start[breaks], x_start[breaks]),
                   c(y_bottom_left+y_length, y_bottom_left+y_length-(y_length/5)),default.units = 'native',gp=gpar(col='black'))
      }
      
      
      if(j==1) grid.text(label[breaks],x=x_start[breaks],y=y_bottom_left+1.5*y_length,gp = gpar(fontsize = 6), default.units = 'native')
      
    }
    
    grid.polygon(c(x_bottom_left,x_bottom_left,x_bottom_left+x_length,x_bottom_left+x_length),
                 c(y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length,y_bottom_left),
                 default.units = 'native',gp=gpar(fill=NA, col='black'))
    
    if(row%in%c(1,3) & column ==1 & j==1) grid.text('No iso', x=x_bottom_left+x_length+0.13, y=y_bottom_left+(y_length/2), default.units = 'native', gp = gpar(fontsize = 6))
    if(row%in%c(1,3) & column ==1 & j==2) grid.text('Iso D8', x=x_bottom_left+x_length+0.13, y=y_bottom_left+(y_length/2), default.units = 'native', gp = gpar(fontsize = 6))
    if(row%in%c(1,3) & column ==1 & j==3) grid.text('Iso D4', x=x_bottom_left+x_length+0.13, y=y_bottom_left+(y_length/2), default.units = 'native', gp = gpar(fontsize = 6))
    
    if(row==2 & column ==1 & j==1) grid.text('20%', x=x_bottom_left+x_length+0.1, y=y_bottom_left+(y_length/2), default.units = 'native', gp = gpar(fontsize = 6))
    if(row==2 & column ==1 & j==2) grid.text('50%', x=x_bottom_left+x_length+0.1, y=y_bottom_left+(y_length/2), default.units = 'native', gp = gpar(fontsize = 6))
    if(row==2 & column ==1 & j==3) grid.text('80%', x=x_bottom_left+x_length+0.1, y=y_bottom_left+(y_length/2), default.units = 'native', gp = gpar(fontsize = 6))
    
    
  }
  
  # grid.text(bquote('|'~mu[s[2]]-mu[s[1]]~'|'),x=x_bottom_left+0.2, y=y_bottom_left_start[1]+2.5*y_length, 
  #           default.units = 'native', gp = gpar(fontsize = 8))
  grid.text('|Incub difference|',x=x_bottom_left, y=y_bottom_left_start[1]+2.5*y_length, 
            default.units = 'native', just = 'left', gp = gpar(fontsize = 8))

  
}

png('figure/power_pathogen.png',height=24,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=3)))

paneller(1,1,type='gen')
paneller(1,2,type='si')
paneller(1,3,type='dgi')

paneller(2,1,type='gen')
paneller(2,2,type='si')
paneller(2,3,type='dgi')

paneller(3,1,type='gen')
paneller(3,2,type='si')
paneller(3,3,type='dgi')

popViewport()
popViewport()
dev.off()

# rm(paneller)

paneller=function(row = 1,column=1, type=NULL)
{
  xlm=c(0,2); ylm=c(0,1) 
 
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # data
  if(type=='gen') row.pair[,`:=`(diff=diff.mean.gen,
                                 pow=pow.gen)]
  if(type=='si') row.pair[,`:=`(diff=diff.mean.si,
                                pow=pow.si)]
  if(type=='dgi') row.pair[,`:=`(diff=diff.mean.si,
                                pow=pow.dgi)]
  
  if(row == 1) data = row.pair[fig=='d']
  if(row == 2) data = row.pair[fig=='e']
  
  
  # set colour and symbol
  pch.sample.size = c(16,17)
  col.pts.hue = list(c('#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15', '#67000d'),
                     c('#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b'),
                     c('#feec9a', '#fee472', '#fddc49', '#fdd421', '#f2c602', '#caa502', '#a18402'))
  
  col.lines = c('#fb6a4a', '#6baed6', '#fddc49')
  
  col.lines=sapply(1:length(col.lines), function(x) lightup(col.lines[x], alpha=0.25))
  
  val=seq(0,1,length.out = 7)
  
  # gridlines
  for(y.gl in seq(0,1,0.05)){
    grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
    
  for(x.gl in seq(0,2,0.1)){
    grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
  
  # 1: power to detect difference in GI due to changes in hh or nhh contact patterns, all else equal
  # 2: power to detect difference in GI due to changes in hh or hh weekly contact patterns, all else equal
  data[,diff.p.trans:=abs(p.trans.overall.alt-p.trans.overall.ref)]
  size=sort(unique(data$sample.size), decreasing = T)
  
  for(s in 1:length(size)){
    
    data.spline=data[sample.size==size[s]]
    data.points=data[sample.size==size[s] & diff<=1.99 & diff>=0.05 & pow<=0.99]
    
    for(g in 1:3){
      
      # fit spline
      if(row==2 & g==1 & type=='gen'){
        fit=predict(smooth.spline(data.spline[grp==g]$diff, data.spline[grp==g]$pow), seq(0.8,2,0.01))
      } else if(row==2 & g==1 & type%in%c('si','dgi')){
        fit=predict(smooth.spline(data.spline[grp==g]$diff, data.spline[grp==g]$pow), seq(0.9,2,0.01))
      } else if(row==1 & g==3 & s==1 & type=='dgi'){
        fit=predict(smooth.spline(data.spline[grp==g]$diff, data.spline[grp==g]$pow), seq(0,0.4,0.01))
      } else{
        fit=predict(smooth.spline(data.spline[grp==g]$diff, data.spline[grp==g]$pow), seq(0,2,0.01))
      }
      
      
      grid.lines(fit$x, fit$y, default.units = 'native', gp=gpar(col=col.lines[g]))
      
      # assign colour to points
      pal=gradient_n_pal(colours = col.pts.hue[[g]], values = val)
      data.points[grp==g, diff.p.trans.hex:=pal(diff.p.trans)]
      
      
      for(i in 1:data.points[grp==g,.N]){
        
        grid.points(data.points[grp==g]$diff[i],
                    data.points[grp==g]$pow[i],
                    default.units = 'native', pch=pch.sample.size[s],
                    gp=gpar(col=data.points[grp==g]$diff.p.trans.hex[i], cex=0.3))
        
      }
    }
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  grid.xaxis(at=seq(0,2,0.5),label=seq(0,2,0.5),gp=gpar(fontsize=unit(7,'pt')))
  grid.yaxis(at=seq(0,1,0.2),label=seq(0,100,20),gp=gpar(fontsize=unit(7,'pt')))
  
  # top label
  if(row ==1 & column == 1) grid.text('Theoretical GI',y=unit(14,'lines'),gp = gpar(fontsize = 10,fontface='bold'))
  if(row ==1 & column == 2) grid.text('Observed SI',y=unit(14.7,'lines'),gp = gpar(fontsize = 10,fontface='bold'))
  if(row ==1 & column == 2) grid.text('(Lower limit of inferred GI)',y=unit(14,'lines'),gp = gpar(fontsize = 10,fontface='bold'))
  if(row ==1 & column == 3) grid.text('(Upper limit of inferred GI)',y=unit(14,'lines'),gp = gpar(fontsize = 10,fontface='bold'))
  
  # fig labels
  if(row == 1 & column == 1) grid.text('A',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 3) grid.text('C',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 1) grid.text('D',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 2) grid.text('E',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 3) grid.text('F',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 3 & column == 1) grid.text('G',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 3 & column == 2) grid.text('H',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 3 & column == 3) grid.text('I',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  
  
 
  # if(column == 1) grid.text(bquote('Difference in '~omega~' (day)'),y=unit(-2,'lines'))
  # if(column == 2) grid.text(bquote('Difference in '~sigma~' (day)'),y=unit(-2,'lines'))
  if(column == 1) grid.text('Difference in GI (day)',y=unit(-2,'lines'))
  if(column == 2) grid.text('Difference in SI (day)',y=unit(-2,'lines'))
  if(column == 3) grid.text('Difference in GI (day)',y=unit(-2,'lines'))
  
  grid.text('Power (%)',x=unit(-3,'lines'),rot=90)
  
  if(row==1) grid.text('Daily HH vs', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==1 & column==1) grid.text(bquote('non-HH contact, '~omega), x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==1 & column==2) grid.text(bquote('non-HH contact, '~sigma), x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  if(row==1) grid.text('non-HH contact', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==1 & column==2) grid.text('non-HH contact', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  
  if(row==2) grid.text('Daily HH vs', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==2 & column==1) grid.text(bquote('once wkly HH contact, '~omega), x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==2 & column==2) grid.text(bquote('once wkly HH contact, '~sigma), x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  if(row==2) grid.text('once wkly HH contact', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  # if(row==2 & column==2) grid.text('once wkly HH contact, SI', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  
  # symbols
  if(row==1 & column==1) {
    
    grid.text('sample size', x=0.05, y=0.475, just='left', default.units = 'native', gp = gpar(fontsize = 8))
    
    grid.points(0.075,0.425,default.units = 'native', pch=16,gp=gpar(col='grey60', cex=0.5))
    grid.points(0.075,0.375,default.units = 'native', pch=17,gp=gpar(col='grey60', cex=0.5))
    
    grid.text('100', x=0.15, y=0.425, just='left', default.units = 'native', gp = gpar(fontsize = 6))
    grid.text('25', x=0.15, y=0.375, just='left', default.units = 'native', gp = gpar(fontsize = 6))
    
  }
  
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  
  
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # x=1.3, y=0.2
  if(row == 1 & column == 1){ colourbar.contact(x_bottom_left = 0.05, y_bottom_left = 0.7, y_length = 0.050, x_length = 0.6, palette=col.pts.hue, num=1, row=1, column=1, type='gen')}
  if(row == 2 & column == 1){ colourbar.contact(x_bottom_left = 0.05, y_bottom_left = 0.7, y_length = 0.050, x_length = 0.6, palette=col.pts.hue, num=2, row=2, column=1, type='gen')}
  
  
  
  popViewport()
  popViewport()
  
}

colourbar.contact <- function(x_bottom_left = 1.1, y_bottom_left = 0.4, y_length = 0.2, x_length = 0.25, palette=NULL, num=NULL, row=1, column=1, type='gen'){
  
  y_bottom_left_start = c(y_bottom_left, y_bottom_left-1.5*y_length, y_bottom_left-3*y_length)
  
  for(j in 1:length(y_bottom_left_start)){
    y_bottom_left=y_bottom_left_start[j]
    
    colour_hex = palette[[j]]
    cols <- colorRampPalette(colour_hex)
    
    x_start = seq(x_bottom_left, x_bottom_left+x_length,length.out = 126)[-126]
    x_end = seq(x_bottom_left, x_bottom_left+x_length,length.out = 126)[-1]
    
    for(i in 1:125){
      grid.polygon(c(x_start[i],x_end[i],x_end[i],x_start[i]),
                   c(y_bottom_left,y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length),
                   default.units = 'native',gp=gpar(fill=cols(125)[i], col=cols(125)[i]))
    }
    
    
    if(num%in%c(1:2)){
      x_start <-  seq(x_bottom_left, x_bottom_left+x_length,length.out=3)
      label <- c('0', '0.5', '1')
    } 
    
    for(breaks in 1:length(x_start)){
      
      if(breaks %in% c(2:(length(x_start)-1))){
        grid.lines(c(x_start[breaks], x_start[breaks]),
                   c(y_bottom_left, y_bottom_left+(y_length/5)),default.units = 'native',gp=gpar(col='black'))
        grid.lines(c(x_start[breaks], x_start[breaks]),
                   c(y_bottom_left+y_length, y_bottom_left+y_length-(y_length/5)),default.units = 'native',gp=gpar(col='black'))
      }
      
      
      if(j==1) grid.text(label[breaks],x=x_start[breaks],y=y_bottom_left+1.5*y_length,gp = gpar(fontsize = 6), default.units = 'native')
      
    }
    
    grid.polygon(c(x_bottom_left,x_bottom_left,x_bottom_left+x_length,x_bottom_left+x_length),
                 c(y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length,y_bottom_left),
                 default.units = 'native',gp=gpar(fill=NA, col='black'))
    
    if(row%in%c(1,2) & column ==1 & j==1) grid.text('No iso', x=x_bottom_left+x_length+0.13, y=y_bottom_left+(y_length/2), default.units = 'native', gp = gpar(fontsize = 6))
    if(row%in%c(1,2) & column ==1 & j==2) grid.text('Iso D8', x=x_bottom_left+x_length+0.13, y=y_bottom_left+(y_length/2), default.units = 'native', gp = gpar(fontsize = 6))
    if(row%in%c(1,2) & column ==1 & j==3) grid.text('Iso D4', x=x_bottom_left+x_length+0.13, y=y_bottom_left+(y_length/2), default.units = 'native', gp = gpar(fontsize = 6))
    
   
    
  }
  
  if(num%in%c(1,2)) grid.text(bquote('|'~p[2]-p[1]~'|'),x=x_bottom_left+0.18, y=y_bottom_left_start[1]+2.5*y_length, 
                       default.units = 'native', gp = gpar(fontsize = 8))
  
}

png('figure/power_contact.png',height=16,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=3)))

paneller(1,1,type='gen')
paneller(1,2,type='si')
paneller(1,3,type='dgi')

paneller(2,1,type='gen')
paneller(2,2,type='si')
paneller(2,3,type='dgi')


popViewport()
popViewport()
dev.off()



paneller=function(row = 1,column=1, type=NULL)
{
  xlm=c(0,2); ylm=c(0,1) 
  
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # data
  if(type=='gen') row.pair[,`:=`(diff=diff.mean.gen,
                                 pow=pow.gen)]
  if(type=='si') row.pair[,`:=`(diff=diff.mean.si,
                                pow=pow.si)]
  
  if(column==1) data = row.pair[fig=='f']
  if(column==2) data = row.pair[fig=='g']
  
  # set colour and symbol
  pch.sample.size = c(16,17)
  col.pts.hue = list(c('#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15', '#67000d'),
                     c('#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b'),
                     c('#feec9a', '#fee472', '#fddc49', '#fdd421', '#f2c602', '#caa502', '#a18402'))
  
  col.lines = c('#fb6a4a', '#6baed6', '#fddc49')
  val=seq(0,5,length.out = 7)
  
  
  # gridlines
  for(y.gl in seq(0,1,0.05)){
    grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
    
  for(x.gl in seq(0,2,0.1)){
    grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
  

  data[,diff.incub:=abs(mean.incub.alt-mean.incub.ref)]
  size=sort(unique(data$sample.size), decreasing = T)
  
  for(s in 1:length(size)){
    
    data.spline=data[sample.size==size[s]]
    data.points=data[sample.size==size[s] & diff<=1.99 & diff>=0.05 & pow<=0.99]
    
    for(g in 1:3){
      
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
      }
    }
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  grid.xaxis(at=seq(0,2,0.5),label=seq(0,2,0.5),gp=gpar(fontsize=unit(7,'pt')))
  grid.yaxis(at=seq(0,1,0.2),label=seq(0,100,20),gp=gpar(fontsize=unit(7,'pt')))
  
 
  # labels
  if(row == 1 & column == 1) grid.text('A',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 1) grid.text('C',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 2) grid.text('D',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  
  
  if(row == 1) grid.text('Difference in GI (day)',y=unit(-2,'lines'))
  if(row == 2) grid.text('Difference in SI (day)',y=unit(-2,'lines'))
  
  grid.text('Power (%)',x=unit(-3,'lines'),rot=90)
  
  if(row==1 & column==1) grid.text('Diff incub period & isolate status', x=0.05, y=0.95, just='left', default.units = 'native', gp = gpar(fontsize = 8))
  # if(row==1 & column==1) grid.text(bquote('Wild SARS-CoV-2 '~omega), x=0.05, y=0.90, just='left', default.units = 'native', gp = gpar(fontsize = 8))
  if(row==1 & column==1) grid.text('Wild SARS-CoV-2, GI', x=0.05, y=0.90, just='left', default.units = 'native', gp = gpar(fontsize = 8))
  
  if(row==1 & column==2) grid.text('Diff incub & isolate status', x=0.05, y=0.95, just='left', default.units = 'native', gp = gpar(fontsize = 8))
  # if(row==1 & column==2) grid.text(bquote('Ferretti et al model '~omega), x=0.05, y=0.90, just='left', default.units = 'native', gp = gpar(fontsize = 8))
  if(row==1 & column==2) grid.text('Ferretti et al model, GI', x=0.05, y=0.90, just='left', default.units = 'native', gp = gpar(fontsize = 8))
  
  # if(row==2 & column==1) grid.text(bquote('Wild SARS-CoV-2 '~sigma), x=0.05, y=0.95, just='left', default.units = 'native', gp = gpar(fontsize = 8))
  # if(row==2 & column==2) grid.text(bquote('Ferretti et al '~sigma), x=0.05, y=0.95, just='left', default.units = 'native', gp = gpar(fontsize = 8))
  if(row==2 & column==1) grid.text('Wild SARS-CoV-2, SI', x=0.05, y=0.95, just='left', default.units = 'native', gp = gpar(fontsize = 8))
  if(row==2 & column==2) grid.text('Ferretti et al, SI', x=0.05, y=0.95, just='left', default.units = 'native', gp = gpar(fontsize = 8))
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  
  
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(row == 1 & column == 1){ colourbar.pathogen(x_bottom_left = 0.05, y_bottom_left = 0.7, y_length = 0.050, x_length = 0.6, palette=col.pts.hue, num=1)}
  
  
  popViewport()
  popViewport()
  
}

png('figure/power_ferretti_wild_supp.png',height=16,width=16,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=2)))

paneller(1,1,type='gen')
paneller(1,2,type='gen')

paneller(2,1,type='si')
paneller(2,2,type='si')

popViewport()
popViewport()
dev.off()



paneller=function(row = 1,column=1)
{
  xlm=c(0,25); ylm=c(0,0.16)
  
  
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # data
  if(column==1){row.index = c(80, 112); space = c(0,0.0064*2)} #reference, alternate
  if(column==2){row.index = c(280, 312); space = c(0,0.0064*2)}
  if(column==3){row.index = c(480, 508); space = c(0,0.0064*2)}
 
  # set colour and symbol
  col.lines = list(c('#a50f15', lightup('#fb6a4a',0.5)),
                   c('#08519c', lightup('#6baed6',0.5)),
                   c('#caa502', lightup('#fddc49',0.5)))
  
  
  
  # gridlines
  for(y.gl in seq(0,0.16,0.0064)){
      grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
    
  for(x.gl in seq(0,25,1)){
      grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
  

  # hist plots 
  for(i in 1:2){
    data=hist(gen.dist[row.index[i], ], breaks=c(0:25), plot = F)
    
    for(b in 1:(length(data$breaks)-1)){
      
      grid.polygon(c(data$breaks[b], data$breaks[b+1], data$breaks[b+1], data$breaks[b]),
                   c(0,0,data$counts[b]/sum(data$counts),data$counts[b]/sum(data$counts)), default.units = 'native', gp=gpar(fill=col.lines[[column]][i], col=col.lines[[column]][i]))
      
    }
    
    mean=round(mean(gen.dist[row.index[i], ]), 1)
    sd=round(sd(gen.dist[row.index[i], ]), 1)
    # grid.text(bquote(omega[.(i)] == .(mean) ~ ' (' ~ .(sd) ~')'), x=15, y=0.15-space[i], just='left', default.units = 'native', gp = gpar(fontsize = 8))
    grid.text(bquote('GI'[.(i)] == .(mean) ~ ' (' ~ .(sd) ~')'), x=15, y=0.15-space[i], just='left', default.units = 'native', gp = gpar(fontsize = 8))
    
  }
  
  
  
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  grid.xaxis(at=seq(0,25,5),label=seq(0,25,5),gp=gpar(fontsize=unit(7,'pt')))
  grid.yaxis(at=seq(0,0.15,0.05),label=seq(0,15,5),gp=gpar(fontsize=unit(7,'pt')))
  
  
  
  # labels
  if(row == 1 & column == 1) grid.text('A',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 3) grid.text('C',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  

  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}

png('figure/gen_dist_supp.png',height=8,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3)))

paneller(1,1)
paneller(1,2)
paneller(1,3)


popViewport()
popViewport()
dev.off()


# plot schematic
load('output/param_pair_disease/gen.dist.RData')
load('output/param_pair_disease/si.dist.RData')
load('output/param_pair_disease/param.pair.disease.RData')
load('output/param_cluster_schematic/cluster.results.RData')

results.three.gen = results[,max(generation), by=.(sim,set.no)]
results.three.gen = results.three.gen[V1==3,]

set.seed(123)
paneller=function(row = 1,column=1)
{
  
  if(row==1){xlm=c(-0.5,10.5);ylm=c(0.5,4)}
  if(row==2 & column==1) {xlm=c(0,1); ylm=c(-5,25)}
  if(row==2 & column==2) {xlm=c(0,15); ylm=c(0,15)}
  
  innermargins = c(2,2,2,2)
  
  if(row==1) pushViewport(viewport(layout.pos.col=c(1,2),layout.pos.row=row))
  if(row==2) pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # gridlines
  if(row==2 & column==1){
    for(y.gl in seq(-5,25,1.5)){
      grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
    for(x.gl in seq(0,1,0.05)){
      grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
  }
  
  if(row==2 & column==2){
    for(y.gl in seq(0,15,1)){
      grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
    for(x.gl in seq(0,15,1)){
      grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
  }
  
  if(row==2){
    # set colour 
    col.hue = c('#8da0cb', '#a6d854', '#fc8d62', '#ffd92f', '#e5c494', '#66c2a5', '#e78ac3')
    shape.pts =c(15,16,17,3,4,18,25)
    
    # incubation period, serial interval, attack rate of each pathogen
    data.diseases = data.table(pathogen=c('covid.wild', 'covid.delta', 'smallpox', 'measles', 'h1n1pdm09', 'sars', 'ebola'),
                               inf.func.name=c('spline.covid.wild', 'spline.covid.delta', 'spline.smallpox', 'spline.measles',  'spline.flu', NA, NA),
                               inc.func.name=c('lnorm', 'weibull', 'lnorm', 'lnorm', 'lnorm', 'lnorm', 'gamma'),
                               inc.pars.1=c(1.62, 1.83, 2.510476, 2.525726, 1.42947, 1.3862935, 5.53404),
                               inc.pars.2=c(0.42, 4.93, 0.1746595, 0.2097952, 0.241453, 0.5762278, 0.921569),
                               si.func.name=c('normal','normal','normal','normal','normal','normal', 'gamma'),
                               si.pars.1=c(7.8, 4.1, 16, 11.7, 2.8, 8.4, 2.706555671),
                               si.pars.2=c(5.2, 2.8, 4, 0.5, 0.5, 3.8, 5.652941176),
                               si.min = c(3.8,0,10,9,1.7,NA,NA), 
                               si.max = c(8,7,17.7,13,3.7,NA,NA),
                               ar.pars.1=c(0.132,0.230,0.6,0.8,0.11,NA,NA),
                               ar.pars.2=c(0.182,0.373,0.9,0.9,0.18,NA,NA))
    
    # serial interval, onset to isolation in Ali et al for first three data, Lesser for last three data
    data.onset.to.iso = data.table(period=c('pre-peak','peak','post-peak', 'jan-feb', 'jan-feb', 'jan-feb'),
                                   si.min=c(5.9,4.5,2.0, 3,5.3,6.2),
                                   si.max=c(8.5,6.0,3.9, 4.2,11.0,9.7),
                                   si.mean=c(7.2,5.2,3.0, 3.6,8.1,8.0),
                                   iso.min=c(5.54,1.90,0.899, 0,3,6),
                                   iso.max=c(10.9,7.6,5.764, 2,5,12.2),
                                   iso.mean=c(8.34,4.75,3.315, 1,4,9.1))
  }

  if(row==2 & column==1){
    
    legend.symbol.y = c(1,-0.25,-1.5,-2.75,-4)
    
    for(i in 1:5){
      
      row.index = which(param$inf.func.name==data.diseases$inf.func.name[i] & is.na(param$mean.iso) & param$fig=='a')
      
      ar.param = extract_param(values=c(data.diseases[i]$ar.pars.1,data.diseases[i]$ar.pars.2))
      ar.x = rbeta(150, shape1 = ar.param[1], shape2 = ar.param[2])
      si.y = random.generate(data.diseases[i,si.func.name], val=c(data.diseases[i,si.pars.1],data.diseases[i,si.pars.2])) 
      si.y = si.y(150)
      
      # trim extreme points
      rm.pts = which(si.y>ylm[2] | si.y<ylm[1])
      if(length(rm.pts)!=0){
        si.y = si.y[-rm.pts]
        ar.x = ar.x[-rm.pts]
      }
      
      grid.points(ar.x, si.y, default.units = 'native', pch=shape.pts[i], gp=gpar(col=lightup(col.hue[i],0.5), fill=lightup(col.hue[i],0.5), cex=0.3))
      
      
      fit=predict(smooth.spline(param[row.index,]$p.trans.overall, apply(si.dist[row.index,],1,median), nknots=4), seq(0,1,0.1))
      grid.lines(fit$x, fit$y, default.units = 'native', gp=gpar(col=col.hue[i]))
     
      
      # legend
      grid.points(0.6, legend.symbol.y[i], default.units = 'native', pch=shape.pts[i], gp=gpar(col=lightup(col.hue[i],0.5), fill=lightup(col.hue[i],0.5), cex=0.3))
      if(i==1) grid.text('SARS-CoV-2 wild', x=0.63, y=legend.symbol.y[i], just='left', default.units = 'native', gp = gpar(fontsize = 6))
      if(i==2) grid.text('SARS-CoV-2 delta', x=0.63, y=legend.symbol.y[i], just='left', default.units = 'native', gp = gpar(fontsize = 6))
      if(i==3) grid.text('Smallpox', x=0.63, y=legend.symbol.y[i], just='left', default.units = 'native', gp = gpar(fontsize = 6))
      if(i==4) grid.text('Measles', x=0.63, y=legend.symbol.y[i], just='left', default.units = 'native', gp = gpar(fontsize = 6))
      if(i==5) grid.text('H1N1pdm09', x=0.63, y=legend.symbol.y[i], just='left', default.units = 'native', gp = gpar(fontsize = 6))
      
      
    }
  }  
    
  if(row==2 & column==2){
    
    
    for(i in 1:2){
      
      row.index = which(param$inf.func.name==data.diseases$inf.func.name[i] & param$fig=='b')
      
      fit.lwr=predict(smooth.spline(param[row.index,]$mean.iso, apply(si.dist[row.index,],1,quantile, probs=0.25), nknots=4), seq(0,15,1))
      fit.upp=predict(smooth.spline(param[row.index,]$mean.iso, apply(si.dist[row.index,],1,quantile, probs=0.75), nknots=4), seq(0,15,1))
      
      grid.polygon(c(fit.lwr$x, rev(fit.upp$x)),
                   c(fit.lwr$y, rev(fit.upp$y)),
                   default.units = 'native',gp=gpar(fill=lightup(col.hue[i], 0.2),col=NA))
      
      fit=predict(smooth.spline(param[row.index,]$mean.iso, apply(si.dist[row.index,],1,median), nknots=4), seq(0,15,1))
      grid.lines(fit$x, fit$y, default.units = 'native', gp=gpar(col=col.hue[i]))
   
      if(i==1) {
        grid.lines(c(8,9), c(2,2), default.units = 'native', gp=gpar(col=col.hue[i]))
        grid.text('SARS-CoV-2 wild', x=9.25, y=2, just='left', default.units = 'native', gp = gpar(fontsize = 6))
        
      }
      if(i==2) {
        grid.lines(c(8,9), c(1,1), default.units = 'native', gp=gpar(col=col.hue[i]))
        grid.text('SARS-CoV-2 Delta', x=9.25, y=1, just='left', default.units = 'native', gp = gpar(fontsize = 6))
        
      }
      
      if(i==1){
        
        legend.symbol.y=c(13.75,13,12.25,11.5)
        grid.text('Time periods', x=0.85, y=14.5, just='left', default.units = 'native', gp = gpar(fontsize = 6))
        
        
        for(j in 1:6){
          if(j %in% c(1:3)) grid.points(data.onset.to.iso$iso.mean[j], data.onset.to.iso$si.mean[j], default.units = 'native', pch=shape.pts[j], gp=gpar(col=col.hue[i], fill=col.hue[i], cex=0.5))
          if(j %in% c(4:6)) grid.points(data.onset.to.iso$iso.mean[j], data.onset.to.iso$si.mean[j], default.units = 'native', pch=4, gp=gpar(col=col.hue[i], fill=col.hue[i], cex=0.5))
          
          grid.lines(c(data.onset.to.iso$iso.min[j], data.onset.to.iso$iso.max[j]), c(data.onset.to.iso$si.mean[j], data.onset.to.iso$si.mean[j]), default.units = 'native', gp=gpar(col=col.hue[i]))
          grid.lines(c(data.onset.to.iso$iso.mean[j], data.onset.to.iso$iso.mean[j]), c(data.onset.to.iso$si.min[j], data.onset.to.iso$si.max[j]), default.units = 'native', gp=gpar(col=col.hue[i]))
        
          if(j %in% c(1:3)) grid.points(1, legend.symbol.y[j], default.units = 'native', pch=shape.pts[j], gp=gpar(col=col.hue[i], fill=col.hue[i], cex=0.5))
          if(j == 4) grid.points(1, legend.symbol.y[j], default.units = 'native', pch=4, gp=gpar(col=col.hue[i], fill=col.hue[i], cex=0.5))
          
          if(j==1) grid.text('pre-peak, Jan 9-22, 2020 (Ali et al.)', x=1.5, y=legend.symbol.y[j], just='left', default.units = 'native', gp = gpar(fontsize = 6))
          if(j==2) grid.text('peak, Jan 23-29, 2020 (Ali et al.)', x=1.5, y=legend.symbol.y[j], just='left', default.units = 'native', gp = gpar(fontsize = 6))
          if(j==3) grid.text('post-peak, Jan 30-Feb 13, 2020 (Ali et al.)', x=1.5, y=legend.symbol.y[j], just='left', default.units = 'native', gp = gpar(fontsize = 6))
          if(j==4) grid.text('Jan 14-Feb 12, 2020 (Bi et al.)', x=1.5, y=legend.symbol.y[j], just='left', default.units = 'native', gp = gpar(fontsize = 6))
          
        }
        
        
        
        
      }
    
    }
  }  
  
  if(row==1){
    # case data 9, 10, 11, 16, 19
    case.data = results[sim==results.three.gen$sim[19] & set.no==results.three.gen$set.no[19]]
    case.data
    setorder(case.data, exposure)
    
    # period of overlap
    grid.polygon(c(case.data$exposure[2], case.data$exposure[3], case.data$exposure[3], case.data$exposure[2]),
                 c(0.5,0.5,4,4), default.units = 'native', gp=gpar(fill='#FAFAFA', col='#F5F5F5'))
    
    grid.text('competing infectors', x=case.data$exposure[2] + ((case.data$exposure[3]-case.data$exposure[2])/2), 
              y=1.5, default.units = 'native', just='centre', gp = gpar(fontsize = 6))
    
    
    # infectiousness period distribution
    inf.func.name='spline.covid.delta'
    inf.f = inf.func[[inf.func.name]]
    max.inf.day=30
    dt=5/(24*60)
    
    col.hue = list(c('#fdd0a2', '#fdae6b', '#fd8d3c', '#f16913', '#d94801', '#a63603', '#7f2704'),
                   c('#bfd3e6', '#9ebcda', '#8c96c6', '#8c6bb1', '#88419d', '#810f7c', '#4d004b'),
                   c('#d9f0a3', '#addd8e', '#78c679', '#41ab5d', '#238443', '#006837', '#004529'))
    val=seq(0,0.0005,length.out=7)
    
    
    for(i in 1:3){
      
      timestep = 1:(max.inf.day/dt)*dt
      
      index = case.data$inc[i]+18
      index = which(timestep< (case.data$inc[i]+18))
      
      case.inf = inf.f(case.data$inc[i], max.inf.day, dt)*0.0005
      case.inf = case.inf[index]
      
      timestep = timestep[index]
      timestep = timestep + case.data$exposure[i]
      
      pal=gradient_n_pal(colours = col.hue[[i]], values = val)
      case.inf.hex = pal(case.inf)
      
      start=seq(9.5,10.5,0.25)
      stop=seq(9.625,10.625,0.25)
      count=1
      
      for(j in 1:length(timestep)){ 
        
        if(timestep[j]<=9.5){
          if(j==1){
            grid.lines(c(timestep[j]-dt,timestep[j]), c(ylm[2]-i,ylm[2]-i), 
                       default.units = 'native', gp=gpar(col=case.inf.hex[j], lwd=2))
          } else{
            grid.lines(c(timestep[j-1],timestep[j]), c(ylm[2]-i,ylm[2]-i), 
                       default.units = 'native', gp=gpar(col=case.inf.hex[j], lwd=2))
          }
        }else if(timestep[j]<=9.5 & timestep[j+1]>=9.5){
          grid.lines(c(timestep[j],9.5), c(ylm[2]-i,ylm[2]-i), 
                     default.units = 'native', gp=gpar(col=case.inf.hex[j], lwd=2))
          
        }
      
        
        if(timestep[j]>=9.5 & timestep[j]<=10.48){
          
          if(timestep[j]>=start[count] & timestep[j]<=stop[count]){
            grid.lines(c(timestep[j-1],timestep[j]), c(ylm[2]-i,ylm[2]-i), 
                       default.units = 'native', gp=gpar(col=case.inf.hex[j], lwd=2))
          }else if(timestep[j]>=stop[count]){
            count=count+1
          }
        }
  
      }
      
      grid.text(bquote(~i[.(i)]), x=case.data$exposure[i]-0.12, y=ylm[2]-i+0.15, just='left', default.units = 'native', gp = gpar(fontsize = 6))
      grid.text(bquote(~o[.(i)]), x=case.data$onset[i]-0.12, y=ylm[2]-i+0.15, just='left', default.units = 'native', gp = gpar(fontsize = 6))
      
      
      # serial interval, generation interval, incubation period
      if(i==1){
        
        grid.lines(c(case.data$exposure[i], case.data$exposure[i]+0.35), c(ylm[2]-i+0.5,ylm[2]-i+0.5),
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$onset[i], case.data$onset[i]-0.35), c(ylm[2]-i+0.5,ylm[2]-i+0.5),
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$exposure[i], case.data$exposure[i]), c(ylm[2]-i+0.45,ylm[2]-i+0.55), 
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$onset[i], case.data$onset[i]), c(ylm[2]-i+0.45,ylm[2]-i+0.55), 
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.text('incubation', x=case.data$onset[i]/2, y=ylm[2]-i+0.6, just='centre', default.units = 'native', gp = gpar(fontsize = 6))
        grid.text('period, s', x=case.data$onset[i]/2, y=ylm[2]-i+0.4, just='centre', default.units = 'native', gp = gpar(fontsize = 6))
        
        grid.lines(c(case.data$onset[i], case.data$onset[i]+0.4), c(ylm[2]-i+0.5,ylm[2]-i+0.5),
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$exposure[i+1], case.data$exposure[i+1]-0.4), c(ylm[2]-i+0.5,ylm[2]-i+0.5),
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$onset[i], case.data$onset[i]), c(ylm[2]-i+0.45,ylm[2]-i+0.55), 
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$exposure[i+1], case.data$exposure[i+1]), c(ylm[2]-i+0.45,ylm[2]-i+0.55), 
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.text('onset to transmission,', x=case.data$onset[i]+(case.data$exposure[i+1]-case.data$onset[i])/2, y=ylm[2]-i+0.7, just='centre', default.units = 'native', gp = gpar(fontsize = 6))
        grid.text(bquote(tau), x=case.data$onset[i]+(case.data$exposure[i+1]-case.data$onset[i])/2, y=ylm[2]-i+0.5, just='centre', default.units = 'native', gp = gpar(fontsize = 6))
        
        
        
        grid.lines(c(case.data$exposure[i], case.data$exposure[i]+0.43), c(ylm[2]-i-0.33,ylm[2]-i-0.33),
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$exposure[i+1], case.data$exposure[i+1]-0.43), c(ylm[2]-i-0.33,ylm[2]-i-0.33),
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$exposure[i], case.data$exposure[i]), c(ylm[2]-i-0.38,ylm[2]-i-0.28), 
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$exposure[i+1], case.data$exposure[i+1]), c(ylm[2]-i-0.38,ylm[2]-i-0.28), 
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.text(bquote('generation interval, '~omega), x=case.data[i]$exposure+(case.data[i+1]$exposure-case.data[i]$exposure)/2, y=ylm[2]-i-0.33, just='centre', default.units = 'native', gp = gpar(fontsize = 6))
        
        
        grid.lines(c(case.data$onset[i], case.data$onset[i]+0.75), c(ylm[2]-i-0.66,ylm[2]-i-0.66),
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$onset[i+1], case.data$onset[i+1]-0.75), c(ylm[2]-i-0.66,ylm[2]-i-0.66),
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$onset[i], case.data$onset[i]), c(ylm[2]-i-0.71,ylm[2]-i-0.61), 
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.lines(c(case.data$onset[i+1], case.data$onset[i+1]), c(ylm[2]-i-0.71,ylm[2]-i-0.61), 
                   default.units = 'native', gp=gpar(col='black', lwd=1))
        grid.text(bquote('serial interval, '~sigma), x=case.data[i]$onset+(case.data[i+1]$onset-case.data[i]$onset)/2, y=ylm[2]-i-0.66, just='centre', default.units = 'native', gp = gpar(fontsize = 6))
        
        
      }
    }
  }
   
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # axis
  if(row==2 & column==1){
    grid.xaxis(at=seq(0,1,0.25),label=seq(0,1,0.25),gp=gpar(fontsize=unit(7,'pt')))
    grid.yaxis(at=seq(-5,25,5),label=seq(-5,25,5),gp=gpar(fontsize=unit(7,'pt')))
  }
  if(row==2 & column==2){
    grid.xaxis(at=seq(0,15,5),label=seq(0,15,5),gp=gpar(fontsize=unit(7,'pt')))
    grid.yaxis(at=seq(0,15,5),label=seq(0,15,5),gp=gpar(fontsize=unit(7,'pt')))
  }
  if(row==1){
    grid.xaxis(at=seq(0,10,1),label=seq(0,10,1),gp=gpar(fontsize=unit(7,'pt')))
  }
  
 
  # labels
  if(row==2 & column == 1) {
    grid.text('B',x=unit(-2,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(12,'pt')))
    grid.text('Serial interval (day)',x=unit(-2.5,'lines'),rot=90)
    grid.text('Infection probability (%)',y=unit(-2,'lines'))
  }
  if(row==2 & column == 2) {
    grid.text('C',x=unit(-2,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(12,'pt')))
    grid.text('Serial interval (day)',x=unit(-2.5,'lines'),rot=90)
    grid.text('Onset to isolation (days)',y=unit(-2,'lines'))
  }
  if(row == 1) {
    grid.text('A',x=unit(-2,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(12,'pt')))
    grid.text('Day',y=unit(-2,'lines'))
    grid.text('Cases',x=unit(-2,'lines'),rot=90)
  }
  
  
   
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  
  
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  popViewport()
  popViewport()
  
}


extract_param <- function(type,values,distribution){
  
  # DEBUG:
  # type = "median"; values = c(0.6,0.7); distribution = "lnorm"
  
  # Validate inputs
  #if(type=="median" & length(values)!=3){stop("Need 'values' to be vector length 3")}
  if(length(values)!=2){stop("Need 'values' to be vector length 2")}
  
  # Extract distribution parameters using optimisation
  param = c(a = 1, b = 1)
  result2 = optim(param, fit_function, method="L-BFGS-B", val=values, lower=c(0,0), hessian=FALSE) #, control=list(trace=1))
  
  # Output parameters
  result2$par
  
}



fit_function <- function(param, val) {
  
  (pbeta(val[1], shape1 = param[["a"]],shape2 = param[["b"]]) - 0.025)^2 + (pbeta(val[2], shape1 = param[["a"]],shape2 = param[["b"]]) - 0.975)^2
  
}

pdf.generate <- function(func.name, val){
  
  if(func.name=='normal') return(function(x) dnorm(x, mean = val[1], sd = val[2]))
  if(func.name=='lnorm') return(function(x) dlnorm(x, meanlog = val[1], sdlog = val[2]))
  if(func.name=='weibull') return(function(x) dweibull(x, shape = val[1], scale = val[2]))
  if(func.name=='gamma') return(function(x) dgamma(x, shape = val[1], scale = val[2]))
  
}

random.generate <- function(func.name, val){
  
  if(func.name=='normal') return(function(x) rnorm(x, mean = val[1], sd = val[2]))
  if(func.name=='lnorm') return(function(x) rlnorm(x, meanlog = val[1], sdlog = val[2]))
  if(func.name=='weibull') return(function(x) rweibull(x, shape = val[1], scale = val[2]))
  if(func.name=='gamma') return(function(x) rgamma(x, shape = val[1], scale = val[2]))
  
  
}


png('figure/schematic.png',height=16,width=16,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=2)))

paneller(1,1)

paneller(2,1)
paneller(2,2)

popViewport()
popViewport()
dev.off()





paneller=function(row = 1,column=1)
{
  xlm=c(0,30);ylm=c(0,0.002)
  
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  max.inf.day=30
  dt=5/(24*60)
  time=1:(max.inf.day/dt)*dt
  
  beta=c(1e-05,1e-04,1e-03,1e-02)
  inc=rep(4.145, 4)
  
  inf.func.name=c('spline.covid.delta')
  inf.f = inf.func[[inf.func.name]]
  inf.matrix = sapply(1:length(beta), function(x){ inf.f(inc[x], max.inf.day, dt) })
  
  inf.matrix = sapply(1:length(beta), function(x){ inf.matrix[,x]*beta[x] })
  
  p.inf = sapply(1:length(beta), function(x){
    # P(not infected)
    p.non.inf = exp(-sum(inf.matrix[,x]))
    # P(infected at interval | not infected previously)
    p.inf.interval = (1-exp(-inf.matrix[,x]))*c(1,exp(-cumsum(inf.matrix[,x])[1:(length(inf.matrix[,x])-1)]))
    
    
    p.inf.interval = p.inf.interval/(1-p.non.inf)
    
    return(c(p.inf.interval))  
  })
  
  col.line = c('#8da0cb', '#a6d854', '#fc8d62', '#ffd92f', '#e5c494', '#66c2a5', '#e78ac3')
  
  for(i in 1:length(beta)){
    grid.lines(time, p.inf[,i], default.units = 'native', gp=gpar(col=col.line[i], lwd=2))
    
    y.label = 0.001875-(4-i)*0.000125
    beta.label = c('0.00001','0.0001','0.001','0.01')
    
    grid.lines(c(20,21), c(y.label,y.label), default.units = 'native', gp=gpar(col=col.line[i]))
    grid.text(bquote(beta == .(beta.label[i])), x=21.5, y=y.label, just='left', default.units = 'native', gp = gpar(fontsize = 6))

   
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # axis
  grid.xaxis(at=seq(0,30,5),label=seq(0,30,5),gp=gpar(fontsize=unit(7,'pt')))
  grid.yaxis(at=seq(0,0.002,0.0005),label=c('0', '5e-04', '1e-03', '1.5e-03', '2e-03'),gp=gpar(fontsize=unit(7,'pt')))
  
  
  # labels
  grid.text('Probability (%)',x=unit(-3,'lines'),rot=90)
  grid.text('Time (unit)',y=unit(-2,'lines'))
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  
  
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  popViewport()
  popViewport()
  
}

png('figure/approx_supp.png',height=8,width=8,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=1)))

paneller(1,1)

popViewport()
popViewport()
dev.off()


# plot schematic
load('output/20230227/gen.dist.RData')
load('output/20230227/si.dist.RData')
load('output/20230227/param.pair.disease.RData')

paneller=function(row = 1,column=1)
{
  xlm=c(0,15); ylm=c(0,15)
  
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # gridlines
  for(y.gl in seq(0,15,1)){
    grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
  for(x.gl in seq(0,15,1)){
    grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
  }
  
  
  # set colour 
  col.hue = c('#8da0cb', '#a6d854', '#fc8d62', '#ffd92f', '#e5c494', '#66c2a5', '#e78ac3')
  shape.pts =c(15,16,17,3,4,18,25)
    
  # incubation period, serial interval, attack rate of each pathogen
  data.diseases = data.table(pathogen=c('covid.wild', 'covid.delta', 'smallpox', 'measles', 'h1n1pdm09', 'sars', 'ebola'),
                             inf.func.name=c('spline.covid.wild', 'spline.covid.delta', 'spline.smallpox', 'spline.measles',  'spline.flu', NA, NA),
                             inc.func.name=c('lnorm', 'weibull', 'lnorm', 'lnorm', 'lnorm', 'lnorm', 'gamma'),
                             inc.pars.1=c(1.62, 1.83, 2.510476, 2.525726, 1.42947, 1.3862935, 5.53404),
                             inc.pars.2=c(0.42, 4.93, 0.1746595, 0.2097952, 0.241453, 0.5762278, 0.921569),
                             si.func.name=c('normal','normal','normal','normal','normal','normal', 'gamma'),
                             si.pars.1=c(7.8, 4.1, 16, 11.7, 2.8, 8.4, 2.706555671),
                             si.pars.2=c(5.2, 2.8, 4, 0.5, 0.5, 3.8, 5.652941176),
                             si.min = c(3.8,0,10,9,1.7,NA,NA), 
                             si.max = c(8,7,17.7,13,3.7,NA,NA),
                             ar.pars.1=c(0.132,0.230,0.6,0.8,0.11,NA,NA),
                             ar.pars.2=c(0.182,0.373,0.9,0.9,0.18,NA,NA))
    
  # serial interval, onset to isolation in Ali et al for first three data, Lesser for last three data
  data.onset.to.iso = data.table(period=c('pre-peak','peak','post-peak', 'jan-feb', 'jan-feb', 'jan-feb'),
                                 si.min=c(5.9,4.5,2.0, 3,5.3,6.2),
                                 si.max=c(8.5,6.0,3.9, 4.2,11.0,9.7),
                                 si.mean=c(7.2,5.2,3.0, 3.6,8.1,8.0),
                                 iso.min=c(5.54,1.90,0.899, 0,3,6),
                                 iso.max=c(10.9,7.6,5.764, 2,5,12.2),
                                 iso.mean=c(8.34,4.75,3.315, 1,4,9.1))
  
  
  for(i in 1:2){
    
    row.index = which(param$inf.func.name==data.diseases$inf.func.name[i] & param$fig=='b')
      
    fit.lwr=predict(smooth.spline(param[row.index,]$mean.iso, apply(si.dist[row.index,],1,quantile, probs=0.25), nknots=4), seq(0,15,1))
    fit.upp=predict(smooth.spline(param[row.index,]$mean.iso, apply(si.dist[row.index,],1,quantile, probs=0.75), nknots=4), seq(0,15,1))
      
    grid.polygon(c(fit.lwr$x, rev(fit.upp$x)),
                 c(fit.lwr$y, rev(fit.upp$y)),
                 default.units = 'native',gp=gpar(fill=lightup(col.hue[i], 0.2),col=NA))
      
    fit=predict(smooth.spline(param[row.index,]$mean.iso, apply(si.dist[row.index,],1,median), nknots=4), seq(0,15,1))
    grid.lines(fit$x, fit$y, default.units = 'native', gp=gpar(col=col.hue[i]))
      
    if(i==1) {
      grid.lines(c(8,9), c(2,2), default.units = 'native', gp=gpar(col=col.hue[i]))
      grid.text('SARS-CoV-2 wild', x=9.25, y=2, just='left', default.units = 'native', gp = gpar(fontsize = 6))
        
    }
    if(i==2) {
      grid.lines(c(8,9), c(1,1), default.units = 'native', gp=gpar(col=col.hue[i]))
      grid.text('SARS-CoV-2 Delta', x=9.25, y=1, just='left', default.units = 'native', gp = gpar(fontsize = 6))
        
    }
      
    if(i==1){
        
      legend.symbol.y=c(13.75,13,12.25,11.5)
      grid.text('Time periods', x=0.85, y=14.5, just='left', default.units = 'native', gp = gpar(fontsize = 6))
        
        
      for(j in 1:6){
        if(j %in% c(1:3)) grid.points(data.onset.to.iso$iso.mean[j], data.onset.to.iso$si.mean[j], default.units = 'native', pch=shape.pts[j], gp=gpar(col=col.hue[i], fill=col.hue[i], cex=0.5))
        if(j %in% c(4:6)) grid.points(data.onset.to.iso$iso.mean[j], data.onset.to.iso$si.mean[j], default.units = 'native', pch=4, gp=gpar(col=col.hue[i], fill=col.hue[i], cex=0.5))
          
        grid.lines(c(data.onset.to.iso$iso.min[j], data.onset.to.iso$iso.max[j]), c(data.onset.to.iso$si.mean[j], data.onset.to.iso$si.mean[j]), default.units = 'native', gp=gpar(col=col.hue[i]))
        grid.lines(c(data.onset.to.iso$iso.mean[j], data.onset.to.iso$iso.mean[j]), c(data.onset.to.iso$si.min[j], data.onset.to.iso$si.max[j]), default.units = 'native', gp=gpar(col=col.hue[i]))
          
        if(j %in% c(1:3)) grid.points(1, legend.symbol.y[j], default.units = 'native', pch=shape.pts[j], gp=gpar(col=col.hue[i], fill=col.hue[i], cex=0.5))
        if(j == 4) grid.points(1, legend.symbol.y[j], default.units = 'native', pch=4, gp=gpar(col=col.hue[i], fill=col.hue[i], cex=0.5))
          
        if(j==1) grid.text('pre-peak, Jan 9-22, 2020 (Ali et al.)', x=1.5, y=legend.symbol.y[j], just='left', default.units = 'native', gp = gpar(fontsize = 6))
        if(j==2) grid.text('peak, Jan 23-29, 2020 (Ali et al.)', x=1.5, y=legend.symbol.y[j], just='left', default.units = 'native', gp = gpar(fontsize = 6))
        if(j==3) grid.text('post-peak, Jan 30-Feb 13, 2020 (Ali et al.)', x=1.5, y=legend.symbol.y[j], just='left', default.units = 'native', gp = gpar(fontsize = 6))
        if(j==4) grid.text('Jan 14-Feb 12, 2020 (Bi et al.)', x=1.5, y=legend.symbol.y[j], just='left', default.units = 'native', gp = gpar(fontsize = 6))
          
      }
      
        
      }
      
    }
   

  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # axis
  grid.xaxis(at=seq(0,15,5),label=seq(0,15,5),gp=gpar(fontsize=unit(7,'pt')))
  grid.yaxis(at=seq(0,15,5),label=seq(0,15,5),gp=gpar(fontsize=unit(7,'pt')))
  
  
  # labels
  # grid.text('C',x=unit(-2,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  grid.text('Serial interval (day)',x=unit(-2.5,'lines'),rot=90)
  grid.text('Onset to isolation (days)',y=unit(-2,'lines'))
  
  
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  
  
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  popViewport()
  popViewport()
  
}


png('figure/delta_wild_supp.png',height=8,width=8,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=1)))

paneller()

popViewport()
popViewport()
dev.off()
