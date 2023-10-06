source('code/step0_load_library.R')

lightup = function(c, alpha)
{
  z=col2rgb(c)/255
  return(rgb(z[1],z[2],z[3],alpha))  # 0.125 for col3
}

paneller=function(row = 1,column=1, type=NULL)
{
  if(row==1){xlm=c(0,25); ylm=c(0,0.2)}
  if(row==2){xlm=c(0,2); ylm=c(0,1)}
  
  
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # data
  if(column==1){
    # is.na(mean.iso) & is.na(var.iso)
    row.index.ref = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso == 4 & var.incub==5 & mean.incub==4 & growth_rate==log(2)/3.5, set.no]
    row.index.alt = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso == 4 & var.incub==5 & mean.incub==5 & growth_rate==-log(2)/3.5, set.no]
    row.index = c(row.index.ref, row.index.alt); space = c(0,0.008*2) #reference, alternate
  }
  if(column==2){
    row.index.ref = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso == 4 & var.incub==5 & mean.incub==4 & growth_rate==log(2)/3.5, set.no]
    row.index.alt = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso == 4 & var.incub==5 & mean.incub==5 & is.na(growth_rate), set.no]
    row.index = c(row.index.ref, row.index.alt); space = c(0,0.008*2) #reference, alternate
  }
  if(column==3){
    row.index.ref = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso == 4 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
    row.index.alt = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso == 4 & var.incub==5 & mean.incub==5 & growth_rate==log(2)/3.5, set.no]
    row.index = c(row.index.ref, row.index.alt); space = c(0,0.008*2) #reference, alternate
  }
  
  # power 
  if(row==2){
    if(type=='gen') row.pair[,`:=`(diff=diff.mean.gen,
                                   pow=pow.gen)]
    if(type=='si') row.pair[,`:=`(diff=diff.mean.si,
                                  pow=pow.si)]
  }
  
  if(row==2 & column==1) data = row.pair[fig=='a']
  if(row==2 & column==2) data = row.pair[fig=='b']
  if(row==2 & column==3) data = row.pair[fig=='c']
  
  # set colour and symbol
  if(row==1){
    col.lines = list(c('#a50f15', lightup('#fb6a4a',0.5)),
                     c('#08519c', lightup('#6baed6',0.5)),
                     c('#caa502', lightup('#fddc49',0.5)))
  }
  
  if(row==2){
    pch.sample.size = c(16,17)
    col.pts.hue = list(c('#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15', '#67000d'),
                       c('#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b'),
                       c('#feec9a', '#fee472', '#fddc49', '#fdd421', '#f2c602', '#caa502', '#a18402'))
    
    col.lines = c('#fb6a4a', '#6baed6', '#fddc49')
    
    col.lines=sapply(1:length(col.lines), function(x) lightup(col.lines[x], alpha=0.25))
    val=seq(0,7,length.out = 7)
    
  }
  
  # gridlines
  if(row==1){
    for(y.gl in seq(0,0.2,0.008)){
      grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
    
    for(x.gl in seq(0,25,1)){
      grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
  }
  
  if(row==2){
    for(y.gl in seq(0,1,0.05)){
      grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
    
    for(x.gl in seq(0,2,0.1)){
      grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
  }
  
  
  # hist plots 
  if(row==1){
    for(i in 1:2){
      data=hist(gen.dist[row.index[i], ], breaks=c(0:30), plot = F)
      
      for(b in 1:(length(data$breaks)-1)){
        if(data$breaks[b+1]<=25){
          grid.polygon(c(data$breaks[b], data$breaks[b+1], data$breaks[b+1], data$breaks[b]),
                       c(0,0,data$counts[b]/sum(data$counts),data$counts[b]/sum(data$counts)), default.units = 'native', gp=gpar(fill=col.lines[[column]][i], col=col.lines[[column]][i]))
          
        }
      }
      
      mean=round(mean(gen.dist[row.index[i], ]), 1)
      sd=round(sd(gen.dist[row.index[i], ]), 1)
      grid.text(bquote(omega[.(i)] == .(mean) ~ ' (' ~ .(sd) ~')'), x=15, y=0.19-space[i], just='left', default.units = 'native', gp = gpar(fontsize = 8))
      
    }
  }
  
  # power plots 
  if(row==2){
    data[,diff.incub:=abs(mean.incub.alt-mean.incub.ref)]
    size=sort(unique(data$sample.size), decreasing = T)
    
    for(s in 1:length(size)){
      
      data.spline=data[sample.size==size[s]]
      data.points=data[sample.size==size[s] & diff<=1.99 & diff>=0.05 & pow<=0.99]
      
      for(g in 1:3){
        
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
          
        }
      }
    }
  }
  
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(row==1){
    grid.xaxis(at=seq(0,25,5),label=seq(0,25,5),gp=gpar(fontsize=unit(7,'pt')))
    grid.yaxis(at=seq(0,0.2,0.05),label=seq(0,20,5),gp=gpar(fontsize=unit(7,'pt')))
  }
  if(row==2){
    grid.xaxis(at=seq(0,2,0.5),label=seq(0,2,0.5),gp=gpar(fontsize=unit(7,'pt')))
    grid.yaxis(at=seq(0,1,0.2),label=seq(0,100,20),gp=gpar(fontsize=unit(7,'pt')))
  }
  
  
  
  # labels
  if(row == 1 & column == 1) grid.text('A',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 3) grid.text('C',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  
  if(row == 2 & column == 1) grid.text('D',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 2) grid.text('E',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 3) grid.text('F',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  
  if(row==1){
    grid.text('Probability (%)',x=unit(-2.5,'lines'),rot=90)
    grid.text('Incubation (day)',y=unit(-2,'lines'))
  }
  
  if(row==2){
    grid.text('Power (%)',x=unit(-3,'lines'),rot=90)
    grid.text(bquote('Difference in '~omega~' (day)'),y=unit(-2,'lines'))
  }
  
  if(row==2 & column==1) {
    grid.text('Exp growth in ref', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('Exp decline in alt', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  }
  
  if(row==2 & column==2) {
    grid.text('Constant in ref', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('Constant in alt', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  }
  
  if(row==2 & column==3) {
    grid.text('Constant in ref', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('Exp growth in alt', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  }

  
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
  if(row == 2 & column == 1){ colourbar.pathogen(x_bottom_left = 0.05, y_bottom_left = 0.7, y_length = 0.050, x_length = 0.6, palette=col.pts.hue, num=1, row=1, column=1, type='gen')}
  
  
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
  
  grid.text(bquote('|'~mu[s[2]]-mu[s[1]]~'|'),x=x_bottom_left+0.2, y=y_bottom_left_start[1]+2.5*y_length, 
            default.units = 'native', gp = gpar(fontsize = 8))
  
  
}

png('figure/epi_dynamics.png',height=16,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=3)))

paneller(1,1)
paneller(1,2)
paneller(1,3)

paneller(2,1,type='gen')
paneller(2,2,type='gen')
paneller(2,3,type='gen')


popViewport()
popViewport()
dev.off()




paneller=function(row = 1,column=1, type=NULL)
{
  if(row==1){xlm=c(0,25); ylm=c(0,0.2)}
  if(row==2){xlm=c(0,2); ylm=c(0,1)}
  
  
  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # data
  if(column==1){
    # is.na(mean.iso) & is.na(var.iso)
    row.index.ref = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.002 & mean.iso == 4 & var.incub==5 & mean.incub==4 & growth_rate==log(2)/3.5, set.no]
    row.index.alt = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso == 4 & var.incub==5 & mean.incub==5 & growth_rate==-log(2)/3.5, set.no]
    row.index = c(row.index.ref, row.index.alt); space = c(0,0.008*2) #reference, alternate
  }
  if(column==2){
    row.index.ref = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.002 & mean.iso == 4 & var.incub==5 & mean.incub==4 & growth_rate==log(2)/3.5, set.no]
    row.index.alt = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso == 4 & var.incub==5 & mean.incub==5 & is.na(growth_rate), set.no]
    row.index = c(row.index.ref, row.index.alt); space = c(0,0.008*2) #reference, alternate
  }
  if(column==3){
    row.index.ref = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.002 & mean.iso == 4 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
    row.index.alt = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso == 4 & var.incub==5 & mean.incub==5 & growth_rate==log(2)/3.5, set.no]
    row.index = c(row.index.ref, row.index.alt); space = c(0,0.008*2) #reference, alternate
  }
  
  # power 
  if(row==2){
    if(type=='gen') row.pair[,`:=`(diff=diff.mean.gen,
                                   pow=pow.gen)]
    if(type=='si') row.pair[,`:=`(diff=diff.mean.si,
                                  pow=pow.si)]
  }
  
  if(row==2 & column==1) data = row.pair[fig=='d']
  if(row==2 & column==2) data = row.pair[fig=='e']
  if(row==2 & column==3) data = row.pair[fig=='f']
  
  
  
  # set colour and symbol
  if(row==1){
    col.lines = list(c('#a63603', lightup('#fd8d3c',0.5)),
                     c('#810f7c', lightup('#8c96c6',0.5)),
                     c('#006837', lightup('#78c679',0.5)))
  }
  
  if(row==2){
    pch.sample.size = c(16,17)
    col.pts.hue = list(c('#fdd0a2', '#fdae6b', '#fd8d3c', '#f16913', '#d94801', '#a63603', '#7f2704'),
                       c('#bfd3e6', '#9ebcda', '#8c96c6', '#8c6bb1', '#88419d', '#810f7c', '#4d004b'),
                       c('#d9f0a3', '#addd8e', '#78c679', '#41ab5d', '#238443', '#006837', '#004529'))
    
    col.lines = c('#fd8d3c', '#8c96c6', '#78c679')
    
    col.lines=sapply(1:length(col.lines), function(x) lightup(col.lines[x], alpha=0.25))
    val=seq(0,7,length.out = 7)
    
  }
  
  # gridlines
  if(row==1){
    for(y.gl in seq(0,0.2,0.008)){
      grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
    
    for(x.gl in seq(0,25,1)){
      grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
  }
  
  if(row==2){
    for(y.gl in seq(0,1,0.05)){
      grid.lines(xlm, c(y.gl,y.gl), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
    
    for(x.gl in seq(0,2,0.1)){
      grid.lines(c(x.gl,x.gl), ylm, default.units = 'native',gp=gpar(col='#F0F0F0',lwd=0.75))
    }
  }
  
  
  # hist plots 
  if(row==1){
    for(i in 1:2){
      data=hist(gen.dist[row.index[i], ], breaks=c(0:30), plot = F)
      
      for(b in 1:(length(data$breaks)-1)){
        if(data$breaks[b+1]<=25){
          grid.polygon(c(data$breaks[b], data$breaks[b+1], data$breaks[b+1], data$breaks[b]),
                       c(0,0,data$counts[b]/sum(data$counts),data$counts[b]/sum(data$counts)), default.units = 'native', gp=gpar(fill=col.lines[[column]][i], col=col.lines[[column]][i]))
          
        }
      }
      
      mean=round(mean(gen.dist[row.index[i], ]), 1)
      sd=round(sd(gen.dist[row.index[i], ]), 1)
      grid.text(bquote(omega[.(i)] == .(mean) ~ ' (' ~ .(sd) ~')'), x=15, y=0.19-space[i], just='left', default.units = 'native', gp = gpar(fontsize = 8))
      
    }
  }
  
  # power plots 
  if(row==2){
    data[,diff.incub:=abs(mean.incub.alt-mean.incub.ref)]
    size=sort(unique(data$sample.size), decreasing = T)
    
    for(s in 1:length(size)){
      
      data.spline=data[sample.size==size[s]]
      data.points=data[sample.size==size[s] & diff<=1.99 & diff>=0.05 & pow<=0.99]
      
      for(g in 1:3){
        
        # fit spline
        fit=predict(smooth.spline(data.spline[grp==g]$diff, data.spline[grp==g]$pow), seq(0,2,0.01))
        if(column==1) fit=predict(smooth.spline(data.spline[grp==g]$diff, data.spline[grp==g]$pow, nknots=10), seq(0,2,0.01))
        
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
  }
  
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(row==1){
    grid.xaxis(at=seq(0,25,5),label=seq(0,25,5),gp=gpar(fontsize=unit(7,'pt')))
    grid.yaxis(at=seq(0,0.2,0.05),label=seq(0,20,5),gp=gpar(fontsize=unit(7,'pt')))
  }
  if(row==2){
    grid.xaxis(at=seq(0,2,0.5),label=seq(0,2,0.5),gp=gpar(fontsize=unit(7,'pt')))
    grid.yaxis(at=seq(0,1,0.2),label=seq(0,100,20),gp=gpar(fontsize=unit(7,'pt')))
  }
  
  
  
  # labels
  if(row == 1 & column == 1) grid.text('A',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 3) grid.text('C',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  
  if(row == 2 & column == 1) grid.text('D',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 2) grid.text('E',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 2 & column == 3) grid.text('F',x=unit(-2.5,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  
  if(row==1){
    grid.text('Probability (%)',x=unit(-2.5,'lines'),rot=90)
    grid.text('Incubation (day)',y=unit(-2,'lines'))
  }
  
  if(row==2){
    grid.text('Power (%)',x=unit(-3,'lines'),rot=90)
    grid.text(bquote('Difference in '~omega~' (day)'),y=unit(-2,'lines'))
  }
  
  if(row==2 & column==1) {
    grid.text('Exp growth in ref', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('Exp decline in alt', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  }
  
  if(row==2 & column==2) {
    grid.text('Constant in ref', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('Constant in alt', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  }
  
  if(row==2 & column==3) {
    grid.text('Constant in ref', x=0.05, y=0.95, default.units = 'native', just='left', gp = gpar(fontsize = 8))
    grid.text('Exp growth in alt', x=0.05, y=0.90, default.units = 'native', just='left', gp = gpar(fontsize = 8))
  }
  
  
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
  if(row == 2 & column == 1){ colourbar.pathogen(x_bottom_left = 0.05, y_bottom_left = 0.7, y_length = 0.050, x_length = 0.6, palette=col.pts.hue, num=2, row=2, column=1, type='gen')}
  
  
  popViewport()
  popViewport()
  
}


png('figure/epi_dynamics_supp.png',height=16,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=3)))

paneller(1,1)
paneller(1,2)
paneller(1,3)

paneller(2,1,type='gen')
paneller(2,2,type='gen')
paneller(2,3,type='gen')


popViewport()
popViewport()
dev.off()
