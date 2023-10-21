source('code/step0_load_library.R')

# unadjusted
load('output/param_pair_epi_unadj/gen.dist.RData')
load('output/param_pair_epi_unadj/si.dist.RData')
load('output/param_pair_epi_unadj/param.pair.RData')

param[,set.no:=seq_len(.N)]
ref=alt=list()

# unadjusted
# diff incub, ref growing, alt decline
ref[['a.1']] = param[mean.incub==4, set.no]
alt[['a.1']] = param[mean.incub==5 & growth_rate=='-0.2', set.no]

ref[['a.2']] = param[mean.incub==4, set.no]
alt[['a.2']] = param[mean.incub==5 & growth_rate==0, set.no]

ref[['a.3']] = param[mean.incub==4, set.no]
alt[['a.3']] = param[mean.incub==5 & growth_rate=='0.2', set.no]

# adjusted
# growing ref and declining alt
ref[['a.1']] = param[nsim==1000 & inf.func.name=='spline.covid.delta' , set.no]
alt[['a.1']] = param[nsim==1000 & inf.func.name=='spline.covid.wild' & growth_rate=='-0.2', set.no]

ref[['a.2']] = param[nsim==100 & inf.func.name=='spline.covid.delta', set.no]
alt[['a.2']] = param[nsim==100 & inf.func.name=='spline.covid.wild' & growth_rate=='-0.2', set.no]

ref[['a.3']] = param[nsim==25 & inf.func.name=='spline.covid.delta', set.no]
alt[['a.3']] = param[nsim==25 & inf.func.name=='spline.covid.wild' & growth_rate=='-0.2', set.no]

# constant ref and alt
ref[['b.1']] = param[nsim==1000 & inf.func.name=='spline.covid.delta', set.no]
alt[['b.1']] = param[nsim==1000 & inf.func.name=='spline.covid.wild' & growth_rate==0, set.no]

ref[['b.2']] = param[nsim==100 & inf.func.name=='spline.covid.delta', set.no]
alt[['b.2']] = param[nsim==100 & inf.func.name=='spline.covid.wild' & growth_rate==0, set.no]

ref[['b.3']] = param[nsim==25 & inf.func.name=='spline.covid.delta', set.no]
alt[['b.3']] = param[nsim==25 & inf.func.name=='spline.covid.wild' & growth_rate==0, set.no]

# declining ref and growing alt
ref[['c.1']] = param[nsim==1000 & inf.func.name=='spline.covid.delta', set.no]
alt[['c.1']] = param[nsim==1000 & inf.func.name=='spline.covid.wild' & growth_rate=='0.2', set.no]

ref[['c.2']] = param[nsim==100 & inf.func.name=='spline.covid.delta', set.no]
alt[['c.2']] = param[nsim==100 & inf.func.name=='spline.covid.wild' & growth_rate=='0.2', set.no]

ref[['c.3']] = param[nsim==25 & inf.func.name=='spline.covid.delta', set.no]
alt[['c.3']] = param[nsim==25 & inf.func.name=='spline.covid.wild' & growth_rate=='0.2', set.no]


row.pair=data.table()
for(i in 1:length(ref)){
  
  if(length(ref[[i]])==1){
    row.pair=rbind(row.pair, data.table(ref=rep(ref[[i]], length(alt[[i]])), alt=alt[[i]], 
                                        fig=strsplit(names(ref)[i], '[.]')[[1]][1],
                                        grp=as.numeric(strsplit(names(ref)[i], '[.]')[[1]][2])))
  }else{
    row.pair=rbind(row.pair, data.table(ref=ref[[i]], alt=alt[[i]],
                                        fig=strsplit(names(ref)[i], '[.]')[[1]][1],
                                        grp=as.numeric(strsplit(names(ref)[i], '[.]')[[1]][2])))
  }
}

row.pair[param,nsim:=i.nsim, on=c(ref='set.no')]


# power calculation
diff.mean.gen=diff.mean.si=as.numeric()
pow.gen=pow.si=as.numeric()

for(r in 1:row.pair[,.N]){
  
  set.seed(123)
  gen.ref=gen.dist[row.pair[r,ref],]
  gen.alt=gen.dist[row.pair[r,alt],]
  
  gen.ref=gen.ref[!is.na(gen.ref)]
  gen.alt=gen.alt[!is.na(gen.alt)]
  
  si.ref=si.dist[row.pair[r,ref],]
  si.alt=si.dist[row.pair[r,alt],]
  
  si.ref=si.ref[!is.na(si.ref)]
  si.alt=si.alt[!is.na(si.alt)]
  
  mu.gen.ref=mean(gen.ref)
  mu.gen.alt=mean(gen.alt)
  sd.gen.ref=sd(gen.ref)
  sd.gen.alt=sd(gen.alt)
  
  mu.si.ref=mean(si.ref)
  mu.si.alt=mean(si.alt)
  sd.si.ref=sd(si.ref)
  sd.si.alt=sd(si.alt)
  
  diff.mean.gen=c(diff.mean.gen, mu.gen.alt-mu.gen.ref)
  diff.mean.si=c(diff.mean.si, mu.si.alt-mu.si.ref)
  
  pow.gen=c(pow.gen, power.welch.t.test(n=row.pair[r,nsim], delta=mu.gen.alt-mu.gen.ref, 
                                          sd1=sd(gen.ref), sd2=sd(gen.alt), 
                                          sig.level = 0.05, alternative = 'two.sided')$power)
    
  pow.si=c(pow.si, power.welch.t.test(n=row.pair[r,nsim], delta=mu.si.alt-mu.si.ref, 
                                        sd1=sd(si.ref), sd2=sd(si.alt), 
                                        sig.level = 0.05, alternative = 'two.sided')$power)
    
  
}


# tabulate
row.pair[param,beta.ref:=i.beta, on=c(ref='set.no')]
row.pair[param,mean.incub.ref:=i.mean.incub, on=c(ref='set.no')]
row.pair[param,growth.ref:=i.growth_rate, on=c(ref='set.no')]

row.pair[param,beta.alt:=i.beta, on=c(alt='set.no')]
row.pair[param,mean.incub.alt:=i.mean.incub, on=c(alt='set.no')]
row.pair[param,growth.alt:=i.growth_rate, on=c(alt='set.no')]


row.pair[,`:=`(diff.mean.gen=diff.mean.gen,
               diff.mean.si=diff.mean.si,
               pow.gen=pow.gen,
               pow.si=pow.si)]


# adjusted
load('output/param_pair_epi_adj/mean.gen.RData')
load('output/param_pair_epi_adj/mean.si.RData')
load('output/param_pair_epi_adj/param.pair.RData')


param[,set.no:=rep(rep(1:100, each=14), times=3)]
param[,row_no:=1:.N]

row.pair=data.table()
n_size = c(1000,100,25)

for(n in 1:length(n_size)){
  for(s in 1:100){
    for(g in 1:3){
      ref.i = param[nsim==n_size[n] & set.no==s & inf.func.name=='spline.covid.delta', row_no]
      
      if(g==1){
        alt.i = param[nsim==n_size[n] & set.no==s & inf.func.name=='spline.covid.wild' & growth_rate=='-0.2', row_no]
        fig.i='a'
      }
      if(g==2){
        alt.i = param[nsim==n_size[n] & set.no==s & inf.func.name=='spline.covid.wild' & growth_rate==0, row_no]
        fig.i = 'b'
      }
      if(g==3){
        alt.i = param[nsim==n_size[n] & set.no==s & inf.func.name=='spline.covid.wild' & growth_rate=='0.2', row_no]
        fig.i = 'c'
      }
      
      row.pair=rbind(row.pair, data.table(ref=ref.i, alt=alt.i, fig=fig.i, grp=n, set=s))
                     
    }
  }
}

setorder(row.pair, fig, grp, set)


diff.mean.gen=diff.mean.si=as.numeric()

for(r in 1:row.pair[,.N]){
  
  set.seed(123)
  
  ref.i = row.pair$ref[r]
  alt.i = row.pair$alt[r]
  
  mu.gen.ref=mean.gen[[ref.i]][1]
  mu.gen.alt=mean.gen[[alt.i]][1]
  # sd.gen.ref=mean.gen[[ref.i]][2]
  # sd.gen.alt=mean.gen[[alt.i]][2]
  
  mu.si.ref=mean.si[[ref.i]][1]
  mu.si.alt=mean.si[[alt.i]][1]
  # sd.si.ref=mean.si[[ref.i]][2]
  # sd.si.alt=mean.si[[alt.i]][2]
  
  diff.mean.gen=c(diff.mean.gen, mu.gen.alt-mu.gen.ref)
  diff.mean.si=c(diff.mean.si, mu.si.alt-mu.si.ref)
  
}

# tabulate
row.pair[param,beta.ref:=i.beta, on=c(ref='row_no')]
row.pair[param,mean.incub.ref:=i.mean.incub, on=c(ref='row_no')]
row.pair[param,growth.ref:=i.growth_rate, on=c(ref='row_no')]

row.pair[param,beta.alt:=i.beta, on=c(alt='row_no')]
row.pair[param,mean.incub.alt:=i.mean.incub, on=c(alt='row_no')]
row.pair[param,growth.alt:=i.growth_rate, on=c(alt='row_no')]

row.pair[,`:=`(diff.mean.gen=diff.mean.gen,
               diff.mean.si=diff.mean.si)]

# mu.gen=lapply(mean.gen, `[[`, 1)
# mu.gen=unlist(mu.gen)
# param[, mu.gen:=mu.gen]

row.pair = row.pair[,.(median(diff.mean.gen), quantile(diff.mean.gen, 0.25), quantile(diff.mean.gen, 0.75),
                median(diff.mean.si), quantile(diff.mean.si, 0.25), quantile(diff.mean.si, 0.75)), by=.(fig, grp, growth.ref, growth.alt)]
setnames(row.pair, c('fig', 'grp', 'growth.ref', 'growth.alt',
                     'med.diff.mean.gen', 'lwr.diff.mean.gen', 'upp.diff.mean.gen',
                     'med.diff.mean.gen', 'lwr.diff.mean.gen', 'upp.diff.mean.gen'))

row.pair[grp==1, nsim:=1000]
row.pair[grp==2, nsim:=100]
row.pair[grp==3, nsim:=25]
