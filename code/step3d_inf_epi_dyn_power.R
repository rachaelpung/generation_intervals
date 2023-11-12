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

row.pair.unadj = copy(row.pair)
param.unadj = copy(param)
rm(row.pair, param)



# un/adjusted power
load('output/param_pair_epi_unadj_pow/gen.dist.RData')
load('output/param_pair_epi_unadj_pow/si.dist.RData')
load('output/param_pair_epi_unadj_pow/param.pair.RData')

param[,set.no:=seq_len(.N)]
ref=alt=list()

# diff incub
ref[['a.1']] = param[inf.func.name=='spline.covid.delta' & growth_rate==0.2, set.no]
alt[['a.1']] = param[inf.func.name=='spline.covid.wild' & growth_rate==-0.2, set.no]

ref[['a.2']] = param[inf.func.name=='spline.covid.delta' & growth_rate==0, set.no]
alt[['a.2']] = param[inf.func.name=='spline.covid.wild' & growth_rate==0, set.no]

ref[['a.3']] = param[inf.func.name=='spline.covid.delta' & growth_rate==0, set.no]
alt[['a.3']] = param[inf.func.name=='spline.covid.wild' & growth_rate==0.2, set.no]

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


# power calculation
diff.mean.gen=diff.mean.si=as.numeric()
pow.gen=pow.si=as.numeric()
# size=c(100)
size=c(25, 100) #c(25, 50, 75, 100)

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

  
  for(s in 1:length(size)){
    
    # power for theoretical GI
    pow.gen=c(pow.gen, power.welch.t.test(n=size[s], delta=mu.gen.alt-mu.gen.ref, 
                                          sd1=sd.gen.ref, sd2=sd.gen.alt, 
                                          sig.level = 0.05, alternative = 'two.sided')$power)
    
    # power for observed SI
    pow.si=c(pow.si, power.welch.t.test(n=size[s], delta=mu.si.alt-mu.si.ref, 
                                        sd1=sd.si.ref, sd2=sd.si.alt, 
                                        sig.level = 0.05, alternative = 'two.sided')$power)
    
  }
  
}


# tabulate
row.pair=row.pair[rep(seq_len(nrow(row.pair)), each=length(size)),]
row.pair[,sample.size:=rep(size, times=.N/length(size))]

row.pair[param,beta.ref:=i.beta, on=c(ref='set.no')]
row.pair[param,mean.incub.ref:=i.mean.incub, on=c(ref='set.no')]
row.pair[param,var.incub.ref:=i.var.incub, on=c(ref='set.no')]
row.pair[param,growth.rate.ref:=i.growth_rate, on=c(ref='set.no')]

row.pair[param,beta.alt:=i.beta, on=c(alt='set.no')]
row.pair[param,mean.incub.alt:=i.mean.incub, on=c(alt='set.no')]
row.pair[param,var.incub.alt:=i.var.incub, on=c(alt='set.no')]
row.pair[param,growth.rate.alt:=i.growth_rate, on=c(alt='set.no')]


row.pair[,`:=`(diff.mean.gen=rep(diff.mean.gen, each=length(size)),
               diff.mean.si=rep(diff.mean.si, each=length(size)),
               pow.gen=pow.gen,
               pow.si=pow.si)]

row.pair.unadj.pow = copy(row.pair)
row.pair.adj.pow = copy(row.pair)

