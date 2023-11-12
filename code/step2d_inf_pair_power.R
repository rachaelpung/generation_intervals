source('code/step0_load_library.R')

load('output/param_pair/gen.dist.RData')
load('output/param_pair/si.dist.RData')
load('output/param_pair/param.pair.RData')
load('output/param_pair/inc.dist.RData')

# load('output/param_pair/gen.dist_sensitivity_analysis_peak_inf.RData')
# load('output/param_pair/si.dist_sensitivity_analysis_peak_inf.RData')
# load('output/param_pair/param.pair_sensitivity_analysis_peak_inf.RData')

param[,set.no:=seq_len(.N)]
ref=alt=list()

# diff incub
ref[['a.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4, set.no]
alt[['a.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5, set.no]

ref[['a.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5 & mean.incub==4, set.no]
alt[['a.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5, set.no]

ref[['a.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4, set.no]
alt[['a.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5, set.no]

# diff incub, peak, duration shedding
ref[['b.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4, set.no]
alt[['b.1']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5, set.no]

ref[['b.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.002 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4, set.no]
alt[['b.2']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5, set.no]

ref[['b.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.006 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4, set.no]
alt[['b.3']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5, set.no]

# diff incub, duration shedding
ref[['c.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4, set.no]
alt[['c.1']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5, set.no]

ref[['c.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5 & mean.incub==4, set.no]
alt[['c.2']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5, set.no]

ref[['c.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4, set.no]
alt[['c.3']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5, set.no]

# diff contact, hh vs nhh
ref[['d.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]
alt[['d.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='nhh' & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]

ref[['d.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & mean.iso==8 & var.iso==15 & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]
alt[['d.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='nhh' & mean.iso==8 & var.iso==15 & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]

ref[['d.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]
alt[['d.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='nhh' & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]

# diff contact freq, daily vs weekly
ref[['e.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]
alt[['e.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh.wk.01' & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]

ref[['e.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & mean.iso==8 & var.iso==15 & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]
alt[['e.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh.wk.01' & mean.iso==8 & var.iso==15 & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]

ref[['e.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]
alt[['e.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh.wk.01' & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & beta%in%10^(seq(-5,-1,0.01)), set.no]

# ferretti vs spline
ref[['f.1']] = param[inf.func.name=='spline.covid.wild.ferretti' & ct.list.name=='hh' & beta==1 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4, set.no]
alt[['f.1']] = param[inf.func.name=='spline.covid.wild.ferretti' & ct.list.name=='hh' & beta==1 & is.na(mean.iso) & is.na(var.iso) & var.incub==5, set.no]

ref[['f.2']] = param[inf.func.name=='spline.covid.wild.ferretti' & ct.list.name=='hh' & beta==1 & mean.iso==8 & var.iso==15 & var.incub==5 & mean.incub==4, set.no]
alt[['f.2']] = param[inf.func.name=='spline.covid.wild.ferretti' & ct.list.name=='hh' & beta==1 & mean.iso==8 & var.iso==15 & var.incub==5, set.no]

ref[['f.3']] = param[inf.func.name=='spline.covid.wild.ferretti' & ct.list.name=='hh' & beta==1 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4, set.no]
alt[['f.3']] = param[inf.func.name=='spline.covid.wild.ferretti' & ct.list.name=='hh' & beta==1 & mean.iso==4 & var.iso==5 & var.incub==5, set.no]

ref[['g.1']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4, set.no]
alt[['g.1']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5, set.no]

ref[['g.2']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5 & mean.incub==4, set.no]
alt[['g.2']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5, set.no]

ref[['g.3']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4, set.no]
alt[['g.3']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5, set.no]

# diff incub, peak, duration shedding
ref[['h.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4, set.no]
alt[['h.1']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5, set.no]

ref[['h.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.002 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4, set.no]
alt[['h.2']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5, set.no]

ref[['h.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.006 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4, set.no]
alt[['h.3']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5, set.no]


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
pow.gen=pow.si=pow.dgi=as.numeric()
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
  
  inc.ref=inc.dist[row.pair[r,ref],]
  inc.alt=inc.dist[row.pair[r,alt],]
  
  inc.ref=inc.ref[!is.na(inc.ref)]
  inc.alt=inc.alt[!is.na(inc.alt)]
  
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
  
  # sd of derived GI from observed SI
  var.derived.gen.ref=var(si.ref)-2*var(inc.ref)
  var.derived.gen.alt=var(si.alt)-2*var(inc.alt)
  
  if(var.derived.gen.ref<0) var.derived.gen.ref=0
  if(var.derived.gen.alt<0) var.derived.gen.alt=0
  
  sd.derived.gen.ref=sqrt(var.derived.gen.ref)
  sd.derived.gen.alt=sqrt(var.derived.gen.alt)
  
  for(s in 1:length(size)){
    
    # power for theoretical GI
    pow.gen=c(pow.gen, power.welch.t.test(n=size[s], delta=mu.gen.alt-mu.gen.ref, 
                                          sd1=sd.gen.ref, sd2=sd.gen.alt, 
                                          sig.level = 0.05, alternative = 'two.sided')$power)
    
    # power for observed SI
    pow.si=c(pow.si, power.welch.t.test(n=size[s], delta=mu.si.alt-mu.si.ref, 
                                          sd1=sd.si.ref, sd2=sd.si.alt, 
                                          sig.level = 0.05, alternative = 'two.sided')$power)
    
    # power for derived GI
    pow.dgi=c(pow.dgi, power.welch.t.test(n=size[s], delta=mu.si.alt-mu.si.ref, 
                                        sd1=sd.derived.gen.ref, sd2=sd.derived.gen.alt, 
                                        sig.level = 0.05, alternative = 'two.sided')$power)
    
  }
  
}


# tabulate
row.pair=row.pair[rep(seq_len(nrow(row.pair)), each=length(size)),]
row.pair[,sample.size:=rep(size, times=.N/length(size))]

row.pair[param,beta.ref:=i.beta, on=c(ref='set.no')]
row.pair[param,mean.incub.ref:=i.mean.incub, on=c(ref='set.no')]
row.pair[param,var.incub.ref:=i.var.incub, on=c(ref='set.no')]
row.pair[param,mean.iso.ref:=i.mean.iso, on=c(ref='set.no')]
row.pair[param,var.iso.ref:=i.var.iso, on=c(ref='set.no')]
row.pair[param,p.trans.overall.ref:=i.p.trans.overall, on=c(ref='set.no')]

row.pair[param,beta.alt:=i.beta, on=c(alt='set.no')]
row.pair[param,mean.incub.alt:=i.mean.incub, on=c(alt='set.no')]
row.pair[param,var.incub.alt:=i.var.incub, on=c(alt='set.no')]
row.pair[param,mean.iso.alt:=i.mean.iso, on=c(alt='set.no')]
row.pair[param,var.iso.alt:=i.var.iso, on=c(alt='set.no')]
row.pair[param,p.trans.overall.alt:=i.p.trans.overall, on=c(alt='set.no')]


row.pair[,`:=`(diff.mean.gen=rep(diff.mean.gen, each=length(size)),
               diff.mean.si=rep(diff.mean.si, each=length(size)),
               pow.gen=pow.gen,
               pow.si=pow.si,
               pow.dgi=pow.dgi)]
