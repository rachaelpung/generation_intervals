source('code/step0_load_library.R')

load('output/param_pair_epi_dyn/gen.dist.RData')
load('output/param_pair_epi_dyn/si.dist.RData')
load('output/param_pair_epi_dyn/param.pair.RData')

param[,set.no:=seq_len(.N)]
ref=alt=list()

# diff incub, ref growing, alt decline
ref[['a.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4 & growth_rate==log(2)/3.5, set.no]
alt[['a.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & growth_rate==-log(2)/3.5, set.no]

ref[['a.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5 & mean.incub==4 & growth_rate==log(2)/3.5, set.no]
alt[['a.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5 & growth_rate==-log(2)/3.5, set.no]

ref[['a.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & growth_rate==log(2)/3.5, set.no]
alt[['a.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & growth_rate==-log(2)/3.5, set.no]

# diff incub, ref and alt constant
ref[['b.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['b.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & is.na(growth_rate), set.no]

ref[['b.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['b.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5 & is.na(growth_rate), set.no]

ref[['b.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['b.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & is.na(growth_rate), set.no]

# diff incub, ref constant, alt exp growing
ref[['c.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['c.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & is.na(mean.iso) & is.na(var.iso) & var.incub==5 & growth_rate==log(2)/3.5, set.no]

ref[['c.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['c.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==8 & var.iso==15 & var.incub==5 & growth_rate==log(2)/3.5, set.no]

ref[['c.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['c.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & growth_rate==log(2)/3.5, set.no]

# diff incub, peak, duration shedding, ref growing, alt decline
ref[['d.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & growth_rate==log(2)/3.5, set.no]
alt[['d.1']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & growth_rate==-log(2)/3.5, set.no]

ref[['d.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.002 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & growth_rate==log(2)/3.5, set.no]
alt[['d.2']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & growth_rate==-log(2)/3.5, set.no]

ref[['d.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.006 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & growth_rate==log(2)/3.5, set.no]
alt[['d.3']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & growth_rate==-log(2)/3.5, set.no]

# diff incub, peak, duration shedding, ref and alt constant
ref[['e.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['e.1']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & is.na(growth_rate), set.no]

ref[['e.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.002 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['e.2']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & is.na(growth_rate), set.no]

ref[['e.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.006 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['e.3']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & is.na(growth_rate), set.no]

# diff incub, peak, duration shedding, ref growing, alt constant
ref[['f.1']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['f.1']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & growth_rate==log(2)/3.5, set.no]

ref[['f.2']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.002 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['f.2']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & growth_rate==log(2)/3.5, set.no]

ref[['f.3']] = param[inf.func.name=='spline.covid.delta' & ct.list.name=='hh' & beta==0.006 & mean.iso==4 & var.iso==5 & var.incub==5 & mean.incub==4 & is.na(growth_rate), set.no]
alt[['f.3']] = param[inf.func.name=='spline.covid.wild' & ct.list.name=='hh' & beta==0.0005 & mean.iso==4 & var.iso==5 & var.incub==5 & growth_rate==log(2)/3.5, set.no]

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
    
    pow.gen=c(pow.gen, power.welch.t.test(n=size[s], delta=mu.gen.alt-mu.gen.ref, 
                                          sd1=sd(gen.ref), sd2=sd(gen.alt), 
                                          sig.level = 0.05, alternative = 'two.sided')$power)
    
    pow.si=c(pow.si, power.welch.t.test(n=size[s], delta=mu.si.alt-mu.si.ref, 
                                        sd1=sd(si.ref), sd2=sd(si.alt), 
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
               pow.si=pow.si)]
