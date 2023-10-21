source('code/step0_load_library.R')
source('code/step2a_inf_pair_function.R')
source('code/step2b_inf_pair_param.R')

# load contact list
list.file = dir('data', pattern = 'ct.list')
list.file
ct.list=list()

for(f in 1:length(list.file)){
  file = load(paste('data/', list.file[f], sep = ''))
  ct.list[[f]] = get(file)
  list.file[f] = gsub("ct.list.(.*).10m.RData", "\\1", list.file[f])
 
}

names(ct.list) = list.file
rm(f,file,list.file)
rm(list=ls(pattern = 'ct.list.'))


# trial run
gi.pair(nsim=10, dt=5/(24*60), ct.pair = ct.list[['hh']],
        inc.func.name='lnorm', inf.func.name='spline.covid.delta',
        iso.func.name='weibull',
        inc.pars=c(1.62, 0.42), iso.pars=c(1.62, 0.42),
        max.inf.day = 21, beta=0.01)

# setup storage 
nsim=1000
# param = param.pair()
# param.pair.disease()
param = param.pair.epi.dyn.adj()
gen.dist=matrix(NA, nrow=param[,.N], ncol=nsim) 
si.dist=matrix(NA, nrow=param[,.N], ncol=nsim)

set.seed(123)
## set up clusters
cl = makeCluster(detectCores())
clusterExport(cl, as.vector(lsf.str()))
registerDoParallel(cl)

# results = foreach(p = 1:param[,.N], .packages = c('data.table','glogis','EnvStats')) %dopar%  {
foreach(p = 1736:param[,.N], .packages = c('data.table','glogis','EnvStats'), .combine = 'c') %dopar%  {
    
    
  output = gi.pair(nsim=param[p,nsim], dt=5/(24*60), ct.pair = ct.list[[param[p,ct.list.name]]],
                   inc.func.name=param[p,inc.func.name], inf.func.name=param[p,inf.func.name],
                   iso.func.name=param[p,iso.func.name],
                   inc.pars=c(param[p,inc.pars.1], param[p,inc.pars.2]),
                   iso.pars=c(param[p,iso.pars.1], param[p,iso.pars.2]),
                   max.inf.day = param[p,max.inf.day], beta=param[p,beta], growth_rate = param[p,growth_rate])
  
  # # parallel runs
  # p.trans=output$p.trans
  # p.trans.overall=output$p.trans.overall
  # p.presymp.trans=output$p.presymp.trans
  # gen.dist=output$gen
  # si.dist=output$si
  # tost.dist=output$tost
  # 
  # inc.dist=output$inc
  # onset.to.iso.dist=output$onset.to.iso
  # 
  # list(p.trans.overall=p.trans.overall,
  #      p.presymp.trans=p.presymp.trans,
  #      gen.dist=gen.dist,
  #      si.dist=si.dist,
  #      tost.dist=tost.dist,
  #      inc.dist=inc.dist,
  #      onset.to.iso.dist=onset.to.iso.dist,
  #      p.trans=p.trans)
  
  # print(p)
  
  if(p<10){
    save(output,file= paste0("output/20231015/output_000", p, ".rdata"))
  } else if(p<100){
    save(output,file= paste0("output/20231015/output_00", p, ".rdata"))
  } else if(p<1000){
    save(output,file= paste0("output/20231015/output_0", p, ".rdata"))
  } else{
    save(output,file= paste0("output/20231015/output_", p, ".rdata"))
  }
  
}
  
# for parallel runs
p.trans.overall=lapply(results, `[[`, 1)
p.trans.overall=unlist(p.trans.overall)
param[,p.trans.overall:=p.trans.overall]

p.presymp.trans=lapply(results, `[[`, 2)
p.presymp.trans=unlist(p.presymp.trans)
param[,p.presymp.trans:=p.presymp.trans]

gen.dist=lapply(results, `[[`, 3)
gen.dist=unlist(gen.dist)
gen.dist=matrix(gen.dist, ncol=nsim, byrow = T)

si.dist=lapply(results, `[[`, 4)
si.dist=unlist(si.dist)
si.dist=matrix(si.dist, ncol=nsim, byrow = T)

tost.dist=lapply(results, `[[`, 5)
tost.dist=unlist(tost.dist)
tost.dist=matrix(tost.dist, ncol=nsim, byrow = T)

inc.dist=lapply(results, `[[`, 6)
inc.dist=unlist(inc.dist)
inc.dist=matrix(inc.dist, ncol=nsim, byrow = T)

onset.to.iso.dist=lapply(results, `[[`, 7)
onset.to.iso.dist=unlist(onset.to.iso.dist)
onset.to.iso.dist=matrix(onset.to.iso.dist, ncol=nsim, byrow = T)

p.trans=lapply(results, `[[`, 8)
p.trans=unlist(p.trans)
p.trans=matrix(p.trans, ncol=nsim, byrow = T)


# stop clusters
stopCluster(cl) 


save(param, file='output/20231011/param.pair.RData')
save(gen.dist, file='output/20231011/gen.dist.RData')
save(si.dist, file='output/20231011/si.dist.RData')
save(inc.dist, file='output/20231011/inc.dist.RData')
save(onset.to.iso.dist, file='output/20231011/onset.to.iso.dist.RData')
save(p.trans, file='output/20231011/p.trans.RData')
save(tost.dist, file='output/20231011/tost.dist.RData')


index.1000 = which(param$nsim==1000)
index.100 = which(param$nsim==100)
index.25 = which(param$nsim==25)

col_num=8
dist=lapply(results, `[[`, col_num)
dist.1000=unlist(dist[index.1000])
dist.1000=matrix(dist.1000, ncol=1000, byrow = T)

dist.100=unlist(dist[index.100])
dist.100=matrix(dist.100, ncol=100, byrow = T)
dist.100=cbind(dist.100, matrix(NA, ncol=1000-100, nrow=length(index.100)))

dist.25=unlist(dist[index.25])
dist.25=matrix(dist.25, ncol=25, byrow = T)
dist.25=cbind(dist.25, matrix(NA, ncol=1000-25, nrow=length(index.25)))

dist = rbind(dist.1000, dist.100, dist.25)
p.trans = copy(dist)
save(p.trans, file='output/20231015/p.trans.RData')

list.file = dir('output/20231015')
list.file
length(list.file)
results = foreach(f = 1:length(list.file), .packages = c('data.table','glogis','EnvStats')) %dopar%  {
  file = load(paste('output/20231015/', list.file[f], sep = ''))
  test = get(file)
  
  output=copy(test)
  
  # parallel runs
  p.trans=output$p.trans
  p.trans.overall=output$p.trans.overall
  p.presymp.trans=output$p.presymp.trans
  gen.dist=output$gen
  si.dist=output$si
  tost.dist=output$tost

  inc.dist=output$inc
  onset.to.iso.dist=output$onset.to.iso

  list(p.trans.overall=p.trans.overall,
       p.presymp.trans=p.presymp.trans,
       gen.dist=gen.dist,
       si.dist=si.dist,
       tost.dist=tost.dist,
       inc.dist=inc.dist,
       onset.to.iso.dist=onset.to.iso.dist,
       p.trans=p.trans)
  
}

test = lapply(1:4200, function(x){
  
  mean.si = mean(si.dist[x,], na.rm=T)
  sd.si = sd(si.dist[x,], na.rm=T)
  
  return(c(mean.si, sd.si))
})

mean.si=copy(test)

save(mean.si, file='output/20231015/mean.si.RData')
