source('code/step0_load_library.R')
source('code/step2b_inf_pair_param.R')
source('code/step2a_inf_pair_function.R')
source('code/step3x_inf_cluster_function.R')

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


# set parameters
param = param.cluster()


set.seed(123) 
## set up clusters
cl = makeCluster(detectCores())
clusterExport(cl, as.vector(lsf.str()))
registerDoParallel(cl)

foreach(p = 1:param[,.N], .packages = c('dplyr', 'purrr', 'data.table', 'glogis', 'EnvStats'), .combine = 'c') %dopar%  {
  
  param.set = param[p,]
  
  sim.output <- gi.cluster(set.no=param.set$set.no, nsim=param.set$nsim,
                           net=ct.list[[param.set$ct.list.name]], 
                           mean.incub=param.set$mean.incub, var.incub=param.set$var.incub,
                           inc.func.name=param.set$inc.func.name,
                           max.inf.day=param.set$max.inf.day, beta=param.set$beta, 
                           inf.func.name=param.set$inf.func.name,
                           mean.iso=param.set$mean.iso, var.iso=param.set$var.iso,
                           iso.func.name=param.set$iso.func.name,
                           hh.size=param.set$hh.size)
  
  if(p<10){
    save(sim.output,file= paste0("output/20220831/sim_output_00", p, ".rdata"))
  } else if(p<100){
    save(sim.output,file= paste0("output/20220831/sim_output_0", p, ".rdata"))
  } else{
    save(sim.output,file= paste0("output/20220831/sim_output_", p, ".rdata"))
  }
  
  
}


# stop clusters
stopCluster(cl) 

# combine cluster outputs
list.folder = dir('output/20220831/cluster output', pattern = '')
list.folder
results=list()

for(f in 1:length(list.folder)){
  
  output = load(paste('output/20220815/cluster output/', list.folder[f], sep = ''))
  results[[f]] = get(output)
  
  
}

results = rbindlist(results)

setorder(results, set.no, sim)
save(results, file='output/20220831/cluster_20220831_w.RData')

