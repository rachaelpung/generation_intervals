#' setup functions to convert sample mean and variance to distribution parameters
#' @param func.name: function name
#' @param mu: mean
#' @param var: variance
param.convert<-function(func.name, mu, var){
  
  if(func.name=='lnorm' & !is.na(mu) & !is.na(var)){
    param.lnorm(mu,var)
  } else if(func.name=='gamma' & !is.na(mu) & !is.na(var)){
    param.gamma(mu,var)
  } else if(func.name=='weibull' & !is.na(mu) & !is.na(var)){
    param.weibull(mu,var)
  } else if(func.name=='normal' & !is.na(mu) & !is.na(var)){
    param.normal(mu,var)
  } else if(func.name=='uniform' & !is.na(mu) & !is.na(var)){
    c(mu,var)
  }
  
}

param.lnorm<-function(mu,var){
  
  sdlog=sqrt(log((var/(mu^2))+1))
  meanlog=log(mu)-0.5*(sdlog^2)
  
  list(meanlog=meanlog,
       sdlog=sdlog)
}

param.gamma<-function(mu,var){
  
  shape=(mu^2)/var
  scale=var/mu
  
  list(shape=shape,
       scale=scale)
}

param.weibull<-function(mu,var){
  
  shape=uniroot(function(k) log(factorial(2/k))-2*log(factorial(1/k))-log(var+mu^2)+log(mu^2), lower=0.1, upper=20)$root
  scale=mu/factorial(1/shape)
  
  list(shape=shape,
       scale=scale)
}

param.normal<-function(mu,var){
  
  mean=mu
  sd=sqrt(var)
  
  list(mean=mu,
       sd=sd)
}


#' parameter combination for pairwise transmission
param.pair<-function(){
  
  param = list()
  nsim=1000; max.inf.day=30
  
  #sars-cov-2 delta, varying incub
  param[[1]] = data.table(nsim, mean.incub=seq(0.05,10,0.05), var.incub=5, inc.func.name='lnorm',
                          max.inf.day, beta=0.0005, inf.func.name='spline.covid.delta', 
                          iso.func.name='weibull', ct.list.name='hh') 
  param[[1]] = param[[1]][rep(seq_len(nrow(param[[1]])), times=3)]
  param[[1]][,`:=`(mean.iso=rep(c(NA,8,4), each=.N/3),
                var.iso=rep(c(NA,15,5), each=.N/3))]
  
  #sars-cov-2 wild, varying incub
  param[[2]] = data.table(nsim, mean.incub=seq(0.05,15,0.05), var.incub=5, inc.func.name='lnorm',
                          max.inf.day, beta=0.0005, inf.func.name='spline.covid.wild', 
                          iso.func.name='weibull', ct.list.name='hh') 
  param[[2]] = param[[2]][rep(seq_len(nrow(param[[2]])), times=3)]
  param[[2]][,`:=`(mean.iso=rep(c(NA,8,4), each=.N/3),
                   var.iso=rep(c(NA,15,5), each=.N/3))]
  
  #sars-cov-2 delta, 20%, 50%, 80% prob of infection
  param[[3]] = data.table(nsim, mean.incub=4, var.incub=5, inc.func.name='lnorm',
                          max.inf.day, beta=c(0.0005, 0.002, 0.006), inf.func.name='spline.covid.delta', 
                          mean.iso=4, var.iso=5, iso.func.name='weibull', ct.list.name='hh') 
  
  param[[4]] = data.table(nsim, mean.incub=4, var.incub=5, inc.func.name='lnorm',
                          max.inf.day, beta=c(0.0005, 0.002, 0.006), inf.func.name='spline.covid.delta', 
                          mean.iso=NA, var.iso=NA, iso.func.name='none', ct.list.name='hh') 
  
  #sars-cov-2 delta, varying beta, nhh
  param[[5]] = data.table(nsim, mean.incub=4, var.incub=5, inc.func.name='lnorm',
                          max.inf.day, beta=10^(seq(-5,-1,0.01)), inf.func.name='spline.covid.delta', 
                          iso.func.name='weibull') 
  param[[5]] = param[[5]][rep(seq_len(nrow(param[[5]])), times=3)]
  param[[5]][,`:=`(mean.iso=rep(c(NA,8,4), each=.N/3),
                   var.iso=rep(c(NA,15,5), each=.N/3))]
  param[[5]] = param[[5]][rep(seq_len(nrow(param[[5]])), times=2)]
  param[[5]][,`:=`(ct.list.name=rep(c('hh','nhh'), each=.N/2))]
  
  #sars-cov-2 delta, varying beta, weekly
  param[[6]] = copy(param[[5]])
  param[[6]][,`:=`(ct.list.name=rep(c('hh','hh.wk.01'), each=.N/2))]
  
  #sars-cov-2 wild, ferretti
  param[[7]] = data.table(nsim, mean.incub=seq(0.05,10,0.05), var.incub=5, inc.func.name='lnorm',
                          max.inf.day, beta=1, inf.func.name='spline.covid.wild.ferretti', 
                          iso.func.name='weibull', ct.list.name='hh') 
  param[[7]] = param[[7]][rep(seq_len(nrow(param[[7]])), times=3)]
  param[[7]][,`:=`(mean.iso=rep(c(NA,8,4), each=.N/3),
                   var.iso=rep(c(NA,15,5), each=.N/3))]
  
  param=rbindlist(param, use.names = T, fill = T)
  
  param[,inc.pars.1:=param.convert(inc.func.name, mean.incub, var.incub)[1], by = seq_len(nrow(param))]
  param[,inc.pars.2:=param.convert(inc.func.name, mean.incub, var.incub)[2], by = seq_len(nrow(param))]
  param[,iso.pars.1:=param.convert(iso.func.name, mean.iso, var.iso)[1], by = seq_len(nrow(param))]
  param[,iso.pars.2:=param.convert(iso.func.name, mean.iso, var.iso)[2], by = seq_len(nrow(param))]
  
  param[is.na(mean.iso) & is.na(var.iso), iso.func.name:='none']
  param=param[(!is.na(mean.iso) & !is.na(var.iso))| (is.na(mean.iso) & is.na(var.iso))]
  param=unique(param)
  
}

#' parameter combination for pairwise transmission for respective diseases
param.pair.disease<-function(){
  
  param = list()
  nsim=1000; max.inf.day=30
  
  param[[1]] = data.table(nsim, fig='a', ct.list.name='hh', max.inf.day, 
                          inf.func.name=c('spline.covid.wild', 'spline.covid.delta',
                                          'spline.measles', 'spline.smallpox', 'spline.flu'),
                          mean.incub=c(5.519, 4.145, 14, 12, 2),
                          var.incub=c(5.876, 3.8607, 2.25, 1, 0.25),
                          inc.func.name=c('lnorm', 'weibull', 'normal', 'normal', 'normal'),
                          iso.func.name='weibull', index=1:5)

  param[[1]] = param[[1]][rep(1:.N, each=36)]
  param[[1]][,beta:=rep(10^(seq(-5,-1.5,0.1)), times=5)]
 
  param[[1]] = param[[1]][rep(1:.N, times=3)]
  param[[1]][,`:=`(mean.iso=rep(c(NA,8,4), each=.N/3),
                   var.iso=rep(c(NA,15,5), each=.N/3))]
  
  setorder(param[[1]], index)
  param[[1]][,index:=NULL]
  
  param[[2]] = data.table(nsim, fig='b', ct.list.name='hh', max.inf.day, 
                          inf.func.name=c('spline.covid.wild', 'spline.covid.delta'),
                          mean.incub=c(5.519, 4.145),
                          var.incub=c(5.876, 3.8607),
                          inc.func.name=c('lnorm', 'weibull'),
                          beta=c(0.0005011872, 0.0005011872),
                          iso.func.name='weibull')
  
  param[[2]] = param[[2]][rep(1:.N, each=35)]
  param[[2]][,`:=`(mean.iso=rep(seq(0.1,15,length.out=35), times=.N/35),
                   var.iso=5)]
  
  param=rbindlist(param, use.names = T, fill = T)
  
  param[,inc.pars.1:=param.convert(inc.func.name, mean.incub, var.incub)[1], by = seq_len(nrow(param))]
  param[,inc.pars.2:=param.convert(inc.func.name, mean.incub, var.incub)[2], by = seq_len(nrow(param))]
  param[,iso.pars.1:=param.convert(iso.func.name, mean.iso, var.iso)[1], by = seq_len(nrow(param))]
  param[,iso.pars.2:=param.convert(iso.func.name, mean.iso, var.iso)[2], by = seq_len(nrow(param))]
  
  param[is.na(mean.iso) & is.na(var.iso), iso.func.name:='none']
  param=param[(!is.na(mean.iso) & !is.na(var.iso))| (is.na(mean.iso) & is.na(var.iso))]
  param=unique(param)
  
}

#' parameter combination for cluster transmission
param.cluster<-function(){
  
  nsim=1000; max.inf.day=30
  
  param = data.table(nsim,mean.incub=c(2,4,6),var.incub=5, 
                     inc.func.name='lnorm',max.inf.day,
                     beta=0.0005,inf.func.name='spline.covid.delta',
                     iso.func.name='weibull',ct.list.name='hh',
                     hh.size=3)
                     
  param=param[rep(seq_len(nrow(param)), times=91)]
  param[,`:=`(mean.iso=rep(seq(1,10,0.1), each=.N/91),
                var.iso=5)]
  setorder(param, mean.incub)
  
  param[,inc.pars.1:=param.convert(inc.func.name, mean.incub, var.incub)[1], by = seq_len(nrow(param))]
  param[,inc.pars.2:=param.convert(inc.func.name, mean.incub, var.incub)[2], by = seq_len(nrow(param))]
  param[,iso.pars.1:=param.convert(iso.func.name, mean.iso, var.iso)[1], by = seq_len(nrow(param))]
  param[,iso.pars.2:=param.convert(iso.func.name, mean.iso, var.iso)[2], by = seq_len(nrow(param))]
  
  param[is.na(mean.iso) & is.na(var.iso), iso.func.name:='none']
  param=param[(!is.na(mean.iso) & !is.na(var.iso))| (is.na(mean.iso) & is.na(var.iso))]
  param=unique(param)
  param[,set.no:=1:.N]
}

#' parameter combination for cluster transmission
param.cluster.schematic<-function(){
  
  nsim=1000; max.inf.day=30
  
  param = data.table(nsim,mean.incub=4,var.incub=5, 
                     inc.func.name='lnorm',max.inf.day,
                     beta=0.0005,inf.func.name='spline.covid.delta',
                     iso.func.name='weibull',ct.list.name='hh',
                     hh.size=3)
  
  param=param[rep(seq_len(nrow(param)), times=19)]
  param[,`:=`(mean.iso=rep(seq(1,10,0.5), each=.N/19),
              var.iso=5)]
  setorder(param, mean.incub)
  
  param[,inc.pars.1:=param.convert(inc.func.name, mean.incub, var.incub)[1], by = seq_len(nrow(param))]
  param[,inc.pars.2:=param.convert(inc.func.name, mean.incub, var.incub)[2], by = seq_len(nrow(param))]
  param[,iso.pars.1:=param.convert(iso.func.name, mean.iso, var.iso)[1], by = seq_len(nrow(param))]
  param[,iso.pars.2:=param.convert(iso.func.name, mean.iso, var.iso)[2], by = seq_len(nrow(param))]
  
  param[is.na(mean.iso) & is.na(var.iso), iso.func.name:='none']
  param=param[(!is.na(mean.iso) & !is.na(var.iso))| (is.na(mean.iso) & is.na(var.iso))]
  param=unique(param)
  param[,set.no:=1:.N]
}

# parameter sensitivity analysis epidemic dynamics
# param.pair.epi.dyn<-function(){
# 
#   param = list()
#   nsim=1000; max.inf.day=30
# 
#   #sars-cov-2 delta, varying incub
#   param[[1]] = data.table(nsim, mean.incub=seq(0.05,3.95,0.05), var.incub=5, inc.func.name='lnorm',
#                           max.inf.day, beta=0.0005, inf.func.name='spline.covid.delta',
#                           iso.func.name='weibull', ct.list.name='hh')
#   param[[1]] = param[[1]][rep(seq_len(nrow(param[[1]])), times=3)]
#   param[[1]][,`:=`(mean.iso=rep(c(NA,8,4), each=.N/3),
#                    var.iso=rep(c(NA,15,5), each=.N/3))]
#   param[[1]] = param[[1]][rep(seq_len(nrow(param[[1]])), times=3)]
#   param[[1]][, growth_rate:=rep(c(-log(2)/3.5, NA, log(2)/3.5), each=.N/3)]
# 
#   #sars-cov-2 wild, varying incub
#   param[[2]] = data.table(nsim, mean.incub=seq(0.05,3.95,0.05), var.incub=5, inc.func.name='lnorm',
#                           max.inf.day, beta=0.0005, inf.func.name='spline.covid.wild',
#                           iso.func.name='weibull', ct.list.name='hh')
#   param[[2]] = param[[2]][rep(seq_len(nrow(param[[2]])), times=3)]
#   param[[2]][,`:=`(mean.iso=rep(c(NA,8,4), each=.N/3),
#                    var.iso=rep(c(NA,15,5), each=.N/3))]
#   param[[2]] = param[[2]][rep(seq_len(nrow(param[[2]])), times=3)]
#   param[[2]][, growth_rate:=rep(c(-log(2)/3.5, NA, log(2)/3.5), each=.N/3)]
# 
# 
#   # #sars-cov-2 delta, 20%, 50%, 80% prob of infection
#   # param[[3]] = data.table(nsim, mean.incub=4, var.incub=5, inc.func.name='lnorm',
#   #                         max.inf.day, beta=c(0.0005, 0.002, 0.006), inf.func.name='spline.covid.delta',
#   #                         mean.iso=4, var.iso=5, iso.func.name='weibull', ct.list.name='hh')
#   # param[[3]] = param[[3]][rep(seq_len(nrow(param[[3]])), times=3)]
#   # param[[3]][, growth_rate:=rep(c(-log(2)/3.5, NA, log(2)/3.5), each=.N/3)]
#   #
#   # param[[4]] = data.table(nsim, mean.incub=4, var.incub=5, inc.func.name='lnorm',
#   #                         max.inf.day, beta=c(0.0005, 0.002, 0.006), inf.func.name='spline.covid.delta',
#   #                         mean.iso=NA, var.iso=NA, iso.func.name='none', ct.list.name='hh')
#   # param[[4]] = param[[4]][rep(seq_len(nrow(param[[4]])), times=3)]
#   # param[[4]][, growth_rate:=rep(c(-log(2)/3.5, NA, log(2)/3.5), each=.N/3)]
#   #
# 
# 
#   param=rbindlist(param, use.names = T, fill = T)
# 
#   param[,inc.pars.1:=param.convert(inc.func.name, mean.incub, var.incub)[1], by = seq_len(nrow(param))]
#   param[,inc.pars.2:=param.convert(inc.func.name, mean.incub, var.incub)[2], by = seq_len(nrow(param))]
#   param[,iso.pars.1:=param.convert(iso.func.name, mean.iso, var.iso)[1], by = seq_len(nrow(param))]
#   param[,iso.pars.2:=param.convert(iso.func.name, mean.iso, var.iso)[2], by = seq_len(nrow(param))]
# 
#   param[is.na(mean.iso) & is.na(var.iso), iso.func.name:='none']
#   param=param[(!is.na(mean.iso) & !is.na(var.iso))| (is.na(mean.iso) & is.na(var.iso))]
#   param=unique(param)
# 
# }


param.pair.epi.dyn.unadj<-function(){
  
  param = list()
  nsim=1000; max.inf.day=30
  
  #sars-cov-2 delta, varying incub
  param[[1]] = data.table(nsim, mean.incub=4, var.incub=5, inc.func.name='lnorm',
                          max.inf.day, beta=0.0005, inf.func.name='spline.covid.delta', 
                          iso.func.name='weibull', ct.list.name='hh', mean.iso=4, var.iso=5) 
  param[[1]] = param[[1]][rep(seq_len(nrow(param[[1]])), times=11)]
  param[[1]][, growth_rate:=rep(seq(-0.5,0.5,0.1), each=.N/11)]

  #sars-cov-2 wild, varying incub
  param[[2]] = data.table(nsim, mean.incub=5, var.incub=5, inc.func.name='lnorm',
                          max.inf.day, beta=0.0005, inf.func.name='spline.covid.wild', 
                          iso.func.name='weibull', ct.list.name='hh', mean.iso=4, var.iso=5, 
                          growth_rate=c(-0.2,0,0.2)) 
  
  param=rbindlist(param, use.names = T, fill = T)
  
  param[,inc.pars.1:=param.convert(inc.func.name, mean.incub, var.incub)[1], by = seq_len(nrow(param))]
  param[,inc.pars.2:=param.convert(inc.func.name, mean.incub, var.incub)[2], by = seq_len(nrow(param))]
  param[,iso.pars.1:=param.convert(iso.func.name, mean.iso, var.iso)[1], by = seq_len(nrow(param))]
  param[,iso.pars.2:=param.convert(iso.func.name, mean.iso, var.iso)[2], by = seq_len(nrow(param))]
  
  param[is.na(mean.iso) & is.na(var.iso), iso.func.name:='none']
  param=param[(!is.na(mean.iso) & !is.na(var.iso))| (is.na(mean.iso) & is.na(var.iso))]
  param=unique(param)
  
}

# parameter sensitivity analysis epidemic dynamics
param.pair.epi.dyn.adj<-function(){

  # param = list()
  # max.inf.day=30
  # 
  # #sars-cov-2 delta, varying incub
  # param[[1]] = data.table(nsim=c(1000,100,25), mean.incub=4, var.incub=5, inc.func.name='lnorm',
  #                         max.inf.day, beta=0.0005, inf.func.name='spline.covid.delta',
  #                         iso.func.name='weibull', ct.list.name='hh', mean.iso=4, var.iso=5)
  # 
  # param[[1]] = param[[1]][rep(seq_len(nrow(param[[1]])), times=3)]
  # param[[1]][, growth_rate:=rep(c(-log(2)/3.5, NA, log(2)/3.5), each=.N/3)]
  # 
  # #sars-cov-2 wild, varying incub
  # param[[2]] = data.table(mean.incub=seq(0.05,10,0.05), var.incub=5, inc.func.name='lnorm',
  #                         max.inf.day, beta=0.0005, inf.func.name='spline.covid.wild',
  #                         iso.func.name='weibull', ct.list.name='hh', mean.iso=4, var.iso=5)
  # param[[2]] = param[[2]][rep(seq_len(nrow(param[[2]])), times=3)]
  # param[[2]][, nsim:=rep(c(1000,100,25), each=.N/3)]
  # 
  # param[[2]] = param[[2]][rep(seq_len(nrow(param[[2]])), times=3)]
  # param[[2]][, growth_rate:=rep(c(-log(2)/3.5, NA, log(2)/3.5), each=.N/3)]
  # 
  # 
  # param=rbindlist(param, use.names = T, fill = T)
  # 
  # param[,inc.pars.1:=param.convert(inc.func.name, mean.incub, var.incub)[1], by = seq_len(nrow(param))]
  # param[,inc.pars.2:=param.convert(inc.func.name, mean.incub, var.incub)[2], by = seq_len(nrow(param))]
  # param[,iso.pars.1:=param.convert(iso.func.name, mean.iso, var.iso)[1], by = seq_len(nrow(param))]
  # param[,iso.pars.2:=param.convert(iso.func.name, mean.iso, var.iso)[2], by = seq_len(nrow(param))]
  # 
  # param[is.na(mean.iso) & is.na(var.iso), iso.func.name:='none']
  # param=param[(!is.na(mean.iso) & !is.na(var.iso))| (is.na(mean.iso) & is.na(var.iso))]
  # param=unique(param)
  
  param = param.pair.epi.dyn.unadj()
  param = param[rep(seq_len(nrow(param)), times=300)]
  param[, nsim:=rep(c(1000,100,25), each=.N/3)]

}
