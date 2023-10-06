# load outputs
load('output/param_cluster/gen.dist.RData')
load('output/param_cluster/si.dist.RData')
load('output/param_cluster/param.pair.RData')
load('output/param_cluster/inc.dist.RData')
load('output/param_cluster/cluster.RData')
load('output/param_cluster/tost.dist.RData')


# reformat results
setnames(results, old=c('id', 'exposure', 'onset'), new=c('infectee', 'infectee.exposure', 'infectee.onset'))
results[results, infector.exposure:=i.infectee.exposure, on=c(infector='infectee', sim='sim', set.no='set.no')]
results[results, infector.onset:=i.infectee.onset, on=c(infector='infectee', sim='sim', set.no='set.no')]
results[,gen.int:=infectee.exposure-infector.exposure]
results[,si.int:=infectee.onset-infector.onset]


row.cluster=data.table(set.no=1:273)
row.cluster[param, mean.incub:=i.mean.incub, on=c(set.no='set.no')]
row.cluster[param, mean.iso:=i.mean.iso, on=c(set.no='set.no')]
row.cluster[param, beta:=i.beta, on=c(set.no='set.no')]

# power calculation
# diff.mean.gen=diff.mean.si=as.numeric()
mean.gen.x=mean.gen.y=mean.si.x=mean.si.y=var.gen.x=var.gen.y=var.si.x=var.si.y=as.numeric()
pow.gen=pow.si=as.numeric()
size=c(25,100) #c(25, 50, 75, 100)

for(r in 1:row.cluster[,.N]){
  
  set.seed(123)
  
  set=row.cluster[r,set.no]
  gen.x=results[set.no%in%set & generation !=1 & !is.na(gen.int),gen.int]
  gen.y=gen.dist[set,]
  gen.y=gen.y[!is.na(gen.y)]
  
  gen.x=gen.x[1:1000]
  gen.y=gen.y[1:1000]
  
  si.x=results[set.no%in%set & generation !=1 & !is.na(si.int),si.int]
  si.y=si.dist[set,]
  si.y=si.y[!is.na(si.y)]
  
  si.x=si.x[1:1000]
  si.y=si.y[1:1000]
  
  
  mu.gen.x=mean(gen.x)
  mu.gen.y=mean(gen.y)
  sd.gen.x=sd(gen.x)
  sd.gen.y=sd(gen.y)
  
  mu.si.x=mean(si.x)
  mu.si.y=mean(si.y)
  sd.si.x=sd(si.x)
  sd.si.y=sd(si.y)
  
  mean.gen.x=c(mean.gen.x, mu.gen.x)
  mean.gen.y=c(mean.gen.y, mu.gen.y)
  mean.si.x=c(mean.si.x, mu.si.x)
  mean.si.y=c(mean.si.y, mu.si.y)
  
  var.gen.x=c(var.gen.x,sd.gen.x)
  var.gen.y=c(var.gen.y,sd.gen.y)
  var.si.x=c(var.si.x,sd.si.x)
  var.si.y=c(var.si.y,sd.si.y)
  
  for(s in 1:length(size)){
    
    pow.gen=c(pow.gen, power.welch.t.test(n=size[s], delta=mu.gen.y-mu.gen.x,
                                          sd1=sd(gen.x), sd2=sd(gen.y),
                                          sig.level = 0.05, alternative = 'two.sided')$power)
    
    pow.si=c(pow.si, power.welch.t.test(n=size[s], delta=mu.si.y-mu.si.x,
                                        sd1=sd(si.x), sd2=sd(si.y),
                                        sig.level = 0.05, alternative = 'two.sided')$power)
    
  }
  
}


# tabulate
row.cluster[,`:=`(mean.gen.x=mean.gen.x,
                  mean.gen.y=mean.gen.y,
                  mean.si.x=mean.si.x,
                  mean.si.y=mean.si.y,
                  var.gen.x=var.gen.x^2,
                  var.gen.y=var.gen.y^2,
                  var.si.x=var.si.x^2,
                  var.si.y=var.si.y^2)]

row.cluster[,`:=`(diff.mean.gen=mean.gen.y-mean.gen.x,
                  diff.mean.si=mean.si.y-mean.si.x)]
n=1000
row.cluster[,diff.mean.gen.ci.low:=diff.mean.gen-1.96*sqrt((var.gen.x/n)+(var.gen.y/n))]
row.cluster[,diff.mean.gen.ci.upp:=diff.mean.gen+1.96*sqrt((var.gen.x/n)+(var.gen.y/n))]
row.cluster[,diff.mean.si.ci.low:=diff.mean.si-1.96*sqrt((var.si.x/n)+(var.si.y/n))]
row.cluster[,diff.mean.si.ci.upp:=diff.mean.si+1.96*sqrt((var.si.x/n)+(var.si.y/n))]

row.cluster=row.cluster[rep(seq_len(nrow(row.cluster)), each=length(size)),]
row.cluster[,sample.size:=rep(size, times=.N/length(size))]
row.cluster[,`:=`(pow.gen=pow.gen,
                  pow.si=pow.si)]


# analysis for cluster and pair wise difference
param[mean.incub==4 & mean.iso==4,set.no]
summary(tost.dist[122,]); quantile(tost.dist[122,], probs = c(0.025,0.975))
tost.cluster = results[set.no==122]$gen.int - (results[set.no==122]$infector.onset-results[set.no==122]$infector.exposure)
tost.cluster = tost.cluster[which(!is.na(tost.cluster))]
summary(tost.cluster); quantile(tost.cluster, probs = c(0.025,0.975))


param[mean.incub==4 & mean.iso==8,set.no]
summary(tost.dist[162,]); quantile(tost.dist[162,], probs = c(0.025,0.975))
tost.cluster = results[set.no==162]$gen.int - (results[set.no==162]$infector.onset-results[set.no==162]$infector.exposure)
tost.cluster = tost.cluster[which(!is.na(tost.cluster))]
summary(tost.cluster); quantile(tost.cluster, probs = c(0.025,0.975))


View(results)
setorder(results, set.no, sim, infectee.exposure)
results[,row.num:=1:.N]

N=results[,.N]
set=param[mean.incub ==4 & mean.iso==4,set.no]

row.ter.case = seq(3,N,3) 
row.ter.case = results[row.num %in% row.ter.case & generation==3, row.num]
inc.ter.case = results[row.num %in% row.ter.case & set.no%in%set]$infector.onset-results[row.num %in% row.ter.case & set.no%in%set]$infector.exposure
summary(inc.ter.case)


inc.sec.case = results[row.num %in% (row.ter.case-1) & set.no%in%set]$infector.onset-results[row.num %in% (row.ter.case-1) & set.no%in%set]$infector.exposure
summary(inc.sec.case)

t.test(inc.sec.case,inc.ter.case)
summary(inc.ter.case-inc.sec.case)
