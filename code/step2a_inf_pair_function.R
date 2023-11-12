#' generation intervals observed from one infector with breaks in contact
#' @param nsim: number of simulations
#' @param dt: size of time interval in days (i.e. 5/(24*60) corresponds to 5 mins intervals)
#' @param ct.pair: contact pairs
#' @param inc.func.name: name of incubation period function
#' @param inf.func.name: name of infectious period function
#' @param inc.pars: mean and sd log normal incubation period
#' @param iso.pars: mean and sd log normal onset to isolation period
#' @param max.inf.day: max day of infectiousness
#' @param beta: to scale the infectiousness
#' 
gi.pair <- function(nsim=10, dt=5/(24*60), ct.pair = NULL,
                    inc.func.name='lnorm',
                    inf.func.name='spline',
                    iso.func.name='weibull',
                    inc.pars=c(1.62, 0.42), 
                    iso.pars=c(1.62, 0.42),
                    growth_rate=NULL,
                    max.inf.day = 30, beta=1) {
  # 5/(24*60)
  
  # set.seed(seed)
  
  # number of intervals in a day
  num.int = 1/dt

  # incubation period 
  if(is.null(growth_rate) | is.na(growth_rate) | growth_rate==0){
    inc.f = inc.func[[inc.func.name]]
    inc = inc.f(nsim, inc.pars)
  }else{
    
    # sensitivity analysis
    # sample backward incubation period
    t=seq(0,30,0.0001) #0.000001
    inc.f = dlnorm(t, meanlog = inc.pars[1], sdlog = inc.pars[2])
    inc.f = inc.f/sum(inc.f)
    inc.b = exp(-growth_rate*t)*inc.f/sum(exp(-growth_rate*t)*inc.f)
    # inc.b = sample(t, nsim, replace=TRUE, prob=inc.b)
    inc = sample(t, nsim, replace=TRUE, prob=inc.b)
    
    # # sensitivity analysis
    # # adjust incubation period
    # inc.b = data.table(table(inc.b))
    # setnames(inc.b, c('t','N'))
    # inc.b[, P:=N/sum(N)]
    # inc.b[, t:=as.numeric(t)]
    # inc.f = exp(growth_rate*inc.b$t)*inc.b$P/sum(exp(growth_rate*inc.b$t)*inc.b$P)
    # inc = sample(inc.b$t, nsim, replace=T, prob = inc.f)
    
  }
  

  # shift in viral peak relative to symptoms onset
  # vp = runif(nsim, min=1, max=5)
  vp = rep(0,nsim)
  
  # infectiousness period
  inf.f = inf.func[[inf.func.name]]
  inf.matrix = sapply(1:nsim, function(x){ inf.f(inc[x]+vp[x], max.inf.day, dt) })
  # 8640 (intervals) x 10 (nsim) matrix 
 
  # contact pattern
  ct.int.start = gen.int.start = sample(85:276, nsim, replace=T) 
  
  ct.int.end = (ct.int.start+(24*60/5*max.inf.day-1))
  ct.col = sample(1:ncol(ct.pair), nsim, replace=T)
  ct.matrix = sapply(1:nsim, function(x){
    ct.pair[ct.int.start[x]:ct.int.end[x],ct.col[x]]
  })
  # 8640 (intervals) x 10 (nsim) matrix
  
  # start interval for infectiousness and symptoms
  inc.int = ceiling(inc/dt)
  si.int.start = inc.int + gen.int.start
  
  # time from symptoms onset to isolation
  onset.to.iso.f = onset.to.iso.func[[iso.func.name]]
  onset.to.iso = onset.to.iso.f(nsim, iso.pars)
  onset.to.iso.int = ceiling(onset.to.iso)/dt
  iso.int = inc.int + onset.to.iso.int
  iso.int[iso.int>max.inf.day*num.int] = max.inf.day*num.int
  
  # amend contact sequences
  iso.seq = matrix(c(iso.int, rep(max.inf.day*num.int, length(iso.int))), byrow = T, nrow=2)
  iso.seq[2,] = iso.seq[2,]-iso.seq[1,]
  iso.seq = as.vector(iso.seq)
  iso.seq = rep(rep(c(1,0),length(iso.seq)/2), times = iso.seq)
  iso.matrix = matrix(iso.seq, ncol=nsim)
  
  ct.matrix = ct.matrix*iso.matrix
  
  # multiply inf and ct matrix by their respective elements
  inf.matrix = inf.matrix*ct.matrix*beta 
  # beta scales the probability of infection between a susceptible and infected pair
  # beta equals 1 when simulating Ferretti et al skew logistic model as 
  # pdf is the distribution on times of infection not the total force of infection 
 
  if(inf.func.name != 'spline.covid.wild.ferretti'){
    # P(infection | contact)
    p.inf = sapply(1:nsim, function(x){
    # P(not infected)
    p.non.inf = exp(-sum(inf.matrix[,x]))
    # P(infected at interval | not infected previously)
    p.inf.interval = (1-exp(-inf.matrix[,x]))*c(1,exp(-cumsum(inf.matrix[,x])[1:(length(inf.matrix[,x])-1)]))
    
    return(c(p.inf.interval,p.non.inf))  
    })
    # 8641 (intervals + NA) x 10 (nsim) matrix 
  }
  
  if(inf.func.name == 'spline.covid.wild.ferretti'){
    p.inf = inf.matrix
    p.inf = rbind(p.inf, rep(0,ncol(p.inf)))
  }
  
  
  
  # infection outcomes
  gen = sapply(1:nsim, function(x){

    # sample(c(1:(max.inf.day/dt)*dt,NA), 1, prob = p.inf[,x], replace = T)
    if(p.inf[nrow(p.inf),x]==1){
      NA
    }else{
      sample(c(1:(max.inf.day/dt)*dt), 1, prob = p.inf[-nrow(p.inf),x], replace = T)
    }

  })
  
  
  gen.int.end = gen.int.start + gen/dt
  
  
  
  # incubation period of secondary case
  inc.sec = inc.func[[inc.func.name]]
  inc.sec = inc.sec(nsim, inc.pars)
  
  #while(any(inc.sec>max.inf.day)){ inc.sec[inc.sec>max.inf.day] = inc.f(length(inc.sec[inc.sec>max.inf.day]), inc.pars) }
  
  tost = gen-inc
  si = inc.sec+tost
  si.int.end = ceiling(si/dt) + si.int.start 
  
  
  # prob presymptomatic transmission
  p.presymp.trans=sum(gen.int.end<si.int.start)/nsim
  
  # prob transmission
  # p.trans = sum(!is.na(gen))
  # p.trans = p.trans/nsim
  p.trans.overall = 1-sum(p.inf[nrow(p.inf),])/nsim
  p.trans = 1-p.inf[nrow(p.inf),]
  
  # msi = mean(si, na.rm = T)
  # mean generation time per infector
  # mgen = mean(gen, na.rm = T)
  
  # track difference in the generation distribution by days contact tracing
  # track difference in the serial interval distribution by days contact tracing
  # find the respective day interval and subtract
  # gen.day = ceiling(gen.int.end/(1/dt))-ceiling(gen.int.start/(1/dt))
  # si.day = ceiling(si.int.end/(1/dt))-ceiling(si.int.start/(1/dt))
 
  list(# mgen=mgen,
       # msi=msi,
       gen=gen,
       si=si,
       tost=tost,
       # gen.day=gen.day,
       # si.day=si.day,
       p.trans=p.trans,
       p.trans.overall=p.trans.overall,
       p.presymp.trans=p.presymp.trans,
       inc=inc,
       onset.to.iso=onset.to.iso)
}


#' incubation period function
#' @param x: number of simulations
#' @param pars: parameters
inc.func<-list(
  lnorm=function(x,pars){
    rlnorm(x, meanlog = pars[1], sdlog = pars[2])
  },
  lnorm.mix=function(x,pars){
    (rlnorm(x,meanlog = 1.621, sdlog = 0.418)+
     rlnorm(x,meanlog = 1.425, sdlog = 0.669)+
     rlnorm(x,meanlog = 1.57, sdlog = 0.65)+
     rlnorm(x,meanlog = 1.53, sdlog = 0.464)+
     rlnorm(x,meanlog = 1.611, sdlog = 0.472)+
     rlnorm(x,meanlog = 1.54, sdlog = 0.47)+
     rlnorm(x,meanlog = 1.857, sdlog = 0.547))/7
  },
  weibull=function(x,pars){
    rweibull(x, shape = pars[1], scale = pars[2])
  },
  normal=function(x,pars){
    rnormTrunc(x, mean = pars[1], sd = pars[2], min=1)
  }
)


#' infectiousness period function
#' @param peak: day of peak infectiousness
#' @param max: day of maximum infectiousness period
inf.func<-list(
  spline.covid.wild=function(peak, max, dt){
    spline.density = splinefunH(c(0,peak,peak+10), c(0,1,0), c(0,0,0))
    spline.density(1:(max/dt)*dt)
  },
  spline.covid.delta=function(peak, max, dt){
    spline.density = splinefunH(c(0,peak,peak+18), c(0,1,0), c(0,0,0))
    spline.density(1:(max/dt)*dt)
  },
  spline.measles=function(peak, max, dt){
    spline.density = splinefunH(c(0,peak-4,peak,peak+4), c(0,0,1,0), c(0,0,0,0))
    spline.density(1:(max/dt)*dt)
  },
  spline.smallpox=function(peak, max, dt){
    spline.density = splinefunH(c(0,peak-dt,peak,peak+3,peak+14), c(0,0,0.25,1,0), c(0,0,0,0,0))
    spline.density(1:(max/dt)*dt)
  },
  spline.flu=function(peak, max, dt){
    if(peak>1) spline.density = splinefunH(c(0,peak-1,peak,peak+6), c(0,0,1,0), c(0,0,0,0))
    if(peak==1) spline.density = splinefunH(c(peak-1,peak,peak+6), c(0,1,0), c(0,0,0))
    spline.density(1:(max/dt)*dt)
   
  },
  normal=function(peak, max, dt){
    spline.density = dnorm(1:(max/dt)*dt, mean=peak, sd=1)
    spline.density = spline.density/max(spline.density)
  },
  spline.covid.wild.ferretti=function(peak, max, dt){
   
    pars.skewlog.wt=c(-3.6658073,  0.5441884,  1.7480858) # mu, sigma, alpha, derived by fitting model to transmission pairs 
    
    t=seq(-peak, max-peak, dt) # 5 mins intervals for 30 days before and 30 days after
    spline.density =(t>(-peak) & t<0)*ferretti(t/(peak/5.42),pars.skewlog.wt)+(t>=0)*ferretti(t,pars.skewlog.wt)
    spline.density = spline.density/max(spline.density)
    spline.density = spline.density[-1]
    spline.density
  }
)

#' Feretti's skew logistic model, asymmetrical scaling, with presymptomatic period scaled to allow for 
#' individuals with longer incubation period to have proportionally earlier and longer presymptomatic infectious period. 
#' @param x: number of simulations
#' @param pars: parameters
ferretti=function(x,pars){dglogis(x, location = pars[1], scale = exp(pars[2]), shape = exp(pars[3]))}


#' onset to isolation period function
#' @param x: number of simulations
#' @param pars: parameters
onset.to.iso.func<-list(
  weibull=function(x,pars){
    rweibull(x, shape = pars[1], scale = pars[2])
  },
  none=function(x,pars){
    rep((10^6),x)
  },
  uniform=function(x,pars){
    runif(x, min = pars[1], max = pars[2])
  }
)
