#' Run a specified number of simulations with identical parameters
#' @author Rachael Pung
#' @param set.no: parameter set number
#' @param nsim: number of simulations
#' @param net: network
#' @param inc.func.name: name of incubation period function  
#' @param max.inf.day: maximum duration of infectiousness
#' @param beta: scale factor for peak infectiousness
#' @param inf.func.name: name of infectiousness profile function 
#' @param mean.iso: mean delay from symptoms onset to isolation 
#' @param var.iso: variance in delay from symptoms onset to isolation 
#' @param iso.func.name: name of onset to isolation delay  function 
#' @param hh.size: household size
#'
#' @importFrom purrr safely map
#' @importFrom dplyr bind_rows mutate
#' @return data.frame of cases, isolations, quarantines and tests for each simulation 
gi.cluster <- function(set.no, nsim, net, 
                       mean.incub, var.incub, inc.func.name,
                       max.inf.day, beta, inf.func.name,
                       mean.iso, var.iso, iso.func.name, 
                       hh.size){
  
  # run n.sim number of model runs and put them all together in a big data.frame
  res = 1:nsim %>% 
    purrr::map(function(i) gi.cluster.run(net = net, dt=5/(24*60),
                                          mean.incub=mean.incub, var.incub=var.incub,
                                          inc.func.name=inc.func.name,
                                          max.inf.day=max.inf.day, beta=beta, 
                                          inf.func.name=inf.func.name,
                                          mean.iso=mean.iso, var.iso=var.iso,
                                          iso.func.name=iso.func.name,
                                          hh.size=hh.size))
  
  # bind output together and add simulation index
  res = rbindlist(res) 
  res[,`:=` (sim = rep(1:nsim, each=hh.size),
             set.no = set.no)]
  
  res = res[which(!is.na(res$generation)),]
  
  return(res)
}



#' Run a single instance of the branching process model
#' @param dt: size of time step in days
#' @return data.frame of weekly cases or daily cases or the raw case data
#' @export
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr filter mutate group_by summarise arrange
#' @importFrom tibble tibble
#' @importFrom tidyr replace_na
gi.cluster.run <- function(net, dt,
                           mean.incub, var.incub, inc.func.name,
                           max.inf.day, beta, inf.func.name,
                           mean.iso, var.iso, iso.func.name, 
                           hh.size){
  
  # number of intervals in a day
  num.int = 1/dt                    
  
  # incubation period distribution
  inc.pars = unlist(param.convert(inc.func.name, mean.incub, var.incub))
  inc.f = inc.func[[inc.func.name]]
  
  # infectiousness period distribution
  inf.f = inf.func[[inf.func.name]]
  
  # onset to isolation distribution 
  if(!is.na(mean.iso)) iso.pars=unlist(param.convert(iso.func.name, mean.iso, var.iso))
  onset.to.iso.f = onset.to.iso.func[[iso.func.name]]
  
  
  # initial setup
  popsize = hh.size    
  case.data = data.table(id=1:hh.size,
                         exposure = c(0,rep(NA,hh.size-1)),
                         generation = c(1,rep(NA,hh.size-1)),
                         infector = 0,
                         onset = c(inc.f(1,inc.pars),rep(NA,hh.size-1)),
                         isolated.time = c(onset.to.iso.f(1,iso.pars),rep(NA,hh.size-1)),
                         status = c('I',rep('S',hh.size-1)),
                         inc=0)
  case.data[isolated.time!=10^6,isolated.time:=onset+isolated.time]
  
  setkey(case.data, id)
  
  # set initial values for loop indices
  total.cases = 1
  extinct = FALSE
  new.inf.time=0
  
  
  # model loop, ensure 100% attack rate
  while (total.cases < popsize & !extinct) {
    
    i.id = case.data[status=='I' & inc==0, id] # newly infected case of the next generation
    s.id = case.data[status=='S', id]
    
    case.data[id %in% i.id, inc:=onset-exposure]
    inc = case.data[id %in% i.id,inc]
    
    # duration at large since start of infectiousness
    iso = case.data[id %in% i.id,isolated.time]-case.data[id %in% i.id,exposure] 
    iso[iso>max.inf.day] = max.inf.day
    
    
    if(!exists('inf.pairs')){
      inf.pairs = data.table(expand.grid(i.id, s.id))
      colnames(inf.pairs) <- c('inf.id','sec.id')
    }else{
      inf.pairs = inf.pairs[sec.id %in% s.id]
      inf.pairs = rbind(inf.pairs, setnames(data.table(expand.grid(i.id, s.id)), c('inf.id','sec.id')), fill=T)
    }
   
    
    # infectiousness profile for new infector(s)
    inf.matrix = sapply(1:inf.pairs[inf.id %in% i.id,.N], function(x){ 
      inc.x=inc[match(inf.pairs[inf.id %in% i.id]$inf[x], i.id)]
      inf.f(inc.x, max.inf.day, dt) 
    })
    
    # contact pattern for new potential transmission pair(s)
    ct.matrix = sapply(1:inf.pairs[inf.id %in% i.id,.N], function(x){
      
      ct.int.start = sample(85:276, 1, replace=T)
      ct.int.end = (ct.int.start+(24*60/5*max.inf.day-1))
      ct.col = sample(1:ncol(net), 1, replace=T)
      net[ct.int.start:ct.int.end,ct.col]
      
    })
    
    # early isolation
    iso.matrix = sapply(1:inf.pairs[inf.id %in% i.id,.N], function(x){
      
      iso.seq = c(iso[match(inf.pairs[inf.id %in% i.id]$inf[x], i.id)], max.inf.day - iso[match(inf.pairs[inf.id %in% i.id]$inf[x], i.id)])
      iso.seq = round(iso.seq/max.inf.day*8640)
      iso.seq = rep(c(1,0), times = iso.seq)
      
    })
    
    # force of infection 
    inf.matrix = inf.matrix*ct.matrix*iso.matrix*beta
    
    # tabulate force of infection at each time step
    if(!exists('inf.pairs.tab')){
      inf.pairs.tab = data.table(inf.pairs[rep(which(inf.id %in% inf.id), each=max.inf.day/dt)])
      inf.pairs.tab[case.data, inf.exposure:=i.exposure, on=c(inf.id='id')]
      inf.pairs.tab[, time:=inf.exposure+rep(1:(max.inf.day/dt)*dt, times=inf.pairs[inf.id %in% inf.id,.N])]
      inf.pairs.tab[, inf:=c(inf.matrix)]
     
    }else{
      inf.pairs.tab = inf.pairs.tab[sec.id %in% s.id]
      inf.pairs.tab = rbind(inf.pairs.tab, data.table(inf.pairs[rep(which(inf.id %in% i.id), each=max.inf.day/dt)]), fill=T)
      inf.pairs.tab[case.data, inf.exposure:=i.exposure, on=c(inf.id='id')]
      inf.pairs.tab[inf.id %in% i.id, time:=inf.exposure+rep(1:(max.inf.day/dt)*dt, times=inf.pairs[inf.id %in% i.id,.N])]
      inf.pairs.tab[inf.id %in% i.id, inf:=c(inf.matrix)]
      
    }
    
    sec.data = lapply(1:length(s.id), function(x){
      
      event = inf.pairs.tab[sec.id==s.id[x]]
      event = event[sample(1:.N, .N, replace=F)] # jumble event sequence as susceptible may encounter two infectors at the same time
      setorder(event, time) 
      
      # probability of infection for each susceptible at each time step 
      # P(infected at timestep|not infected previously)
      event[,p.inf := (1-exp(-event$inf))*c(1,exp(-cumsum(event$inf)[1:(event[,.N]-1)]))]
      event=event[p.inf!=0 & time>new.inf.time]
      
      
      if(event[,.N]!=0) event = event[sample(1:.N, 1, event$p.inf, replace=F)]

    })
    
  
    sec.data = rbindlist(sec.data)
    sec.data = sec.data[order(time)]
    sec.data=sec.data[1,]
    
    new.inf.time = min(sec.data$time)
    sec.data[case.data, inf.gen:=i.generation, on=c(inf.id='id')]
    
    case.data[sec.data, infector:=i.inf.id, on=c(id='sec.id')]
    case.data[sec.data, exposure:=i.time, on=c(id='sec.id')]
    case.data[sec.data, generation:=i.inf.gen+1, on=c(id='sec.id')]
    case.data[id %in% sec.data$sec, status:='I']
    case.data[id %in% sec.data$sec, onset:=inc.f(.N,inc.pars)]
    case.data[id %in% sec.data$sec, onset:=onset+exposure]
    
    case.data[id %in% sec.data$sec, isolated.time:=onset.to.iso.f(.N,iso.pars)] 
    case.data[id %in% sec.data$sec & isolated.time!=10^6, isolated.time:=isolated.time+onset] 
    
    total.cases = sum(!is.na(case.data$generation))
    extinct = all(case.data$status!='S',na.rm = TRUE)
    
  }
  
  case.data[status=='I' & inc==0, inc:=onset-exposure]
  rm(inf.pairs, inf.pairs.tab)
  case.data
  
  # return
  return(case.data)
  
}

