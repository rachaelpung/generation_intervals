load('data/el.10m.RData')

# generate contact sequences
ct.pair=unique(el[ct.grp==1 & ct==1, .(id.x, id.y)])
max.ct.day=50
set.seed(123)
ct.list=sapply(1:100, function(x){
  
  # simulate full contact
  # ct=rep(1, (max.ct.day+1)*288)
  
  ct=0
  while(all(ct==0)){

    row = sample(1:ct.pair[,.N], 1, replace=T)

    day = sample(seq(Sys.Date(),(Sys.Date()+6), by = "day"),1, replace=T)
    day = weekdays(seq(day, day+max.ct.day, 1), abbreviate=T)

    el.pair = data.table(ct.pair[row], day)
    el.pair[day %in% c('Sat', 'Sun'), day:='Sat']
    el.pair[day != 'Sat', day:=sample(c('Thu', 'Fri'),1,replace=T),
            by=seq_len(nrow(el.pair[day != 'Sat']))]

    int.list = list(Thu = 1:288, Fri = 289:576, Sat = 577:864)
    int.list = sapply(1:nrow(el.pair), function(x){ int.list[[ el.pair[x,day] ]]   })
    int.list = as.vector(int.list)

    el.pair = el.pair[rep(seq_len(nrow(el.pair)), each=288)]
    el.pair[, interval:=int.list]
    el.pair[el, ct:=i.ct, on=c(id.x='id.x', id.y='id.y', day='day', interval='interval')]


    # for contacts that meet on one or more days a week at fixed intervals
    # ct.days = rep(0,14)
    # ct.days[sample(1:14,3,replace=T)] = 1
    # el.pair$ct = el.pair$ct*rep(rep_len(ct.days, 51), each=288)
    
    # ct.days = rep(0,7)
    # ct.days[sample(1:7,6,replace=F)] = 1
    # ct.days=c(ct.days,ct.days)
    # ct.days[sample(which(ct.days==0),1,replace = T)] = 1
    
    # for contacts that meet in the mid week
    ct.days = rep(0,7)
    ct.days[4] = 1
    
    el.pair$ct = el.pair$ct*rep(rep_len(ct.days, 51), each=288)
    

    ct=el.pair$ct
  }
  ct
})

dim(ct.list)
# 14688 (intervals) x 1000 (nsim) matrix
# 51*288

ct.list.hh.wk.mid = copy(ct.list)
save(ct.list.hh.wk.mid, file='data/ct.list.hh.wk.mid.10m.RData')

rm(ct.list, ct.list.hh.wk.mid)
dim(ct.list)
