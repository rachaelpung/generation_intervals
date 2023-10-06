source('code/step0_load_library.R')

# load edgelist and timelist
tl <- fread('data/haslemere_time.csv')
el <- fread('data/haslemere_edge.csv')

setnames(tl, c('interval', 'date.time'))
setnames(el, c('interval', 'id.x', 'id.y', 'dist'))
setorder(el, id.x, id.y, interval)

# distance threshold for a contact
dist.thres=10

# identify potential household contacts
tl[, hour.min:= as.numeric(paste(substr(date.time, 17, 18), substr(date.time, 20, 21), sep=''))]
tl[, day:=substr(date.time, 1, 3)]

tl[hour.min<800, loc:='hh.day']
tl[hour.min>=2000, loc:='hh.night']
tl[!(hour.min<800 | hour.min>=2000), loc:='nhh.day']

el[tl, loc:=i.loc, on=c(interval='interval')]
el[tl, date.time:=i.date.time, on=c(interval='interval')]
el[, day:=substr(date.time, 1, 3)]

el.hh=el[loc %in% c('hh.day', 'hh.night') & dist<=dist.thres, uniqueN(loc), by=.(id.x,id.y,day)]
el.hh=el.hh[,sum(V1), by=.(id.x,id.y)]
table(el.hh$V1)

# contact group, 1: household, 2: non-household
el.hh[V1>=5, ct.grp:=1]
el.hh[V1<5, ct.grp:=2]

# undirected edgelist
el.hh=rbind(el.hh,el.hh)
el.hh[(.N/2+1):.N, `:=` (id.x=el.hh[1:(.N/2),id.y],
                         id.y=el.hh[1:(.N/2),id.x])]

el[el.hh, ct.grp:=i.ct.grp, on=c(id.x='id.x', id.y='id.y')]
el[is.na(ct.grp), ct.grp:=2]
el[dist<=dist.thres, ct:=1]

# amend tl to add in intervals from 11pm to 7am 
tmp=data.table(hour.min=c(seq(2300, 2355, 5),
                          rep(seq(0,600,100), each=12) + rep(seq(0,55,5))),
               loc=c(rep('hh.night', 12), rep('hh.day',84)))
tmp=tmp[rep(seq_len(.N), times=3)]
tmp[,day:=rep(c('Thu','Fri','Sat'), each=.N/3)]


tl=rbindlist(list(tl,tmp), use.names = T, fill = T)
tl[,day:=factor(day, levels=c('Thu', 'Fri', 'Sat'))]
setorder(tl, day, hour.min)
tl[,interval:=seq_len(.N)]
rm(tmp)

# amend el to add in intervals from 11pm to 7am
el[tl, interval:=i.interval, on=c(date.time='date.time')]

tmp=unique(el[,.(id.x, id.y,ct.grp)])
tmp=tmp[rep(seq_len(.N), each=nrow(tl))]
tmp[,interval:=rep(tl$interval, times=.N/nrow(tl))]

el = merge(el, tmp, all = T, by=c('id.x'='id.x', 'id.y'='id.y', 'interval'='interval', 'ct.grp'='ct.grp'))
el[is.na(ct), ct:=0] # no hh and non hh contact from 11pm-7am, assume to be asleep
el=el[,.(id.x, id.y, interval, dist, ct, ct.grp, day)]
el[tl, day:=i.day, on=c(interval='interval')]
el[is.na(dist), dist:=0]
rm(tmp, el.hh)

# find number of graph components and size of each component
nl.hh = data.table(id=sort(unique(c(el[ct.grp==1, id.x],el[ct.grp==1, id.y]))))

graph = graph_from_data_frame(unique(el[ct.grp==1, .(id.x, id.y)]), directed = TRUE, vertices = nl.hh)
graph = as.undirected(graph, mode ='collapse', edge.attr.comb="first")

table(components(graph)$csize)
# 2  3  4  5  6  7  8 12 
# 39 17  8  2  4  1  1  2 

components(graph)$no
# 54

save(el, file='data/el.10m.RData')



