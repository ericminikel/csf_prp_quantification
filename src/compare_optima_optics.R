options(stringsAsFactors=FALSE)
library(reshape2)
library(sqldf)
library(plyr)
setwd('~/d/sci/src/csf_prp_quantification/')

folder = 'data/hb/raw/'

# same but for the Fluostar Optima platereader
process_optima_data = function(path) {
  lines_to_skip = grep("Well Row,Well Col",readLines(path))-1 # figure out where the data actually starts
  rawdata = read.table(path,skip=lines_to_skip,header=TRUE,sep=',') # read it in
  rawdata = rawdata[,c(-3,-(ncol(rawdata)))] # remove "Content" column and empty column at end
  colnames(rawdata)[c(1,2)] = c('row','col')
  for (i in 3:ncol(rawdata)) {
    colnames(rawdata)[i] = paste('a',gsub('[A-Za-z\\.]+','',gsub('BP.*','',colnames(rawdata)[i])),sep='')
  }
  if (all(c('a450','a620') %in% colnames(rawdata))) {
    rawdata$a450_620 = rawdata$a450 - rawdata$a620
    rawdata = rawdata[,-which(colnames(rawdata) %in% c('a450','a620'))]
  }
  return (rawdata)
}

a620_new = process_optima_data(paste(folder,'hbELISA00002-a620-newoptic.CSV',sep=''))
a620_old = process_optima_data(paste(folder,'hbELISA00002-a620-oldoptic.CSV',sep=''))
a450_old1 = process_optima_data(paste(folder,'hbELISA00002-oldoptic-scan1.CSV',sep=''))
a450_old2 = process_optima_data(paste(folder,'hbELISA00002-oldoptic-scan2.CSV',sep=''))
a450_new1 = process_optima_data(paste(folder,'hbELISA00002-scan1.CSV',sep=''))
a450_new2 = process_optima_data(paste(folder,'hbELISA00002-scan2.CSV',sep=''))

# very tight same-well agreement between two scans on new optic
plot(a450_new1$a450, a450_new2$a450, pch=20, xlim=c(0,3), ylim=c(0,3))

# and old optic
plot(a450_old1$a450, a450_old2$a450, pch=20, xlim=c(0,3), ylim=c(0,3))

# slightly less good between-optic agreement
plot(a450_new1$a450, a450_old1$a450, pch=20, xlim=c(0,3), ylim=c(0,3))

mean_cv = function(x,y) {
  stopifnot(length(x)==length(y))
  means = numeric(length(x))
  sds = numeric(length(x))
  for (i in 1:length(x)) {
    means[i] = mean(c(x[i],y[i]))
    sds[i] = sd(c(x[i],y[i]))
  }
  return (mean(sds/means))
}

# within-optic CV:
mean_cv(a450_new1$a450, a450_new2$a450) # 1.9%
mean_cv(a450_old1$a450, a450_old2$a450) # 1.8%
mean_cv(a450_new1$a450, a450_old1$a450) # 13.4% - ouch
mean_cv(a450_new1$a450, a450_old2$a450) # 14.4%
mean_cv(a450_new2$a450, a450_old1$a450) # 11.8%
mean_cv(a450_new2$a450, a450_old2$a450) # 12.8%

# and what about A620?
plot(a620_new$a620, a620_old$a620, pch=20, xlim=c(0,.5), ylim=c(0,.5))
# hardly any agreement at all!

# maybe controlling for A620 fixes problems between optics?
plot(a450_new1$a450, a450_old1$a450, pch=20, xlim=c(0,3), ylim=c(0,3))
plot(a450_new1$a450 - a620_new$a620, a450_old1$a450 - a620_old$a620, pch=20, xlim=c(0,3), ylim=c(0,3))
# looks slightly better.
# the mean_cv turns out not to be meaningful here though because after subtracting A620 some values are negative
# or very close to zero, which causes the CV to blow up

# in any case, the new optic is at least as _internally_ consistent as the old one was, so no reason not to just use the 
# new one from now on.

# also try the A750 readings from total protein:
folder = 'data/dc/raw/'
a750_new1 = process_optima_data(paste(folder,'totprot4-2016-11-17.CSV',sep=''))
a750_new2 = process_optima_data(paste(folder,'totprot4-read2-2016-11-17.CSV',sep=''))
plot(a750_new1$a750, a750_new2$a750, pch=20, xlim=c(0,1), ylim=c(0,1))
mean_cv(a750_new1$a750, a750_new2$a750)
# mean CV only 1.5% - looks good



