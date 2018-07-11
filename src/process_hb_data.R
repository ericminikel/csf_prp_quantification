options(stringsAsFactors=FALSE)
library(reshape2)
library(sqldf)
library(drc)
setwd('~/d/sci/src/csf_prp_quantification/')

stdcol = '#9D1309'
csfcol = '#3333FF'


source('src/shared_functions.R')

convert_to_ngml = function(absorbance, dilution_factor=1) {
  estimated_conc = exp(slope*log(absorbance) + intercept) * dilution_factor
  return (estimated_conc)
}

plates = read.table('data/hb/meta/plates.tsv',sep='\t',header=TRUE)


for (plate in plates$plateno) {
  
  datafile = plates$datafile[plates$plateno==plate]
  metafile = plates$metafile[plates$plateno==plate]
  
  datafile_fullpath = paste('data/hb/raw/',datafile,sep='')
  
  if(grepl('\\.txt',datafile,ignore.case=TRUE)) {
    rawdata = process_spectramax_data(datafile_fullpath)
  } else if (grepl('\\.csv',datafile,ignore.case=TRUE)) {
    rawdata = process_optima_data(datafile_fullpath)
  } else {
    stop('unknown file extension: ',datafile)
  }
  
  meta = read.table(paste('data/hb/meta/',metafile,sep=''),sep='\t',header=T)
  
  if (nrow(meta) == 0) {
    next # failed plates get an empty meta file
  }
  
  rawdata$plate = plate
  meta$plate = plate
  
  if (castable_to_integer(plate)) {
    plate_prefix = formatC(as.integer(plate),width=2,flag='0')
  } else {
    plate_prefix = plate
  }
  
  data = merge(rawdata,meta,by=c("plate","row","col"),sort=FALSE)
  data = subset(data, detail != '') # remove unused wells
  # merge's sort option sorts as if they were characters, so 1, 10, 11, 12, 2, etc.
  # so I set sort=FALSE in the merge above and then sort separately:
  data = data[order(data$row, data$col),]
  
  # the log-log model fit is pretty bad for hemoglobin; the kit instructions recommend
  # a 4-point curve fit, but I tried that and it is even worse(!) so I am just doing 
  # linear interpolation with approx(). 
  # I also lop off the .27 ng/mL standard because it is very close
  # to the A450 for the blank, whereas the rest of the standard curve (at least on most plates)
  # has some kind of dynamic range.
  fitdata = data[data$stype=='standard', c('detail', 'a450')]
  fitdata$ngml = as.numeric(fitdata$detail)
  fitsum = sqldf("select ngml, avg(a450) a450 from fitdata where ngml > .27 group by 1 order by 1;")
  
  subtitle = paste('plate ',plate,': ',datafile,sep='')

  # standard curve
  png(paste('data/hb/qc/',formatC(plate,width=2,flag='0'),'_standards.png',sep=''),width=600,height=400,res=100)
  par(mar=c(4,4,4,6))
  plot(NA, NA, xlim=c(.1,250), ylim=c(0,3), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i', log='x')
  xats = as.numeric(unique(data$detail[data$stype=='standard']))
  axis(side=1, at=xats, lwd.ticks=1, lwd=0, cex.axis=.8)
  axis(side=2, at=0:4, lwd.ticks=1, lwd=0, las=2)
  mtext(side=1, line=2.5, text='Standard concentration (ng/mL)')
  mtext(side=2, line=2.5, text='Absorbance (arbitrary units)')
  title(main='Standard curve')
  mtext(side=3,line=0.5,text=subtitle,cex=.7)
  abline(h=0,lwd=2)
  abline(v=0,lwd=2)
  # fit model
  x = (1:4000)/1000
  f_of_x = approx(x=fitsum$a450, y=fitsum$ngml, xout=x)$y
  points(x = f_of_x, y = x, type='l', lwd=.5)
  # actual points
  points(as.numeric(data$detail[data$stype=='standard']),data$a450[data$stype=='standard'], pch=20, col=stdcol)
  dev.off()

  llq = 0.82
  ulq = 200
  llq_fluor = fitsum$a450[fitsum$ngml==llq]
  ulq_fluor = fitsum$a450[fitsum$ngml==ulq]
  
  data$ngml = approx(x=fitsum$a450, y=fitsum$ngml, xout=data$a450)$y * data$dilution
  data$ngml[data$a450 < llq_fluor] = llq*data$dilution[data$a450 < llq_fluor]
  data$flag[data$a450 < llq_fluor] = 'LLQ'
  data$ngml[data$a450 > ulq_fluor] = ulq*data$dilution[data$a450 > ulq_fluor]
  data$flag[data$a450 > ulq_fluor] = 'ULQ'
  
  # annoyingly this has to be done in a second step because ngml might still be NA above
  data$ngml[data$ngml/data$dilution < llq] = llq*data$dilution[data$ngml/data$dilution < llq]
  data$flag[data$ngml/data$dilution < llq] = 'LLQ'
  data$ngml[data$ngml/data$dilution > ulq] = ulq*data$dilution[data$ngml/data$dilution > ulq]
  data$flag[data$ngml/data$dilution > ulq] = 'ULQ'  

  data$flag[is.na(data$flag)] = ''
  
  # table with well-level detail
  write.table(data,paste('data/hb/processed/',formatC(plate,width=2,flag='0'),'.tsv',sep=''),row.names=F,col.names=T,sep='\t',quote=F)

  smry = data.frame(sample = unique(data$detail[data$stype=='sample']))
  smry$ngml_av = as.numeric(NA)
  smry$flag = ''
  for (i in 1:nrow(smry)) {
    relevant_rows = data$stype == 'sample' & data$detail == smry$sample[i]
    in_range = relevant_rows & data$flag == ''
    if (sum(in_range) > 0) {
      smry$ngml_av[i] = mean(data$ngml[in_range])
    } else if (all(data$a450[relevant_rows] > ulq_fluor)) {
      smry$ngml_av[i] = ulq * max(data$dilution[relevant_rows])
      smry$flag[i] = 'ULQ'
    } else if (all(data$a450[relevant_rows] < llq_fluor)) {
      smry$ngml_av[i] = llq * min(data$dilution[relevant_rows])
      smry$flag[i] = 'LLQ'
    }
  }
  write.table(smry,paste('data/hb/processed/',formatC(plate,width=2,flag='0'),'_summary.tsv',sep=''),row.names=F,col.names=T,sep='\t',quote=F)
  
}