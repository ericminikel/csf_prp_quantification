options(stringsAsFactors=FALSE)
library(reshape2)
library(sqldf)
setwd('~/d/sci/src/csf_prp_quantification/')

stdcol = '#9D1309'
csfcol = '#3333FF'

source('src/shared_functions.R')

convert_to_mgml = function(absorbance, dilution_factor=1) {
  estimated_conc = (slope*absorbance + intercept) * dilution_factor
  return (estimated_conc)
}

plates = read.table('data/dc/meta/plates.tsv',sep='\t',header=TRUE)

for (plate in plates$plateno) {
  
  datafile = plates$datafile[plates$plateno==plate]
  metafile = plates$metafile[plates$plateno==plate]
  
  datafile_fullpath = paste('data/dc/raw/',datafile,sep='')
  
  if(grepl('\\.txt',datafile,ignore.case=TRUE)) {
    rawdata = process_spectramax_data(datafile_fullpath)
  } else if (grepl('\\.csv',datafile,ignore.case=TRUE)) {
    rawdata = process_optima_data(datafile_fullpath)
  } else {
    stop('unknown file extension: ',datafile)
  }
  
  meta = read.table(paste('data/dc/meta/',metafile,sep=''),sep='\t',header=T)
  
  rawdata$plate = plate
  meta$plate = plate
  
  if (castable_to_integer(plate)) {
    plate_prefix = formatC(as.integer(plate),width=2,flag='0')
  } else {
    plate_prefix = plate
  }
  
  data = merge(rawdata,meta,by=c("plate","row","col"),sort=FALSE)
  # merge's sort option sorts as if they were characters, so 1, 10, 11, 12, 2, etc.
  # so I set sort=FALSE in the merge above and then sort separately:
  data = data[order(data$row, data$col),]
  
  # log-log fit is TERRIBLE - use linear for DC assay
  m = lm(as.numeric(detail) ~ a750, data=subset(data, stype=='standard')) 
  summary(m)
  intercept = coefficients(m)[1]
  slope = coefficients(m)[2]
  
  subtitle = paste('plate ',plate,': ',datafile,sep='')

  # standard curve
  png(paste('data/dc/qc/',formatC(plate,width=2,flag='0'),'_standards.png',sep=''),width=600,height=400,res=100)
  par(mar=c(4,4,4,6))
  plot(NA, NA, xlim=c(0,2), ylim=c(0,.6), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
  xats = as.numeric(unique(data$detail[data$stype=='standard']))
  axis(side=1, at=xats, lwd.ticks=1, lwd=0, cex.axis=.8)
  axis(side=2, at=(0:4)/10, lwd.ticks=1, lwd=0, las=2)
  mtext(side=1, line=2.5, text='Standard concentration (mg/mL)')
  mtext(side=2, line=2.5, text='Absorbance (arbitrary units)')
  title(main='Standard curve')
  mtext(side=3,line=0.5,text=subtitle,cex=.7)
  abline(h=0,lwd=2)
  abline(v=0,lwd=2)
  # fit model
  x = (1:600)/1000
  f_of_x = intercept + slope*x
  points(x = f_of_x, y = x, type='l', lwd=.5)
  # actual points
  points(as.numeric(data$detail[data$stype=='standard']),data$a750[data$stype=='standard'], pch=20, col=stdcol)
  dev.off()
  
  data = subset(data, detail != '') # remove blanks
  
  data$mgml = convert_to_mgml(data$a750, dilution_factor=data$dilution)
  
  llq = .1
  ulq = 1.6
  
  # annoyingly this has to be done in a second step because mgml might still be NA above
  data$mgml[data$mgml/data$dilution < llq] = llq*data$dilution[data$mgml/data$dilution < llq]
  data$flag[data$mgml/data$dilution < llq] = 'LLQ'
  data$mgml[data$mgml/data$dilution > ulq] = ulq*data$dilution[data$mgml/data$dilution > ulq]
  data$flag[data$mgml/data$dilution > ulq] = 'ULQ'  
  
  data$flag[is.na(data$flag)] = ''
  
  # well-level detail table
  write.table(data,paste('data/dc/processed/',formatC(plate,width=2,flag='0'),'.tsv',sep=''),row.names=F,col.names=T,sep='\t',quote=F)
  
  # summary-level info
  standard_curve = suppressWarnings(as.numeric(unique(data$detail[data$stype=='standard'])))
  llq = min(standard_curve[standard_curve > 0])
  ulq = max(standard_curve[standard_curve > 0])
  llq_fluor = suppressWarnings(mean(data$a750[data$stype=='standard' & as.numeric(data$detail)==llq],na.rm=TRUE)) ##
  ulq_fluor = suppressWarnings(mean(data$a750[data$stype=='standard' & as.numeric(data$detail)==ulq],na.rm=TRUE)) ##
  smry = data.frame(sample = unique(data$detail[data$stype=='sample']))
  smry$mgml_av = as.numeric(NA)
  smry$flag = ''
  for (i in 1:nrow(smry)) {
    relevant_rows = data$stype == 'sample' & data$detail == smry$sample[i]
    in_range = relevant_rows & data$a750 > llq_fluor & data$a750 < ulq_fluor
    if (sum(in_range) > 0) {
      smry$mgml_av[i] = mean(data$mgml[in_range])
    } else if (all(data$a750[relevant_rows] > ulq_fluor)) {
      smry$mgml_av[i] = ulq 
      smry$flag[i] = 'ULQ'
    } else if (all(data$a750[relevant_rows] < llq_fluor)) {
      smry$mgml_av[i] = llq 
      smry$flag[i] = 'LLQ'
    }
  }
  write.table(smry,paste('data/dc/processed/',formatC(plate,width=2,flag='0'),'_summary.tsv',sep=''),row.names=F,col.names=T,sep='\t',quote=F)
  
}