options(stringsAsFactors=FALSE)
library(reshape2)
library(sqldf)
library(plyr)
setwd('~/d/sci/src/csf_prp_quantification/')

source('src/shared_functions.R')

convert_to_ngml = function(absorbance, dilution_factor=1) {
  estimated_conc = exp(slope*log(absorbance) + intercept) * dilution_factor
  return (estimated_conc)
}

stdcol = '#9D1309'
concol = '#444456'
csfcol = '#3333FF'
spkcol = '#4CBB17'
lincol = '#777777'

hucol = '#239ACD'
ocol = '#28170D'

espresso = '#28170D'

standard_curve = c(1,2,4,5,10,20)

plates = read.table('data/elisa/meta/plates.tsv',sep='\t',header=TRUE)

for (plate in plates$plateno) {

datafile = plates$datafile[plates$plateno==plate]
metafile = plates$metafile[plates$plateno==plate]

datafile_fullpath = paste('data/elisa/raw/',datafile,sep='')

if(grepl('\\.txt',datafile,ignore.case=TRUE)) {
  rawdata = process_spectramax_data(datafile_fullpath)
} else if (grepl('\\.csv',datafile,ignore.case=TRUE)) {
  rawdata = process_optima_data(datafile_fullpath)
} else {
  stop('unknown file extension: ',datafile)
}

meta = read.table(paste('data/elisa/meta/',metafile,sep=''),sep='\t',header=T)

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

data$dilution[is.na(data$dilution)] = 1

#### QC plots

subtitle = paste('plate ',plate,': ',datafile,sep='')

# plate heatmap
mat = daply(rawdata, .(row, col), function(x) x$a450_620)
png(paste('data/elisa/qc/',plate_prefix,'_heatmap.png',sep=''),width=600,height=400,res=100)
image(x=1:12, y=1:8, z=t(apply(mat, 2, rev)), col=colorRampPalette(c("black", "red"))(10),axes=FALSE, ann=FALSE)
axis(side=1, at=(1:12), labels=1:12, lwd=0, lwd.ticks=0)
axis(side=2, at=(1:8), labels=LETTERS[8:1], lwd=0, lwd.ticks=0, las=2)
title(main='Heatmap')
mtext(side=3,line=0.5,text=subtitle,cex=.7)
dev.off()

# controls
png(paste('data/elisa/qc/',plate_prefix,'_controls.png',sep=''),width=300,height=400,res=100)
par(mar=c(4,4,4,1))
plot(NA, NA, xlim=c(0,3), ylim=c(0,max(data$a450_620,na.rm=T)), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=2, at=0:3, lwd.ticks=1, lwd=0, las=2)
axis(side=1, at=(1:3)-.5, labels=c('pos hi', 'pos lo', 'blank'), lwd.ticks=0, lwd=0)
mtext(side=2, line=2.5, text='Absorbance (arbitrary units)')
title(main='Kit controls')
mtext(side=3,line=0.5,text=subtitle,cex=.7)
abline(h=0,lwd=2)
abline(v=0,lwd=2)
points(rep(.5, sum(data$detail=='poscon hi')), data$a450_620[data$detail=='poscon hi'], pch=20, col=concol)
points(rep(1.5, sum(data$detail=='poscon lo')), data$a450_620[data$detail=='poscon lo'], pch=20, col=concol)
points(rep(2.5, sum(data$detail=='negcon')), data$a450_620[data$detail=='negcon'], pch=20, col=concol)
points(rep(2.5, sum(data$detail=='blocking')), data$a450_620[data$detail=='blocking'], pch=20, col=concol)
dev.off()

# throw out any standard curve points that are less than the negative control (presumptive pipetting errors)
data$above_baseline = TRUE
if (any(c('negcon','blocking') %in% data$detail)) {
  data$above_baseline = data$a450_620 > mean(data$a450_620[data$detail %in% c('negcon','blocking')])
}

# fit model to standard curve
m = lm(log(as.numeric(detail)) ~ log(a450_620), data=subset(data, stype=='standard' & above_baseline))
summary(m)
intercept = coefficients(m)[1]
slope = coefficients(m)[2]

# standard curve
png(paste('data/elisa/qc/',plate_prefix,'_standards.png',sep=''),width=600,height=400,res=100)
par(mar=c(4,4,4,6))
plot(NA, NA, xlim=c(0,21), ylim=c(0,3.5), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=standard_curve, lwd.ticks=1, lwd=0, cex.axis=.8)
axis(side=2, at=0:3, lwd.ticks=1, lwd=0, las=2)
mtext(side=1, line=2.5, text='Standard concentration (ng/mL)')
mtext(side=2, line=2.5, text='Absorbance (arbitrary units)')
title(main='Standard curve')
mtext(side=3,line=0.5,text=subtitle,cex=.7)
abline(h=0,lwd=2)
abline(v=0,lwd=2)
# fit model
x = (1:400)/100
f_of_x = exp(intercept + slope*log(x))
points(x = f_of_x, y = x, type='l', lwd=.5)
# actual points
points(as.numeric(data$detail[data$stype=='standard']),data$a450_620[data$stype=='standard'], pch=20, col=stdcol)
if ('poscon hi' %in% data$detail) {
  abline(h=mean(data$a450_620[data$detail=='poscon hi']), lty=3, lwd=2, col=concol)
  axis(side=4, at=mean(data$a450_620[data$detail=='poscon hi']), labels=c('poscon hi'), lwd.ticks=0, lwd=0, las=2, col.axis=concol)
}
if ('poscon lo' %in% data$detail) {
  abline(h=mean(data$a450_620[data$detail=='poscon lo']), lty=3, lwd=2, col=concol)
  axis(side=4, at=mean(data$a450_620[data$detail=='poscon lo']), labels=c('poscon lo'), lwd.ticks=0, lwd=0, las=2, col.axis=concol)
}
if (any(c('negcon','blocking') %in% data$detail)) {
  abline(h=mean(data$a450_620[data$detail %in% c('negcon','blocking')]), lty=3, lwd=2, col=concol)
  axis(side=4, at=mean(data$a450_620[data$detail %in% c('negcon','blocking')]), labels=c('blank'), lwd.ticks=0, lwd=0, las=2, col.axis=concol)
}
dev.off()

# calculate the inferred concentration for each sample and store
data$ngml = convert_to_ngml(data$a450_620, dilution_factor = data$dilution)

# occasionally because of A620 normalization the A450-A620 may be <0, in which case ngml is NaN.
# set these ngml values to 0.
data$ngml[is.na(data$ngml) & data$a450_620 < 0] = 0.0

llq = 1
ulq = 20
llq_fluor = mean(data$a450_620[data$stype=='standard' & data$detail=='1']) ## LLQ = 1
ulq_fluor = mean(data$a450_620[data$stype=='standard' & data$detail=='20']) ## ULQ = 20

data$ngml_trunc = data$ngml # make a copy of the ngml column, but with range only in [LLQ, ULQ]
data$ngml_trunc[data$a450_620 < llq_fluor] = llq*data$dilution[data$a450_620 < llq_fluor]
data$flag[data$a450_620 < llq_fluor] = 'LLQ'
data$ngml_trunc[data$a450_620 > ulq_fluor] = ulq*data$dilution[data$a450_620 > ulq_fluor]
data$flag[data$a450_620 > ulq_fluor] = 'ULQ'

# annoyingly this has to be done in a second step because ngml might still be NA above
data$ngml_trunc[data$ngml/data$dilution < llq] = llq*data$dilution[data$ngml/data$dilution < llq]
data$flag[data$ngml/data$dilution < llq] = 'LLQ'
data$ngml_trunc[data$ngml/data$dilution > ulq] = ulq*data$dilution[data$ngml/data$dilution > ulq]
data$flag[data$ngml/data$dilution > ulq] = 'ULQ'  

data$flag[is.na(data$flag)] = ''

# table with well-level detail
write.table(data,paste('data/elisa/processed/',plate_prefix,'.tsv',sep=''),row.names=F,col.names=T,sep='\t',quote=F)

# summary table with sample-level detail. average for each sample
# - throw out observations outside dynamic range of assay
# - if all observations outside dynamic range, set value to limit but also set flag to indicate LLQ/ULQ violated
smry = data.frame(sample = unique(data$detail[data$stype=='sample']))
smry$ngml_av = as.numeric(NA)
smry$se_mean = as.numeric(NA)
smry$flag = ''
for (i in 1:nrow(smry)) {
  relevant_rows = data$stype == 'sample' & data$detail == smry$sample[i]
  in_range = relevant_rows & data$a450_620 > llq_fluor & data$a450 < ulq_fluor
  if (sum(in_range) > 0) {
    # handle a special case first: if <100% of replicates are in dynamic range even at the one valid dilution,
    # then base the estimate on averaging the measured value and the LLQ or ULQ as appropriate.
    dilutions_in_range = unique(data$dilution[in_range])
    if (length(dilutions_in_range)==1 & sum(in_range) < sum(relevant_rows)) {
      at_relevant_dilution = relevant_rows & data$dilution %in% dilutions_in_range
      smry$ngml_av[i] = mean(data$ngml_trunc[at_relevant_dilution])
      smry$se_mean[i] = sd(data$ngml_trunc[at_relevant_dilution]) / sqrt(sum(at_relevant_dilution))
    } else { # most of the time, just average all the observations that do fall within dynamic range
      smry$ngml_av[i] = mean(data$ngml[in_range])
      smry$se_mean[i] = sd(data$ngml[in_range]) / sqrt(sum(in_range))
    }
  } else if (all(data$a450_620[relevant_rows] > ulq_fluor)) {
    smry$ngml_av[i] = ulq * max(data$dilution[relevant_rows])
    smry$se_mean[i] = 0
    smry$flag[i] = 'ULQ'
  } else if (all(data$a450_620[relevant_rows] < llq_fluor)) {
    smry$ngml_av[i] = llq * min(data$dilution[relevant_rows])
    smry$se_mean[i] = 0
    smry$flag[i] = 'LLQ'
  }
}
write.table(smry,paste('data/elisa/processed/',plate_prefix,'_summary.tsv',sep=''),row.names=F,col.names=T,sep='\t',quote=F)

}