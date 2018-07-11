options(stringsAsFactors=FALSE)
library(reshape2)
library(sqldf)
setwd('~/d/sci/src/csf_prp_quantification/')

# see http://seananderson.ca/courses/11-multipanel/multipanel.pdf
# use layout() to create unequal size/shape multipanels

expand_range = function(x, by=.5) {
  return ( c(min(x)-by,max(x)+by) )
}

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}

cv = function(x) sd(x)/mean(x)

col25 = '#FF2016'
col50 = '#0001CD'
col100 = '#777777'
csfcol = '#0001CD' # default

samples = read.table('data/samples/samples.tsv',sep='\t',header=TRUE,quote='',comment.char='')

plates = read.table('data/elisa/meta/plates.tsv',sep='\t',header=TRUE)
plates = plates[-which(plates$plate=='6_alt'),] # remove the re-scan of plate 6 done to compare platereaders
plates$plateno = as.integer(plates$plate)
if (exists('elisa')) {
  rm(elisa)
}
for (plate in plates$plateno) {
  filename = paste('data/elisa/processed/',formatC(plate,width=2,flag='0'),'_summary.tsv',sep='')
  if (file.exists(filename)) {
    temp = read.table(filename,sep='\t',header=TRUE,quote='',comment.char='')
    temp$plate = plate
    if (all(is.na(temp$flag))) {
      temp$flag = '' # when all are blank, it reads them as NA which is annoying, so fix it here.
    }
    if (exists('elisa')) {
      elisa = rbind(elisa,temp)
    } else {
      elisa = temp
    }
  }
}
prp_plates = plates
prp_plates$reader = ''
prp_plates$reader[grepl('\\.txt',prp_plates$datafile,ignore.case=TRUE)] = 'SpectraMax'
prp_plates$reader[grepl('\\.csv',prp_plates$datafile,ignore.case=TRUE)] = 'Fluostar Optima'

colnames(elisa)[2:4] = c('prp_ngml', 'prp_se', 'prp_flag')
elisa$prp_flag[is.na(elisa$prp_flag)] = ''

# delete data from failed experiments that were later re-run
elisa = elisa[-which(elisa$plate==3 & elisa$sample %in% samples$id),]



# raw well-by-well values in elisa
if (exists('elisa_raw')) {
  rm(elisa_raw)
}
for (plate in plates$plateno) {
  filename = paste('data/elisa/processed/',formatC(plate,width=2,flag='0'),'.tsv',sep='')
  if (file.exists(filename)) {
    temp = read.table(filename,sep='\t',header=TRUE,quote='',comment.char='')
    temp$plate = plate
    if (exists('elisa_raw')) {
      elisa_raw = rbind(elisa_raw,temp)
    } else {
      elisa_raw = temp
    }
  }
}




plates = read.table('data/dc/meta/plates.tsv',sep='\t',header=TRUE)
rm(dc)
for (plate in plates$plateno) {
  filename = paste('data/dc/processed/',formatC(plate,width=2,flag='0'),'_summary.tsv',sep='')
  if (file.exists(filename)) {
    temp = read.table(filename,sep='\t',header=TRUE,quote='',comment.char='')
    temp$plate = plate
    if (exists('dc')) {
      dc = rbind(dc,temp)
    } else {
      dc = temp
    }
  }
}
dc$flag[is.na(dc$flag)] = ''


plates = read.table('data/hb/meta/plates.tsv',sep='\t',header=TRUE)
rm(hb)
for (plate in plates$plateno) {
  filename = paste('data/hb/processed/',formatC(plate,width=2,flag='0'),'_summary.tsv',sep='')
  if (file.exists(filename)) {
    temp = read.table(filename,sep='\t',header=TRUE,quote='',comment.char='')
    temp$plate = plate
    if (exists('hb')) {
      hb = rbind(hb,temp)
    } else {
      hb = temp
    }
  }
}
hb$flag[is.na(hb$flag)] = ''

# delete data from dilutions we did outside the normal [LLQ, ULQ] range to see how high v1187.3 was, since
# these can't be plotted with the other data
hb = hb[-which(hb$sample=='v1187.3' & hb$plate==6),]


# in case you need to inspect number of distinct plates on which a sample appears
# sqldf("
# select   sample, count(*) n
# from     elisa
# where    sample in (select id from samples)
# group by sample
# order by 2 desc
# ;")


# table with best estimates of everything for each sample
data = data.frame(sample=unique(elisa$sample))

data$prp_ngml = as.numeric(NA)
data$prp_flag = ''

for (i in 1:nrow(data)) {
  this_sample = data$sample[i]
  
  # Handle PrP ELISA
  if (all(elisa$prp_flag[elisa$sample==this_sample] == 'LLQ')) {
    data$prp_flag[i] = 'LLQ'
    data$prp_ngml[i] = min(elisa$prp_ngml[elisa$sample==this_sample])
    data$prp_se[i] = 0
  } else if (all(elisa$prp_flag[elisa$sample==this_sample] == 'ULQ')) {
    data$prp_flag[i] = 'ULQ'
    data$prp_ngml[i] = max(elisa$prp_ngml[elisa$sample==this_sample])
    data$prp_se[i] = 0
  } else {
    data$prp_flag[i] = ''
    data$prp_ngml[i] = mean(elisa$prp_ngml[elisa$sample==this_sample & elisa$prp_flag==''])
    # this next line - taking mean of the SEs - probably isn't the best way to do this, but this affects
    # very few samples
    data$prp_se[i] = mean(elisa$prp_se[elisa$sample==this_sample & elisa$prp_flag==''])
  }
  
  # Handle hb
  if (this_sample %in% hb$sample) {
    if (all(hb$flag[hb$sample==this_sample] == 'LLQ')) {
      data$hb_flag[i] = 'LLQ'
      data$hb_ngml[i] = min(hb$ngml_av[hb$sample==this_sample])
    } else if (all(hb$flag[hb$sample==this_sample] == 'ULQ')) {
      data$hb_flag[i] = 'ULQ'
      data$hb_ngml[i] = max(hb$ngml_av[hb$sample==this_sample])
    } else {
      data$hb_flag[i] = ''
      data$hb_ngml[i] = mean(hb$ngml_av[hb$sample==this_sample & hb$flag==''])
    } 
  } else {
    data$hb_flag[i] = ''
    data$hb_ngml[i] = NA
  }
  
  # Handle dc
  if (this_sample %in% dc$sample) {
    if (all(dc$flag[dc$sample==this_sample] == 'LLQ')) {
      data$dc_flag[i] = 'LLQ'
      data$dc_mgml[i] = min(dc$mgml_av[dc$sample==this_sample])
    } else if (all(dc$flag[dc$sample==this_sample] == 'ULQ')) {
      data$dc_flag[i] = 'ULQ'
      data$dc_mgml[i] = max(dc$mgml_av[dc$sample==this_sample])
    } else {
      data$dc_flag[i] = ''
      data$dc_mgml[i] = mean(dc$mgml_av[dc$sample==this_sample & dc$flag==''])
    } 
  } else {
    data$dc_flag[i] = ''
    data$dc_mgml[i] = NA
  }
  
}

data$prion_category = samples$prion_category[match(data$sample, samples$id)]
data$bftype = samples$bftype[match(data$sample, samples$id)]
data$sample_source = samples$sample_source[match(data$sample, samples$id)]

# make a CSF table for analysis

csf = subset(data, !is.na(bftype) & bftype == 'CSF')

# overwrite data from outside the normal [LLQ, ULQ] range for each assay -- for a handful of samples
# that were at the LLQ or ULQ, we re-did the assay at different dilutions just to see, but if you include
# these observations, then it's not a fair comparison as not all the samples have been treated equally,
# plus you can't include them all on a plot

prp_llq = 1 * 10
prp_ulq = 20 * 50
hb_llq = .82 * 10
hb_ulq = 200 * 100
dc_llq = .1
dc_ulq = 1.6

csf$prp_ngml[csf$prp_ngml < prp_llq] = prp_llq
csf$prp_flag[csf$prp_ngml <= prp_llq] = 'LLQ'
csf$prp_ngml[csf$prp_ngml > prp_ulq] = prp_ulq
csf$prp_flag[csf$prp_ngml >= prp_ulq] = 'ULQ'

csf$hb_ngml[csf$hb_ngml < hb_llq & !is.na(csf$hb_ngml)] = hb_llq
csf$hb_flag[csf$hb_ngml <= hb_llq & !is.na(csf$hb_ngml)] = 'LLQ'
csf$hb_ngml[csf$hb_ngml > hb_ulq & !is.na(csf$hb_ngml)] = hb_ulq
csf$hb_flag[csf$hb_ngml >= hb_ulq & !is.na(csf$hb_ngml)] = 'ULQ'

csf$dc_mgml[csf$dc_mgml < dc_llq & !is.na(csf$dc_mgml)] = dc_llq
csf$dc_flag[csf$dc_mgml <= dc_llq & !is.na(csf$dc_mgml)] = 'LLQ'
csf$dc_mgml[csf$dc_mgml > dc_ulq & !is.na(csf$dc_mgml)] = dc_ulq
csf$dc_flag[csf$dc_mgml >= dc_ulq & !is.na(csf$dc_mgml)] = 'ULQ'

# add some +/- 1.96SE columns for ease of plotting later
csf$prp_l95 = csf$prp_ngml - 1.96*csf$prp_se
csf$prp_u95 = csf$prp_ngml + 1.96*csf$prp_se


# write out the CSF table for quick reference
write.table(format(csf, digits=2, nsmall=2), 'data/processed/csf.tsv',sep='\t',row.names=F,col.names=T,quote=F)

# now bring in the sex and age data which are not being made public
# for privacy reasons we didn't release these variables in the public tables
# so anyone reproducing our analyses at home will see Fig S4D-E as blank
sex_age_available = file.exists('data/samples/sex_age.tsv')
if(sex_age_available) {
  sex_age = read.table('data/samples/sex_age.tsv',sep='\t',header=T)
}


# -----------------------------------------------------------------------------------------------------


# Figure S1. The BetaPrion Human ELISA kit meets FDA technical standards for bioanalytical methods.
# A. Dilution linearity
# B. Spike recovery
# C. Standard curve reproducibility
 

# -----------------------------------------------------------------------------------------------

# Table 1. Some assay performance stats


# I think this is the best figure for within-plate within-dilution technical replicate CV
all_same_plate_summary = sqldf("
                               select   plate, dilution, detail, avg(ngml) mean, stdev(ngml) sd, count(*) n
                               from     elisa_raw
                               where    stype = 'sample' and detail in (select id from samples where bftype = 'CSF')
                               and      ngml / dilution > 1 and ngml / dilution < 20 -- only consider replicates within dynamic range
                               group by 1, 2, 3
                               having   count(*) > 1
                               order by 1, 2, 3
                               ;")
mean(all_same_plate_summary$sd/all_same_plate_summary$mean) # mean within-plate, same-dilution CV

# I think this is the best figure for within-plate technical replicate CV
same_plate_summary = sqldf("
                           select   plate, detail, avg(ngml) mean, stdev(ngml) sd, count(*) n
                           from     elisa_raw
                           where    stype = 'sample' and detail in (select id from samples where bftype = 'CSF')
                           and      ngml / dilution > 1 and ngml / dilution < 20 -- only consider replicates within dynamic range
                           group by 1, 2
                           having   count(*) > 1
                           order by 1, 2
                           ;")
mean(same_plate_summary$sd/same_plate_summary$mean) # mean within-plate, any-dilution CV


# inter-plate CV
ipc = subset(elisa, sample=='v1205.6+C-IPC')

# further dissection of intra- and inter-plate CV in response to Franc+Inga's comments on the manuscript 2018-03-19

elisa_plates = read.table('data/elisa/meta/plates.tsv',sep='\t',header=TRUE)
elisa_plates = elisa_plates[-which(elisa_plates$plate=='6_alt'),] # remove this one for now to avoid much refactoring
elisa_plates$plateno = as.integer(elisa_plates$plate)

ipc$kit_lot = elisa_plates$kit_lot[match(ipc$plate,elisa_plates$plateno)]
ipc$operator = elisa_plates$operator[match(ipc$plate,elisa_plates$plateno)]
elisa_plates$platereader = ''
elisa_plates$platereader[gsub('.*\\.','',elisa_plates$datafile) %in% c('csv','CSV')] = 'Fluostar Optima'
elisa_plates$platereader[gsub('.*\\.','',elisa_plates$datafile) %in% c('txt')] = 'SpectraMax'
ipc$platereader = elisa_plates$platereader[match(ipc$plate,elisa_plates$plateno)]

sd(ipc$prp_ngml) / mean(ipc$prp_ngml)
nrow(ipc)
ipc

# by platereader
sqldf("
select   platereader, avg(prp_ngml) mean, stdev(prp_ngml) sd, stdev(prp_ngml)/avg(prp_ngml) cv, count(*) n
from     ipc
group by 1
order by 1
;")

# by operator
sqldf("
select   operator, avg(prp_ngml) mean, stdev(prp_ngml) sd, stdev(prp_ngml)/avg(prp_ngml) cv, count(*) n
from     ipc
group by 1
order by 1
;")

# by kit lot
sqldf("
select   kit_lot, avg(prp_ngml) mean, stdev(prp_ngml) sd, stdev(prp_ngml)/avg(prp_ngml) cv, count(*) n
from     ipc
group by 1
order by 1
;")

# by all three
sqldf("
      select   kit_lot, operator, platereader, avg(prp_ngml) mean, stdev(prp_ngml) sd, stdev(prp_ngml)/avg(prp_ngml) cv, count(*) n
      from     ipc
      group by 1, 2, 3
      order by 1, 2, 3
      ;")


sqldf("
      select   avg(ten.mean)/avg(fifty.mean) ten_to_fifty, count(*) n
      from     all_same_plate_summary ten, all_same_plate_summary fifty
      where    ten.plate = fifty.plate and ten.detail = fifty.detail
      and      ten.dilution = 10
      and      fifty.dilution = 50
      and      ten.n >= 2
      and      fifty.n >= 2
      ;")



## ratio of LLQ to blank
llq_blank_ratio = sqldf("
                        select   q1.plate, q1.std1_abs, q2.negcon_abs, q1.std1_abs/q2.negcon_abs ratio
                        from (
                        select   plate, avg(a450_620) std1_abs
                        from     elisa_raw
                        where    stype = 'standard' and detail = '1'
                        group by 1
                        order by 1
                        ) q1,
                        (
                        select   plate, avg(a450_620) negcon_abs
                        from     elisa_raw
                        where    stype = 'control' and detail in ('negcon','blocking')
                        group by 1
                        order by 1
                        ) q2
                        where    q1.plate = q2.plate
                        order by 1
                        ;")
llq_blank_ratio
llq_blank_ratio$negcon_abs = pmax(0, llq_blank_ratio$negcon_abs) # get rid of negatives
plates = read.table('data/elisa/meta/plates.tsv', sep='\t', header=T)
plates$fileext = tolower(gsub('.*\\.','',plates$datafile))
llq_blank_ratio$fileext = plates$fileext[match(llq_blank_ratio$plate, plates$plateno)]

sqldf("
select   fileext, median(ratio) median_ratio -- mean is not useful b/c sometimes blank is negative, set to 0 thus ratio is infinite
from     llq_blank_ratio
group by 1
order by 1
;")


# -----------------------------------------------------------------------------------------------

# Figure 1. CSF PrP is highly sensitive to storage and handling conditions.
# A. Storage & handling without CHAPS
# B. Storage & handling with CHAPS
# C. Different materials
# D. PrP - total protein correlation



pdf('figures/manuscript/figure_1.pdf',width=8,height=8)
layout_matrix = matrix(rbind(c(1,1,1,2,2),c(3,3,3,4,4)),nrow=2)
layout(layout_matrix, heights=c(1,1))
par(oma=c(1,2,1,1))

# 1A. Storage & handling without CHAPS
dilution = 50 # 50 for plates 15 & 19 in orig expt, 12 for plates 35 & 36 in replication expt
qc = elisa[elisa$plate==15 & elisa$sample != 'v1205.6+C-IPC',c('sample','prp_ngml')] # 15 is original, 35 is replication
qc$condition = gsub('-.*','',qc$sample)
qc$aliquot = gsub('.+-','',qc$sample)
parms = data.frame(condition=c("control", "large80uL", "small10uL", "RT24h", "FT3x", "Tx1x", "Tx2x", "Tx3x", "tip5x"),
                   y=c(9, 4, 3, 2, 1, 8:6, 5),
                   disp=c("control 40 uL aliquot","larger (80 uL) aliquot","smaller (10 uL) aliquot","24 hours @ RT","freeze/thawed 3X","transferred 1X","transferred 2X","transferred 3X","mixed 10X each with 5 tips"))
qc$y = parms$y[match(qc$condition, parms$condition)]
ylims = c(.5,9.5)
line_breaks = c(8.5, 5.5, 4.5, 2.5, 1.5)

qcs = sqldf("
            select   condition, avg(prp_ngml) mean, stdev(prp_ngml) sd, count(*) n, avg(y) y
            from     qc
            group by 1
            ;")

qcs$se = qcs$sd / sqrt(qcs$n)
qcs$lower95 = qcs$mean - 1.96*qcs$se
qcs$upper95 = qcs$mean + 1.96*qcs$se

# summary stat - how much does 3x transfer reduce PrP?
percent(1 - (mean(qc$prp_ngml[qc$condition=='Tx3x']) / mean(qc$prp_ngml[qc$condition=='control'])))
t.test(qc$prp_ngml[qc$condition=='control'], qc$prp_ngml[qc$condition=='Tx3x'], alternative='two.sided', paired=FALSE)$p.value

par(mar=c(4,12,3,3))
plot(NA, NA, pch=19, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i', xlim=c(0,300), ylim=ylims)
axis(side=1, at=(0:3)*100, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
abline(h=ylims, lwd=2)
abline(h=line_breaks, lwd=.5, col='#777777')
abline(v=0, lwd=2)
abline(v=1*dilution, col='red', lwd=3, lty=2)
mtext(side=1, at=1*dilution, col='red', text='LLQ', line=0.3)
abline(v=20*dilution, col='red', lwd=3, lty=2)
mtext(side=1, at=20*dilution, col='red', text='ULQ', line=0.3)
mtext(side=2, at=parms$y, text=parms$disp, cex=.8, las=2, line=.5)
points(qcs$mean, qcs$y, pch=19)
segments(x0=qcs$lower95, x1=qcs$upper95, y0=qcs$y, y1=qcs$y, col='black', lwd=2)
#mtext(side=3, line=1, font=2, cex=1.2, text='effects of CSF storage and handling', adj=0)

mtext('A', side=3, cex=2, adj = 0.0, line = 0.3)


# QC - storage & handling experiments - with CHAPS
qc_chaps = elisa[elisa$plate==19 & elisa$sample != 'v1205.6+C-IPC',c('sample','prp_ngml')] # 19 is original, 36 is replication
qc_chaps$condition = gsub('-.*','',qc_chaps$sample)
qc_chaps$aliquot = gsub('.+-','',qc_chaps$sample)
parms = data.frame(condition=c("control", "large80uL", "small10uL", "RT24h", "FT3x", "Tx1x", "Tx2x", "Tx3x", "tip5x"),
                   y=c(9, 4, 3, 2, 1, 8:6, 5),
                   disp=c("control 40 uL aliquot","larger (80 uL) aliquot","smaller (10 uL) aliquot","24 hours @ RT","freeze/thawed 3X","transferred 1X","transferred 2X","transferred 3X","mixed 10X each with 5 tips"))
qc_chaps$y = parms$y[match(qc_chaps$condition, parms$condition)]


# summary stat - how much does 3x transfer reduce PrP?
percent(1 - (mean(qc_chaps$prp_ngml[qc$condition=='Tx3x']) / mean(qc_chaps$prp_ngml[qc$condition=='control'])))
t.test(qc_chaps$prp_ngml[qc$condition=='control'], qc_chaps$prp_ngml[qc$condition=='Tx3x'], alternative='two.sided', paired=FALSE)$p.value


qccs = sqldf("
             select   condition, avg(prp_ngml) mean, stdev(prp_ngml) sd, count(*) n, avg(y) y
             from     qc_chaps
             group by 1
             ;")

qccs$se = qccs$sd / sqrt(qccs$n)
qccs$lower95 = qccs$mean - 1.96*qccs$se
qccs$upper95 = qccs$mean + 1.96*qccs$se



par(mar=c(4,0,3,3))
plot(NA, NA, pch=19, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i', xlim=c(0,300), ylim=ylims)
axis(side=1, at=(0:3)*100, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
abline(h=ylims, lwd=2)
abline(h=line_breaks, lwd=.5, col='#777777')
abline(v=0, lwd=2)
abline(v=1*dilution, col='red', lwd=3, lty=2)
mtext(side=1, at=1*dilution, col='red', text='LLQ', line=0.3)
abline(v=20*dilution, col='red', lwd=3, lty=2)
mtext(side=1, at=20*dilution, col='red', text='ULQ', line=0.3)
points(qccs$mean, qccs$y, pch=19)
segments(x0=qccs$lower95, x1=qccs$upper95, y0=qccs$y, y1=qccs$y, col='black', lwd=2)


#mtext(side=2, at=parms$y, text=parms$disp, cex=.8, las=2, line=.5)
#mtext(side=3, line=1, font=2, cex=1.2, text='effects of CSF storage and handling', adj=0)

mtext('B', side=3, cex=2, adj = 0.0, line = 0.3)


materials = elisa[elisa$plate==27,c('sample','prp_ngml')]
materials$aliquot = gsub('.+-','',materials$sample)
materials$txstatus = gsub('.+-','',gsub('-[0-9]$','',materials$sample))
materials$detergent = gsub('v1238.5-','',gsub('-[0-9]x.*','',materials$sample))
materials$condition = paste(materials$detergent, materials$txstatus, sep='-')
parms = data.frame(condition=c("neat-0xTx", "neat-1xPP", "neat-1xPS", "neat-1xPE", "neat-1xGL",
                               "chaps-0xTx", "chaps-1xPP", "chaps-1xPS", "chaps-1xPE", "chaps-1xGL"),
                   y=10:1,
                   disp=rep(c("no transfer","polypropylene","polystyrene","polyethylene","glass"),2)
)
materials$y = parms$y[match(materials$condition, parms$condition)]
ylims = c(.5,10.5)

materials_summarized = sqldf("
                             select   condition, avg(prp_ngml) mean, stdev(prp_ngml) sd, count(*) n, avg(y) y
                             from     materials
                             group by 1
                             ;")

materials_summarized$se = materials_summarized$sd / sqrt(materials_summarized$n)
materials_summarized$lower95 = materials_summarized$mean - 1.96*materials_summarized$se
materials_summarized$upper95 = materials_summarized$mean + 1.96*materials_summarized$se


par(mar=c(4,12,3,3))
plot(materials_summarized$mean, materials_summarized$y, pch=19, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i', xlim=c(0,600), ylim=ylims)
segments(x0=materials_summarized$lower95, x1=materials_summarized$upper95, y0=materials_summarized$y, y1=materials_summarized$y, col='black', lwd=2)
axis(side=1, at=(0:6)*100, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
abline(h=ylims, lwd=2)
abline(v=0, lwd=2)
abline(h=c(5.5), lwd=.5, col='#777777')
abline(v=50, col='red', lwd=3, lty=2)
mtext(side=1, at=50, col='red', text='LLQ', line=0.3)
mtext(side=2, at=parms$y, text=parms$disp, cex=.8, las=2, line=.5)
mtext(side=2, at=c(3,8), text=c('0.03% CHAPS','no detergent'), cex=0.8, font=2, line=8)
#mtext(side=3, line=1, font=2, cex=1.2, text='PrP loss upon 1x transfer to different materials', adj=0)




mtext('C', side=3, cex=2, adj = 0.0, line = 0.3)

par(mar=c(4,4,3,3))
spearman_test = cor.test(csf$prp_ngml, csf$dc_mgml, method='spearman')
n = sum(!is.na(csf$prp_ngml) & !is.na(csf$dc_mgml))
pval = formatC(spearman_test$p.value, digits=2)
rho = formatC(spearman_test$estimate, digits=2)
main_title = 'PrP vs. total protein in CSF'
message = paste('rho = ',rho,', P = ',pval,' (Spearman rank test, N = ',n,')',sep='')
short_message = paste('r = ',rho,', P = ',pval)

m = lm(dc_mgml ~ prp_ngml, data=csf)
intercept = summary(m)$coefficients['(Intercept)','Estimate']
slope = summary(m)$coefficients['prp_ngml','Estimate']

plot(NA, NA, xlim=c(0,600), ylim=c(0,1.7), xlab='', ylab='', axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=(0:6)*100, lwd=0, lwd.ticks=1)
axis(side=2, at=(0:4)/2, las=2, lwd=0, lwd.ticks=1)
abline(h=0, lwd=2)
abline(v=0, lwd=2)
abline(h=c(.1,1.6), lwd=2, lty=2, col='red')
abline(v=c(1*10), lwd=2, lty=2, col='red')
mtext(side=1, at=10, col='red', text='LLQ', line=0.3)
mtext(side=4, at=.1, col='red', text='LLQ', las=2)
mtext(side=4, at=1.6, col='red', text='ULQ', las=2)
abline(a=intercept, b=slope, lwd=0.5, col='#000000')
points(csf$prp_ngml, csf$dc_mgml, pch=19, col = 'black')
#mtext(side=3,line=2,text=main_title,font=2,cex=1.2)
mtext(side=3,line=0.0,text=short_message,cex=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
mtext(side=2, line=3, text='total protein (mg/mL)')


mtext('D', side=3, cex=2, adj = -0.15, line = 0.3)


dev.off()





# -------------------------------------------------------------------------------------------------------

# Figure 2. CSF PrP is brain-derived.
# A. PrP by tissue type
# B. Blood spike-in with CHAPS
# C. PrP vs. hemoglobin



pdf('figures/manuscript/figure_2.pdf',width=8,height=8)
layout_matrix = matrix(rbind(c(1,1),c(2,3)),nrow=2)
layout(layout_matrix, heights=c(1,1))
par(oma=c(1,1,1,1))

# 2A. PrP by tissue type

bftype_params = read.table(textConnection("
bftype|y
blood - plasma|1
blood - buffy coat|2
blood - red cells|3
CSF|4
brain|5
"),header=TRUE,sep='|')

# grab all the non-CSFs from plates 21, 40, 41
prp_by_bftype = elisa[elisa$plate %in% c(21,40,41) & elisa$sample %in% samples$id[samples$bftype != 'CSF'],c('sample','prp_ngml')]
prp_by_bftype = rbind(prp_by_bftype, csf[,c('sample','prp_ngml')])
prp_by_bftype$bftype = samples$bftype[match(prp_by_bftype$sample, samples$id)]
prp_by_bftype$y = bftype_params$y[match(prp_by_bftype$bftype,bftype_params$bftype)]
ylims=c(0.5, 5.5)

par(mar=c(4,8,3,2))
plot(NA,NA,xlim=c(0,2000),ylim=ylims,axes=FALSE,xaxs='i',yaxs='i',ann=FALSE)
axis(side=1, at=(0:4)*500, lwd.ticks=1, lwd=0, las=1)
axis(side=2, at=bftype_params$y, labels=bftype_params$bftype, lwd.ticks=0, lwd=0, las=2)
abline(v=0,lwd=2)
abline(v=(0:4)*500,lwd=.25,col='#777777')
abline(h=ylims,lwd=2)
abline(v=10,lwd=2,lty=2,col='red')
points(prp_by_bftype$prp_ngml, jitter(prp_by_bftype$y,.75), pch=20, col='black')
mtext(side=1,at=10,text='LLQ',col='red')
mtext(side=1,line=2.5,text='[PrP] (ng/mL)')
#title(main='PrP by biofluid type')

mtext('A', side=3, cex=2, adj = 0.0, line = 0.3)


range(prp_by_bftype$prp_ngml[prp_by_bftype$bftype=='blood - buffy coat'])
range(prp_by_bftype$prp_ngml[prp_by_bftype$bftype=='blood - red cells'])
mean(prp_by_bftype$prp_ngml[prp_by_bftype$bftype=='blood - buffy coat'])
mean(prp_by_bftype$prp_ngml[prp_by_bftype$bftype=='blood - red cells'])

# blood spike-in series with CHAPS
bspike = elisa[elisa$plate==20 & grepl('165',elisa$sample),]
bspike$condition = gsub('-.+','',gsub('165\\.2-','',bspike$sample))
parms = data.frame(condition=c('control','1%','0.1%','0.01%','50uM','5uM','500nM'),
                   spiked=c('nothing',rep('blood',3),rep('EDTA',3)),
                   y=c(7,4,5,6,1,2,3),
                   col=c('#000000',rep('#8A0707',3),rep('#C0C0C0',3)),
                   disp=c('control','+1% whole blood', '+0.1% whole blood', '+0.01% whole blood', '+50uM EDTA', '+5uM EDTA', '+500nM EDTA'))
pdat = sqldf("
             select   p.condition, p.y, p.disp, p.col, p.spiked, avg(b.prp_ngml) prp_av, stdev(b.prp_ngml) prp_sd, count(*) n
             from     bspike b, parms p
             where    b.condition = p.condition
             group by 1, 2, 3, 4, 5
             ;")
pdat$l95 = pdat$prp_av - 1.96*pdat$prp_sd / sqrt(pdat$n)
pdat$u95 = pdat$prp_av + 1.96*pdat$prp_sd / sqrt(pdat$n)
ylims = c(.5, 7.5)

par(mar=c(4,8,3,3))
plot(NA, NA, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i', xlim=c(0,300), ylim=ylims)
axis(side=1, at=(0:3)*100, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
abline(h=ylims, lwd=2)
abline(v=0, lwd=2)
abline(v=10, col='red', lwd=3, lty=2)
mtext(side=1, at=10, col='red', text='LLQ')
mtext(side=2, at=parms$y, text=parms$disp, cex=.8, las=2, line=.5)
#mtext(side=3, line=1, font=2, cex=1.2, text='effects of blood spike-in')
segments(x0=pdat$l95, x1=pdat$u95, y0=pdat$y, y1=pdat$y, col='#000000', lwd=2)
points(pdat$prp_av, pdat$y, pch=19, col='#000000')

mtext('B', side=3, cex=2, adj = 0.0, line = 0.3)



spearman_test = cor.test(csf$prp_ngml, csf$hb_ngml, method='spearman')
n = sum(!is.na(csf$prp_ngml) & !is.na(csf$hb_ngml))
pval = formatC(spearman_test$p.value, digits=2)
rho = formatC(spearman_test$estimate, digits=2)
main_title = 'PrP vs. hemoglobin in CSF'
message = paste('rho = ',rho,', P = ',pval,' (Spearman rank test, N = ',n,')',sep='')
short_message = paste('r = ',rho,', P = ',pval)

par(mar=c(4,4,3,3))
plot(NA, NA, xlim=c(0,600), ylim=c(5,50000), log='y', xlab='', ylab='', axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=(0:6)*100, lwd=0, lwd.ticks=1)
axis(side=2, at=10^(1:4), las=2, lwd=0, lwd.ticks=1)
abline(h=5, lwd=2)
abline(v=0, lwd=2)
abline(h=c(hb_llq, hb_ulq), lwd=2, lty=2, col='red')
abline(v=c(prp_llq), lwd=2, lty=2, col='red')
mtext(side=1, at=prp_llq, col='red', text='LLQ')
mtext(side=4, at=hb_llq, col='red', text='LLQ', las=2)
mtext(side=4, at=hb_ulq, col='red', text='ULQ', las=2)
points(csf$prp_ngml, csf$hb_ngml, pch=19, col='black')
#mtext(side=3,line=2,text=main_title,font=2,cex=1.2)
mtext(side=3,line=0.0,text=short_message,cex=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
mtext(side=2, line=3, text='[Hb] (ng/mL)')

mtext('C', side=3, cex=2, adj = 0.0, line = 0.3)

dev.off()



# Figure 3. Test-retest

test_retest = read.table('data/samples/test_retest.tsv',sep='\t',header=T)

pdf('figures/manuscript/figure_3.pdf',width=6,height=3.75)


arnold = subset(test_retest, cohort == 'Arnold')
arnold$x = 1:nrow(arnold)
breaks = which(!duplicated(arnold$indiv_id)) - 1
breaks = c(breaks, nrow(arnold))
mids = sqldf("
             select   indiv_id, avg(x) midx
             from     arnold
             group by 1
             order by 1
             ;")
arnold$id = arnold$sample_id
arnold$prp_ngml = csf$prp_ngml[match(arnold$id, csf$sample)]
arnold$prp_se = csf$prp_se[match(arnold$id, csf$sample)]
arnold$prp_l95 = csf$prp_l95[match(arnold$id, csf$sample)]
arnold$prp_u95 = csf$prp_u95[match(arnold$id, csf$sample)]
arnold$hb_ngml = csf$hb_ngml[match(arnold$id, csf$sample)]
cv_step1 = sqldf("
                 select   indiv_id, avg(prp_ngml) mean_prp, stdev(prp_ngml) sd_prp
                 from     arnold
                 group by indiv_id
                 order by indiv_id
                 ;")
mean_cv = mean(cv_step1$sd_prp/cv_step1$mean_prp)
mean_cv

par(mar=c(4,4,2,1), mfrow=c(1,1))
plot(NA, NA, xlim=c(0,18), ylim=c(0,300), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=mids$midx-.5, labels=1:9, line=0.5, lwd=0, lwd.ticks=0, cex.axis=1, font=2)
axis(side=1, at=arnold$x-.5, labels=arnold$visit, lwd=0, lwd.ticks=0, line=-1, cex.axis=.6)
mtext(side=2, line=0, at=-10, cex=.6, text='visit', las=2)
mtext(side=3, line=0, font=1, cex=1, paste('Mean CV = ',percent(mean_cv),sep=''))
abline(h=0,lwd=2)
abline(v=0,lwd=2)
abline(v=(0:9)*2,lwd=.25)
abline(h=25, lty=2, lwd=2, col='#FF0000')
mtext(side=2, line=.5, at=25, col='#FF0000', text='LLQ', las=2)
segments(x0=as.integer(arnold$id)-.5, x1=as.integer(arnold$id)-.5, y0=arnold$prp_l95, y1=arnold$prp_u95, lwd=2, col='#000000')
points(as.integer(arnold$id)-.5, arnold$prp_ngml, pch=19, col='#000000')
axis(side=2, at=100*(0:5), lwd.ticks=1, lwd=0, las=2, col.axis='#000000', col.ticks='#000000')
mtext(side=2, line=3, text='[PrP] (ng/mL)', col='#000000')
mtext(side=1, line=2.5, text='Individual ID')
dev.off()

# ----------------------------------------------------------------------------------------------------------------
# Figure S1. Technical validation of BetaPrion Human ELISA
# A. Dilution linearity
# C. Standard curve reproducibility


pdf('figures/manuscript/figure_s1.pdf',width=8,height=5)

layout_matrix = matrix(c(1,2,1,3),nrow=2,byrow=T)
layout(layout_matrix, heights=c(1,1))
par(oma=c(1,1,1,1))

# S1A - dilution lineasrity
par(mar=c(4,4,4,3))

dilin = subset(elisa_raw, plate == '11' & detail %in% c('v1187.3-new','165.2'))
dilin$ngml_raw = dilin$ngml / dilin$dilution

plot(NA,NA, xlim=c(0, 1/7), ylim=c(0, 35), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
diparms = data.frame(id=c('v1187.3-new','165.2'), col=c('#901287','#FF9912'), display_name=c('v1187.3','165.2'))
dilin$col = diparms$col[match(dilin$detail, diparms$id)]
points(1/dilin$dilution, dilin$ngml_raw, col=dilin$col, pch=19)
axis(side=1, at=1/c(10,20,50,100), labels=c('1:10','1:20','1:50','1:100'), lwd=0, lwd.ticks=1, cex.axis=.9)
axis(side=2, at=c(0,10,20,30), lwd=0, lwd.ticks=1, las=2)
abline(h=0, lwd=2)
abline(v=0, lwd=2)
mtext(side=1, line=2.5, text='dilution factor')
mtext(side=2, line=3, text='measured [PrP] (ng/mL)')
legend('topleft',diparms$display_name,col=diparms$col,pch=19,text.col=diparms$col,title='CSF sample',title.col='black')
abline(h=c(1,20),lty=2,lwd=3,col='red')
mtext(side=4,at=c(1,20),text=c('LLQ','ULQ'),col='red',las=2)
#mtext(side=3, line=1, cex=1.2, font=2, text='PrP ELISA dilution linearity')


mtext('A', side=3, cex=2, adj = -0.2, line = 0.3)





# standard curve reproducibility plot

stdrep = subset(elisa_raw, plate == '13' & detail %in% c('1','2','4','5','10','20'))
stdrep$prp = as.integer(stdrep$detail)
stdcv = sqldf("
              select   prp, avg(a450_620) mean, stdev(a450_620) sd
              from     stdrep
              group by 1
              order by 1
              ;")
stdcv$cv = stdcv$sd / stdcv$mean

# CV of standard curve points within 1 plate

par(mar=c(0,4,3,3))
plot(NA,NA,xlim=c(0,21),ylim=c(0,2.5),ann=FALSE,axes=FALSE,xaxs='i',yaxs='i')
points(stdrep$prp, stdrep$a450_620, pch=20, cex=.7)
#axis(side=1, at=unique(stdrep$prp), labels=unique(stdrep$prp), cex.axis=.8, lwd=0, lwd.ticks=1)
axis(side=1, at=unique(stdrep$prp), labels=NA, cex.axis=.8, lwd=0, lwd.ticks=1)
abline(h=0)
abline(v=0)
axis(side=2, at=(0:3), labels=(0:3), lwd=0, lwd.ticks=1, las=2)
#mtext(side=1, line=2.5, text='[PrP] (ng/mL) standard curve')
mtext(side=2, line=3, text='A450 - A620')
#mtext(side=3, line=1, text='Reproducibility of BetaPrion Human ELISA standard curve (5X)', font=2, cex=1.2)


mtext('B', side=3, cex=2, adj = -0.2, line = 0.3)

par(mar=c(4,4,0.5,3))
plot(NA,NA,xlim=c(0,21),ylim=c(0,.28),ann=FALSE,axes=FALSE,xaxs='i',yaxs='i')
points(stdcv$prp, stdcv$cv, pch=20)
axis(side=1, at=unique(stdrep$prp), labels=unique(stdrep$prp), cex.axis=.7, lwd=0, lwd.ticks=1, las=2)
abline(h=0)
abline(v=0)
segments(x0=.5, x1=1.5, y0=.25, y1=.25, col='red', lwd=3, lty=2)
segments(x0=1.5, x1=1.5, y0=.25, y1=.20, col='red', lwd=3, lty=2)
segments(x0=1.5, x1=20, y0=.20, y1=.20, col='red', lwd=3, lty=2)
mtext(side=4, at=.20, line=-.5, col='red', las=2, text='FDA limit')
axis(side=2, at=(0:3)/10, labels=percent((0:3)/10), lwd=0, lwd.ticks=1, las=2)
mtext(side=1, line=2.5, text='[PrP] (ng/mL) standard curve')
mtext(side=2, line=3, text='coefficient of variation (CV)')


dev.off() # -- end of Figure S1





# ----------------------------------------------------------------------------------------------------------------
# Figure S2. Plate position effects
# A. Individual loading, raw results
# B. Individual loading, adjusted
# C. Multichannel loading, raw results
# D. Multichannel loading, adjusted

pdf('figures/manuscript/figure_s2.pdf',width=8,height=10)
par(mfrow=c(4,1), mar=c(2,4,1,2), oma=c(4,1,4,1))

# Loading order experiment. Note: Sonia says loading the plate took 29 minutes.
load_order = subset(elisa_raw, plate==22)
load_order$order = 1:nrow(load_order)
load_order$ngml_raw = load_order$ngml / load_order$dilution
parms = data.frame(stype=c('standard','control','sample'),col=c('#9D1309','#444456','#3333FF'))
load_order$col = parms$col[match(load_order$stype, parms$stype)]

m = lm(ngml_raw ~ order, data=subset(load_order, stype=='sample'))
summary(m)
# by what percent does it decline over a whole plate?
summary(m)$coefficients['order','Estimate']*96 / summary(m)$coefficients['(Intercept)','Estimate']
# -29%, though that is from A1 to H12 and sort of interpolated since actual first & last samples are standards/controls
# perhaps more direct to speak of the difference between actual plated samples:
# so what is the average difference between the first and last 10 samples?
first_ten = load_order$ngml[load_order$stype=='sample'][1:10]
last_ten = load_order$ngml[load_order$stype=='sample'][71:80]
mean(last_ten) / mean(first_ten)
# 78%. so the last-loaded samples are 22% lower than the first-loaded.
values = load_order$ngml_raw[load_order$stype=='sample']
# what is the CV?
cv = sd(values) / mean(values)
message = paste('CV: ',percent(cv),sep='')

# -- Begin panel A
plot(load_order$order, load_order$ngml_raw, pch=20, col=load_order$col, axes=FALSE, yaxs='i', ylim=c(0,21), ann=FALSE)
abline(h=0)
axis(side=2, at=c(1,2,4,5,10,20), lwd=0, lwd.ticks=1, las=2)
abline(v=seq(12,96,by=12)+.5,lwd=.25)
mtext(side=2, line=2.5, text='[PrP] (ng/mL)')
mtext(side=3, line=0.5, cex=0.9, text=paste('Single channel loading, no adjustment.',message))

mtext('A', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel A


# can we rescue this problem just analytically, by imputing between the two standard curves?
standards = subset(load_order, stype=='standard')
standards$standard = as.integer(standards$detail)
standards$standard_ratio = standards$standard / standards$ngml_raw
m = lm(ngml_raw ~ order + standard, data=standards)
summary(m)
m = lm(standard_ratio ~ order, data=standards)
summary(m)
order_intercept = summary(m)$coefficients['(Intercept)','Estimate']
order_coeff = summary(m)$coefficients['order','Estimate']
load_order$ngml_adj = load_order$ngml_raw * (order_intercept + order_coeff * load_order$order)

# does this abolish the difference in first and last ten samples?
first_ten = load_order$ngml_adj[load_order$stype=='sample'][1:10]
last_ten = load_order$ngml_adj[load_order$stype=='sample'][71:80]
mean(last_ten) / mean(first_ten)
# 102.5%. so the last-loaded samples are now basically same as first-loaded


values = load_order$ngml_adj[load_order$stype=='sample']
# what is the CV now that it's adjusted?
cv = sd(values) / mean(values)
message = paste('CV: ',percent(cv),sep='')




# -- Begin Panel B
plot(load_order$order, load_order$ngml_adj, pch=20, col=load_order$col, axes=FALSE, yaxs='i', ylim=c(0,21), ann=FALSE)
abline(h=0)
axis(side=2, at=c(1,2,4,5,10,20), lwd=0, lwd.ticks=1, las=2)
abline(v=seq(12,96,by=12)+.5,lwd=.25)
mtext(side=2, line=2.5, text='[PrP] (ng/mL)')
mtext(side=3, line=0.5, cex=0.9, text=paste('Single channel loading, adjusted based on standard curves.',message))

mtext('B', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel B



# Loading order experiment with adjustable-spacing multichannel. Note: Sonia says loading the plate took 10 minutes 50 seconds.
load_order_mc = subset(elisa_raw, plate==23)
load_order_mc$order = 1:nrow(load_order_mc)
load_order_mc$ngml_raw = load_order_mc$ngml / load_order_mc$dilution
parms = data.frame(stype=c('standard','control','sample'),col=c('#9D1309','#444456','#3333FF'))
load_order_mc$col = parms$col[match(load_order_mc$stype, parms$stype)]


m = lm(ngml_raw ~ order, data=subset(load_order_mc, stype=='sample'))
summary(m)
# by what percent does it decline over a whole plate?
summary(m)$coefficients['order','Estimate']*96 / summary(m)$coefficients['(Intercept)','Estimate']
values = load_order_mc$ngml_raw[load_order_mc$stype=='sample']
# what is the CV?
cv = sd(values) / mean(values)
message = paste('CV: ',percent(cv),sep='')

# -- Begin Panel C
plot(load_order_mc$order, load_order_mc$ngml_raw, pch=20, col=load_order_mc$col, axes=FALSE, yaxs='i', ylim=c(0,23), ann=FALSE)
abline(h=0)
axis(side=2, at=c(1,2,4,5,10,20), lwd=0, lwd.ticks=1, las=2)
abline(v=seq(12,96,by=12)+.5,lwd=.25)
mtext(side=2, line=2.5, text='[PrP] (ng/mL)')
mtext(side=3, line=0.5, cex=0.9, text=paste('Multichannel loading, no adjustment.',message))

mtext('C', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel C


# can we rescue this problem just analytically, by imputing between the two standard curves?
standards = subset(load_order_mc, stype=='standard')
standards$standard = as.integer(standards$detail)
standards$standard_ratio = standards$standard / standards$ngml_raw
m = lm(ngml_raw ~ order + standard, data=standards)
summary(m)
m = lm(standard_ratio ~ order, data=standards)
summary(m)
order_intercept = summary(m)$coefficients['(Intercept)','Estimate']
order_coeff = summary(m)$coefficients['order','Estimate']
load_order_mc$ngml_adj = load_order_mc$ngml_raw * (order_intercept + order_coeff * load_order_mc$order)

values = load_order_mc$ngml_adj[load_order_mc$stype=='sample']
# what is the CV now that it's adjusted?
cv = sd(values) / mean(values)
message = paste('CV: ',percent(cv),sep='')


# -- Begin Panel D
plot(load_order_mc$order, load_order_mc$ngml_adj, pch=20, col=load_order_mc$col, axes=FALSE, yaxs='i', ylim=c(0,21), ann=FALSE)
abline(h=0)
axis(side=1, at=seq(12,96,by=12)-6, labels=LETTERS[1:8], lwd=0, lwd.ticks=0, cex.axis=1, line=.5, font=2)
axis(side=1, at=seq(1:96), labels=rep(1:12, 8), lwd=0, lwd.ticks=0, cex.axis=.6, las=2, line=0)
mtext(side=1, line=3, text='plate position (same as loading order)')
axis(side=2, at=c(1,2,4,5,10,20), lwd=0, lwd.ticks=1, las=2)
abline(v=seq(12,96,by=12)+.5,lwd=.25)
mtext(side=2, line=2.5, text='[PrP] (ng/mL)')
mtext(side=3, line=0.5, cex=0.9, text=paste('Multichannel loading, adjusted based on standard curves.',message))


mtext('D', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel D


dev.off() # -- end figure S2













## -----------------------------------------------------------------------------------------------------------------
## Figure S3 - Spike recovery 
# A. spike recovery relative to kit standard curve (plate 13)
# B. in-house recombinant standard curve vs. kit standard curve (plate 30)
# C. spike recovery relative to in-house recombinant standard curve with serial dilution (plate 30)
# D. spike recovery relative to in-house recombinant standard curve minimal handling and blocking buffer (plate 37)
# E. linearity of mixing high and low PrP CSF samples (plate 37)

pdf('figures/manuscript/figure_s3.pdf',width=8,height=8)


layout_matrix = matrix(c(1,2,3,4,5,6),nrow=3,byrow=T)
layout(layout_matrix, heights=c(1,1))
par(oma=c(1,1,1,1))

# A - spike recovery relative to kit standard
par(mar=c(4,4,4,3))

spikerecov = elisa[elisa$plate==28 & elisa$sample != 'v1209.2',c('sample','prp_ngml','prp_se')]
spikerecov$spike_final = as.numeric(gsub('ng/mL','',gsub('.+-','',spikerecov$sample)))
spikerecov$csf = gsub('\\+.*','',spikerecov$sample)
spikeparms = data.frame(id=c('v1205.6','v1202.1'), dilution=c(100,25), col=c('#901287','#FF9912'), display_name=c('high PrP sample - v1205.6 at 1:100','low PrP sample - v1202.1 at 1:25'))
spikerecov$col = spikeparms$col[match(spikerecov$csf, spikeparms$id)]
spikerecov$dilution = spikeparms$dilution[match(spikerecov$csf, spikeparms$id)]
spikerecov$lower95 = spikerecov$prp_ngml - 1.96*spikerecov$prp_se
spikerecov$upper95 = spikerecov$prp_ngml + 1.96*spikerecov$prp_se

par(mar=c(4,4,4,3))

xlims=c(-10,325)
ylims=c(0,2250)
plot(NA,NA,xlim=xlims,ylim=ylims,xaxs='i',yaxs='i',axes=FALSE,ann=FALSE)
mtext(side=1, line=2.5, text='Recombinant PrP spiked in (ng/mL)')
mtext(side=2, line=3, text='Measured [PrP] (ng/mL)')
axis(side=1, at=(0:3)*100, lwd=0, lwd.ticks=1)
abline(h=0)
axis(side=2, at=(0:4)*500, lwd=0, lwd.ticks=1, las=2)
abline(h=0)
points(spikerecov$spike_final * spikerecov$dilution, spikerecov$prp_ngml, col=spikerecov$col, pch=19)
segments(x0=spikerecov$spike_final * spikerecov$dilution, x1=spikerecov$spike_final * spikerecov$dilution, y0=spikerecov$lower95, y1=spikerecov$upper95, col=spikerecov$col, lwd=3)

spikerecov$spike_actual = spikerecov$spike_final * spikerecov$dilution
spikehi = subset(spikerecov, dilution==100)
spikelo = subset(spikerecov, dilution==25)

m = lm(prp_ngml ~ spike_actual, data=spikehi)
summary(m)
hi_intercept = summary(m)$coefficients['(Intercept)','Estimate']
hi_slope = summary(m)$coefficients['spike_actual','Estimate']

abline(a=hi_intercept, b=hi_slope, lwd=1, col=spikehi$col)
abline(a=hi_intercept, b=1, lwd=1, lty=3, col=spikehi$col)

m = lm(prp_ngml ~ spike_actual, data=spikelo)
summary(m)
lo_intercept = summary(m)$coefficients['(Intercept)','Estimate']
lo_slope = summary(m)$coefficients['spike_actual','Estimate']

abline(a=lo_intercept, b=lo_slope, lwd=1, col=spikelo$col)
abline(a=lo_intercept, b=1, lwd=1, lty=3, col=spikelo$col)

legend('topleft',c(spikeparms$display_name,'observed recovery','100% recovery'),col=c(spikeparms$col,'#000000','#000000'),lwd=c(3,3,1,1),lty=c(1,1,1,3),pch=c(20,20,NA,NA),cex=0.9,text.col=c(spikeparms$col,'#000000','#000000'))

mtext('A', side=3, cex=2, adj = -0.2, line = 2.5)


# B - in-house recombinant standard versus kit standard curve
recomb_spike = subset(elisa_raw, plate == 30)
recomb_spike$recomb_added = as.numeric(gsub('.*-','',gsub('ng/mL','',recomb_spike$detail)))
recomb_spike$recomb_added[!grepl('ng/mL',recomb_spike$detail)] = NA

inhouse = subset(recomb_spike, grepl('recomb',detail))
inhouse_col = '#7570B3'
kit = subset(recomb_spike, stype=='standard')
kit_col = '#D95F02'

# fit for kit standard
m = lm(log(as.numeric(detail)) ~ log(a450_620), data=kit)
summary(m)
kitintercept = coefficients(m)[1]
kitslope = coefficients(m)[2]
x = (1:400)/100
f_of_x = exp(intercept + slope*log(x))
points(x = f_of_x, y = x, type='l', lwd=.5, col=kit_col)

# fit for our standard curve
m = lm(log(recomb_added) ~ log(a450_620), data=subset(inhouse, recomb_added < 10 & recomb_added > 0.25))
summary(m)
inhouseintercept = coefficients(m)[1]
inhouseslope = coefficients(m)[2]

# plot the fit of our standard curve
inhouse_standards = unique(inhouse$recomb_added)

plot(NA, NA, xlim=c(0,21), ylim=c(0,4), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=c(1,2,4,5,10,20), lwd.ticks=1, lwd=0, cex.axis=.8)
axis(side=2, at=0:4, lwd.ticks=1, lwd=0, las=2)
#mtext(side=1, line=2.5, text='In-house recombinant standard (ng/mL)')
mtext(side=2, line=2.5, text='Absorbance (arbitrary units)')
abline(h=0,lwd=2)
abline(v=0,lwd=2)

x = (1:400)/100
f_of_x = exp(inhouseintercept + inhouseslope*log(x))
points(x = f_of_x, y = x, type='l', lwd=.5, col=inhouse_col)

# actual points
points(inhouse$recomb_added, inhouse$a450_620, pch=20, col=inhouse_col)

x = (1:400)/100
f_of_x = exp(kitintercept + kitslope*log(x))
points(x = f_of_x, y = x, type='l', lwd=.5, col=kit_col)

points(as.integer(kit$detail), kit$a450_620, pch=20, col=kit_col)

legend('bottomright',c('in-house recombinant','kit standards'),col=c(inhouse_col,kit_col),pch=20,text.col=c(inhouse_col,kit_col))

mtext('B', side=3, cex=2, adj = -0.2, line = 2.5)



# compare in-house curve with spike results
convert_to_ngml = function(absorbance, slope, intercept, dilution_factor=1) {
  estimated_conc = exp(slope*log(absorbance) + intercept) * dilution_factor
  return (estimated_conc)
}
recomb_spike$inhouse_ngml = convert_to_ngml(recomb_spike$a450_620, inhouseslope, inhouseintercept, dilution_factor = recomb_spike$dilution)

spike_ins = subset(recomb_spike, grepl('PrP5',detail))
spike_ins$sample = substr(spike_ins$detail,1,7)
spike_ins$expected = spike_ins$recomb_added
spike_ins$observed = spike_ins$inhouse_ngml / spike_ins$dilution
spikeparms = data.frame(id=c('v1205.6','v1202.1'), dilution=c(100,25), col=c('#901287','#FF9912'), display_name=c('high PrP sample - v1205.6 at 1:100','low PrP sample - v1202.1 at 1:25'))
spike_ins$col = spikeparms$col[match(spike_ins$sample, spikeparms$id)]


plot(spike_ins$expected, spike_ins$observed, pch=20, col = spike_ins$col,
     xlim=c(-0.1,3.5), ylim=c(0,3.5), axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
mtext(side=1, line=2.5, text='Recombinant PrP spiked in (ng/mL) final')
mtext(side=2, line=3, text='Measured [PrP] (ng/mL) final')
axis(side=1, at=(0:3), lwd=0, lwd.ticks=1)
abline(h=0)
axis(side=2, at=(0:4), lwd=0, lwd.ticks=1, las=2)
abline(v=0)

spikelo = subset(spike_ins, sample=='v1202.1')
m = lm(observed ~ expected, data=spikelo)
summary(m)
lo_intercept = summary(m)$coefficients['(Intercept)','Estimate']
lo_slope = summary(m)$coefficients['expected','Estimate']
abline(a=lo_intercept, b=lo_slope, lwd=2, col=spikelo$col)
abline(a=lo_intercept, b=1, lwd=2, lty=3, col=spikelo$col)

spikehi = subset(spike_ins, sample=='v1205.6')
m = lm(observed ~ expected, data=spikehi)
summary(m)
hi_intercept = summary(m)$coefficients['(Intercept)','Estimate']
hi_slope = summary(m)$coefficients['expected','Estimate']
abline(a=hi_intercept, b=hi_slope, lwd=1, col=spikehi$col)
abline(a=hi_intercept, b=1, lwd=1, lty=3, col=spikehi$col)

legend('bottomright',c(spikeparms$display_name,'observed recovery','100% recovery'),col=c(spikeparms$col,'#000000','#000000'),lwd=c(3,3,1,1),lty=c(1,1,1,3),pch=c(20,20,NA,NA),cex=0.9,text.col=c(spikeparms$col,'#000000','#000000'))

mtext('C', side=3, cex=2, adj = -0.2, line = 2.5)


# D - spike recovery with carrier protein - plate 37

recomb_spike = subset(elisa_raw, plate == 37)
recomb_spike$recomb_added = as.numeric(gsub('.*-','',gsub('ng/mL','',recomb_spike$detail)))
recomb_spike$recomb_added[!grepl('ng/mL',recomb_spike$detail)] = NA
inhouse = subset(recomb_spike, grepl('recomb',detail))
inhouse_col = '#7570B3'
# fit for our standard curve
m = lm(log(recomb_added) ~ log(a450_620), data=subset(inhouse, recomb_added < 10 & recomb_added > 0.25))
summary(m)
inhouseintercept = coefficients(m)[1]
inhouseslope = coefficients(m)[2]

# compare in-house curve with spike results
convert_to_ngml = function(absorbance, slope, intercept, dilution_factor=1) {
  estimated_conc = exp(slope*log(absorbance) + intercept) * dilution_factor
  return (estimated_conc)
}
recomb_spike$inhouse_ngml = convert_to_ngml(recomb_spike$a450_620, inhouseslope, inhouseintercept, dilution_factor = recomb_spike$dilution)

spike_ins = subset(recomb_spike, grepl('PrP5',detail))
spike_ins$sample = substr(spike_ins$detail,1,7)
spike_ins$expected = spike_ins$recomb_added
spike_ins$observed = spike_ins$inhouse_ngml / spike_ins$dilution
spikeparms = data.frame(id=c('v1202.1'), dilution=c(25), col=c('#FF9912'), display_name=c('v1202.1 at 1:25'))
spike_ins$col = spikeparms$col[match(spike_ins$sample, spikeparms$id)]


plot(spike_ins$expected, spike_ins$observed, pch=20, col = spike_ins$col,
     xlim=c(-0.1,3.5), ylim=c(0,3.5), axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
mtext(side=1, line=2.5, text='Recombinant PrP spiked in (ng/mL) final')
mtext(side=2, line=3, text='Measured [PrP] (ng/mL) final')
axis(side=1, at=(0:3), lwd=0, lwd.ticks=1)
abline(h=0)
axis(side=2, at=(0:4), lwd=0, lwd.ticks=1, las=2)
abline(v=0)

spikelo = subset(spike_ins, sample=='v1202.1')
m = lm(observed ~ expected, data=spikelo)
summary(m)
lo_intercept = summary(m)$coefficients['(Intercept)','Estimate']
lo_slope = summary(m)$coefficients['expected','Estimate']
abline(a=lo_intercept, b=lo_slope, lwd=2, col=spikelo$col)
abline(a=lo_intercept, b=1, lwd=2, lty=3, col=spikelo$col)

lo_slope

legend('bottomright',c('observed recovery','100% recovery'),col=spikeparms$col[spikeparms$id=='v1202.1'],lwd=c(2,2),lty=c(1,3),cex=1.0)

mtext('D', side=3, cex=2, adj = -0.2, line = 2.5)


# E: mix high and low PrP CSF samples in different proportions


mixhilo_raw = subset(elisa_raw, plate==37 & grepl('%A',detail))
mixhilo_raw$proportion_a = as.numeric(gsub('%A','',mixhilo_raw$detail)) / 100
mixhilo = sqldf("
                select   proportion_a, avg(ngml) mean_ngml, stdev(ngml) sd_ngml, count(*) n
                from     mixhilo_raw
                where    flag = ''
                group by 1
                order by 1
                ;")
mixhilo
mixhilo$l95 = mixhilo$mean_ngml - 1.96 * mixhilo$sd_ngml/sqrt(mixhilo$n)
mixhilo$u95 = mixhilo$mean_ngml + 1.96 * mixhilo$sd_ngml/sqrt(mixhilo$n)

# slope
m = lm(mean_ngml ~ proportion_a, data=mixhilo)
summary(m)
intercept = summary(m)$coefficients['(Intercept)','Estimate']
slope = summary(m)$coefficients['proportion_a','Estimate']


plot(NA, NA, xlim=c(0,1.05), ylim=c(0,600), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=(0:4)/4, labels=percent(0:4/4), lwd=0, lwd.ticks=1)
axis(side=2, at=(0:6)*100, lwd=0, lwd.ticks=1, las=2)
mtext(side=1, line=2.5, text='proportion CSF sample A')
mtext(side=2, line=3, text='[PrP] (ng/mL)')
abline(h=0)
abline(v=0)
abline(a=intercept, b=slope, lwd=1, col=csfcol)
segments(x0=mixhilo$proportion_a, x1=mixhilo$proportion_a, y0=mixhilo$l95, y1=mixhilo$u95, lwd=3, col=csfcol)
points(mixhilo$proportion_a, mixhilo$mean_ngml, pch=19, col=csfcol)

mtext('E', side=3, cex=2, adj = -0.2, line = 2.5)



dev.off() # ----- end Figure S1II







# ----------------------------------------------------------------------------------------------------------------
# Figure S4. 
# PrP levels by diagnosis
# sex, age, etc.
# lumbar-thoracic

pdf('figures/manuscript/figure_s4.pdf',width=8,height=11)

layout_matrix = matrix(c(1,1,2,2,3,3,4,5,6,6,7,7),nrow=6,byrow=T)
layout(layout_matrix, heights=c(1,1))
par(oma=c(1,1,1,1), mar=c(4,4,4,3))

# A - by diagnosis
prion_cohorts = c('Zerr','Parchi')
csf$prion_category = samples$prion_category[match(csf$sample,samples$id)]
parms = data.frame(prion_category=c('non-prion','symptomatic sporadic','symptomatic genetic'),y=c(3,2,1))
top_row = max(parms$y)

diag = subset(csf, sample_source %in% prion_cohorts)
diag$y = parms$y[match(diag$prion_category, parms$prion_category)]

ylims=c(0.5,top_row + 1.5)
xlims=c(0,600)
par(mar=c(4,12,3,3))
plot(NA, NA, xlim=xlims, ylim=ylims, ann=FALSE, axes=FALSE)
axis(side=1, at=(0:6)*100, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
abline(h=ylims, lwd=2)
abline(v=0, lwd=2)
abline(v=10, col='red', lwd=3, lty=2)
mtext(side=1, at=10, col='red', text='LLQ')
mtext(side=2, at=parms$y, text=parms$disp, cex=.8, las=2, line=.5)
axis(side=2, at=parms$y, labels=parms$prion_category, lwd=0, lwd.ticks=0, las=2)
axis(side=2, at=top_row + 1.0, labels='all samples, all cohorts', lwd=0, lwd.ticks=0, las=2)
abline(h=top_row + 0.5, lwd=.5, col='#222222')
par(xpd=TRUE)
segments(x0=0,x1=-200,y0=3.5,y1=3.5,lwd=0.5,col='#222222')
#abline(h=3.5, lwd=.5, col='#222222')
par(xpd=FALSE)
mtext(side=2, at=2, text='surveillance\ncohorts', line=10, font=2, cex=0.6)
points(csf$prp_ngml, jitter(rep(top_row + 1.0, nrow(csf)),.25), pch=20)
points(diag$prp_ngml, jitter(diag$y, .25), pch=20)

sporadic_ratio = mean(diag$prp_ngml[diag$prion_category=='symptomatic sporadic'])/mean(diag$prp_ngml[diag$prion_category=='non-prion'])
sporadic_pval = ks.test(diag$prp_ngml[diag$prion_category=='symptomatic sporadic'],diag$prp_ngml[diag$prion_category=='non-prion'],alternative='two.sided')
genetic_ratio = mean(diag$prp_ngml[diag$prion_category=='symptomatic genetic'])/mean(diag$prp_ngml[diag$prion_category=='non-prion'])
genetic_pval = ks.test(diag$prp_ngml[diag$prion_category=='symptomatic genetic'],diag$prp_ngml[diag$prion_category=='non-prion'],alternative='two.sided')



mtext('A', side=3, cex=2, adj = -0.2, line = 0.5)




# B - by pre/post symptomatic genetic
prion_cohorts = c('Geschwind')
csf$prion_category = samples$prion_category[match(csf$sample,samples$id)]
parms = data.frame(prion_category=c('symptomatic genetic','presymptomatic genetic'),y=c(1,2))

diag = subset(csf, sample_source %in% prion_cohorts)
diag$y = parms$y[match(diag$prion_category, parms$prion_category)]

ylims=c(0.5,2.5)
xlims=c(0,600)
par(mar=c(4,12,3,3))
plot(NA, NA, xlim=xlims, ylim=ylims, ann=FALSE, axes=FALSE)
axis(side=1, at=(0:6)*100, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
abline(h=ylims, lwd=2)
abline(v=0, lwd=2)
abline(v=10, col='red', lwd=3, lty=2)
mtext(side=1, at=10, col='red', text='LLQ')
mtext(side=2, at=parms$y, text=parms$disp, cex=.8, las=2, line=.5)
axis(side=2, at=parms$y, labels=parms$prion_category, lwd=0, lwd.ticks=0, las=2)
par(xpd=TRUE)

#abline(h=3.5, lwd=.5, col='#222222')
par(xpd=FALSE)
mtext(side=2, at=1.5, text='UCSF cohort', line=10, font=2, cex=0.8)
points(diag$prp_ngml, jitter(diag$y, .25), pch=20)

gsymp_ratio = mean(diag$prp_ngml[diag$prion_category=='symptomatic genetic'])/mean(diag$prp_ngml[diag$prion_category=='presymptomatic genetic'])
gsymp_pval = ks.test(diag$prp_ngml[diag$prion_category=='symptomatic genetic'],diag$prp_ngml[diag$prion_category=='presymptomatic genetic'],alternative='two.sided')

# we don't have properly controlled within-cohort comparison for presymptomatic genetic vs. non-prion individuals,
# but just to get an idea of whether there might be a difference:
# ks.test(diag$prp_ngml[diag$prion_category=='presymptomatic genetic'], csf$prp_ngml[csf$prion_category=='non-prion'])
# t.test(diag$prp_ngml[diag$prion_category=='presymptomatic genetic'], csf$prp_ngml[csf$prion_category=='non-prion'])
# mean(diag$prp_ngml[diag$prion_category=='presymptomatic genetic'])
# mean(csf$prp_ngml[csf$prion_category=='non-prion'])


mtext('B', side=3, cex=2, adj = -0.2, line = 0.5)




# reset margins to something normal
par(mar=c(4,5,3,3))

# C


parms = data.frame(cohort=unique(csf$sample_source), x=1:length(unique(csf$sample_source)))
parms$display = c('MIND\n','metformin\ntrial','cognitive\nimpairment','Bologna\n','sapropterin\ntrial','Gottingen\n','UCSF\n')
csf$x = parms$x[match(csf$sample_source, parms$cohort)]
xlims=c(0.5, max(parms$x)+0.5)
ylims=c(0,650)
plot(NA, NA, xlim=xlims, ylim=ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=parms$x, labels=parms$display, lwd.ticks=0, lwd=0, cex.axis=1, las=1)
axis(side=2, at=(0:6)*100, lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=3, text='CSF [PrP] (ng/mL)')
abline(h=0,lwd=2)
abline(v=xlims,lwd=2)
abline(h=10, lty=2, lwd=2, col='red')
mtext(side=4, line=.5, at=10, col='red', text='LLQ', las=2)
subs = !(csf$prion_category %in% c('symptomatic genetic','symptomatic sporadic'))
points(jitter(csf$x[subs], .25), csf$prp_ngml[subs], pch=20)

anova_obj = aov(prp_ngml ~ sample_source, data=subset(csf, !(prion_category %in% c('symptomatic genetic','symptomatic sporadic'))))
anova_pval = summary(anova_obj)[[1]][["Pr(>F)"]][1]

mtext('C', side=3, cex=2, adj = -0.05, line = 1)


## response to a query from Franc & Inga - 2018-03-19
csf$prion_category = samples$prion_category[match(csf$sample,samples$id)]
anova_obj = aov(prp_ngml ~ sample_source + prion_category, data=subset(csf, sample_source %in% c('Parchi','Zerr')))
anova_pval = summary(anova_obj)[[1]][["Pr(>F)"]][1]
summary(anova_obj)
anova_obj = aov(prp_ngml ~ sample_source, data=subset(csf, sample_source %in% c('Parchi','Zerr','Geschwind') & prion_category == 'symptomatic genetic' ))
summary(anova_obj)
m = lm(prp_ngml ~ prion_category, data=subset(csf, sample_source %in% c('Parchi','Zerr','Geschwind')))
summary(m)


if (sex_age_available) {
  
  csf$age = sex_age$age[match(csf$sample,sex_age$id)]
  csf$sex = sex_age$sex[match(csf$sample,sex_age$id)]
  
  age_spearman = cor.test(csf$age, csf$prp_ngml, method='spearman', alternative='two.sided')
  age_rho = age_spearman$estimate
  age_pval = age_spearman$p.value
  
  m = lm(prp_ngml ~ sample_source, data=csf)
  summary(m)
  m = lm(prp_ngml ~ sample_source + age, data=csf)
  summary(m)
  m = lm(prp_ngml ~ age, data=csf)
  summary(m)
  m = lm(prp_ngml ~ sample_source + prion_category + age, data=csf)
  summary(m)
  m = lm(age ~ sample_source, data=csf)
  summary(m)
  
  # example of how age is confounded - for Figure S4 legend
  ks.test(csf$age[csf$sample_source %in% c('Parchi', 'Zerr') & csf$prion_category=='symptomatic genetic'], csf$age[csf$sample_source %in% c('Parchi', 'Zerr') & csf$prion_category=='symptomatic sporadic'], alternative='two.sided')
  mean(csf$age[csf$sample_source %in% c('Parchi', 'Zerr') & csf$prion_category=='symptomatic genetic'])
  mean(csf$age[csf$sample_source %in% c('Parchi', 'Zerr') & csf$prion_category=='symptomatic sporadic'])
  m = lm(prp_ngml ~ age, data=subset(csf, sample_source %in% c('Parchi', 'Zerr')))
  summary(m)
  m = lm(prp_ngml ~ prion_category + age, data=subset(csf, sample_source %in% c('Parchi', 'Zerr')))
  summary(m) 
}

xlims=c(0,100)
ylims=c(0,600)
plot(NA, NA, xlim=xlims, ylim=ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=2, at=(0:6)*100, lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=3, text='CSF [PrP] (ng/mL)')
abline(h=0,lwd=2)
abline(v=0,lwd=2)
abline(h=10, lty=2, lwd=2, col='red')
mtext(side=4, line=.5, at=10, col='red', text='LLQ', las=2)
axis(side=1, at=(1:10)*10, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='age')
if (sex_age_available) {
  points(csf$age, csf$prp_ngml, pch=20)
}
m = lm(prp_ngml ~ age, data=csf)
intercept = summary(m)$coefficients['(Intercept)','Estimate']
slope = summary(m)$coefficients['age','Estimate']
abline(a=intercept, b=slope, lwd=0.5, col='red')

mtext('D', side=3, cex=2, adj = -0.1, line = 1)




xlims=c(0.5,2.5)
ylims=c(0,600)


plot(NA, NA, xlim=xlims, ylim=ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=2, at=(0:6)*100, lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=3, text='CSF [PrP] (ng/mL)')
abline(h=0,lwd=2)
abline(v=xlims,lwd=2)
abline(h=10, lty=2, lwd=2, col='red')
mtext(side=4, line=.5, at=10, col='red', text='LLQ', las=2)
axis(side=1, at=c(1,2), labels=c('male','female'), lwd=0, lwd.ticks=1)
subs = !(csf$prion_category %in% c('symptomatic genetic','symptomatic sporadic'))
if (sex_age_available) {
  csf$x_sex = as.numeric(NA)
  csf$x_sex[csf$sex=='m'] = 1
  csf$x_sex[csf$sex=='f'] = 2
  points(jitter(csf$x_sex[subs],0.25), csf$prp_ngml[subs], pch=20)
  sex_ks_p = ks.test(csf$prp_ngml[csf$sex=='m'], csf$prp_ngml[csf$sex=='f'], alternative='two.sided')
  
  subs = !(csf$prion_category %in% c('symptomatic genetic','symptomatic sporadic'))
  sex_ks_p = ks.test(csf$prp_ngml[csf$sex=='m' & subs], csf$prp_ngml[csf$sex=='f' & subs], alternative='two.sided')
  
}



mtext('E', side=3, cex=2, adj = -0.1, line = 1)



lumtho_drains = subset(elisa, plate==31 & sample != 'v1205.6+C-IPC')
lumtho_drains$u95 = lumtho_drains$prp_ngml + 1.96 * lumtho_drains$prp_se
lumtho_drains$l95 = lumtho_drains$prp_ngml - 1.96 * lumtho_drains$prp_se
lumtho_drains$indiv = substr(lumtho_drains$sample,1,5)
lumtho_drains$syringe = as.integer(substr(lumtho_drains$sample,7,7))
lumtho_drains$x = (as.integer(as.factor(lumtho_drains$indiv))-1)*4 + as.integer(lumtho_drains$syringe)

indivs = unique(lumtho_drains$indiv)
breaks = (0:3)*4 + 0.5
mids = (1:3)*4 - 1.5


plot(NA, NA, xlim=c(0.5, 12.5), ylim=c(0, 210), axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
abline(h=0, lwd=2)
abline(v=breaks, lwd=1)
abline(v=range(breaks), lwd=2) # double thick on outer borders because half will be cut off by xaxs='i'

axis(side=1, at=lumtho_drains$x, labels=lumtho_drains$syringe, lwd=0, lwd.ticks=0, line=0)
axis(side=1, at=mids, labels=indivs, lwd=0, lwd.ticks=0, line=1)
axis(side=2, at=(0:4)*50, labels=(0:4)*50, lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='CSF [PrP] (ng/mL)')
segments(x0=lumtho_drains$x, x1=lumtho_drains$x, y0=lumtho_drains$l95, y1=lumtho_drains$u95, lwd=2, col=csfcol)
points(x=lumtho_drains$x, y=lumtho_drains$prp_ngml, pch=19, col=csfcol)


lumtho_drain_rel = sqldf("
select   l1.indiv, lother.syringe, lother.prp_ngml/l1.prp_ngml rel_prp
from     lumtho_drains l1, lumtho_drains lother
where    l1.indiv = lother.indiv
and      l1.syringe = 1
order by 1, 2
;")
m = lm(rel_prp ~ syringe, data=lumtho_drain_rel)
summary(m)


mtext('F', side=3, cex=2, adj = -0.05, line = 1)


lumtho_syringe = subset(elisa, plate==41 & grepl('^v',sample) & sample != 'v1205.6+C-IPC')
lumtho_syringe$u95 = lumtho_syringe$prp_ngml + 1.96 * lumtho_syringe$prp_se
lumtho_syringe$l95 = lumtho_syringe$prp_ngml - 1.96 * lumtho_syringe$prp_se
lumtho_syringe$indiv = substr(lumtho_syringe$sample,1,5)
lumtho_syringe$syringe = substr(lumtho_syringe$sample,7,7)
lumtho_syringe$rel_syringe = rep(1:4,5)
lumtho_syringe$x = 1:20 #(as.integer(as.factor(lumtho_syringe$indiv))-1)*4 + as.integer(lumtho_syringe$syringe)

indivs = unique(lumtho_syringe$indiv)
breaks = (0:5)*4 + 0.5
mids = (1:5)*4 - 1.5

plot(NA, NA, xlim=c(0.5, 20.5), ylim=c(0, 500), axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
abline(h=0, lwd=2)
abline(v=breaks, lwd=1)
abline(v=range(breaks), lwd=2) # double thick on outer borders because half will be cut off by xaxs='i'

axis(side=1, at=lumtho_syringe$x, labels=lumtho_syringe$rel_syringe, lwd=0, lwd.ticks=0, line=0)
axis(side=1, at=mids, labels=indivs, lwd=0, lwd.ticks=0, line=1)
axis(side=2, at=(0:5)*100, labels=(0:5)*100, lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='CSF [PrP] (ng/mL)')
segments(x0=lumtho_syringe$x, x1=lumtho_syringe$x, y0=lumtho_syringe$l95, y1=lumtho_syringe$u95, lwd=2, col=csfcol)
points(x=lumtho_syringe$x, y=lumtho_syringe$prp_ngml, pch=19, col=csfcol)


mtext('G', side=3, cex=2, adj = -0.05, line = 1)


lumtho_syr_rel = sqldf("
select   l1.indiv, lother.rel_syringe, lother.prp_ngml/l1.prp_ngml rel_prp
from     lumtho_syringe l1, lumtho_syringe lother
where    l1.indiv = lother.indiv
and      l1.rel_syringe = 1
order by 1, 2
;")
# plot(NA, NA, xlim=c(0.5, 4.5), ylim=c(.5,1.5), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
# axis(side=1, at=1:4, lwd=0, lwd.ticks=1)
# abline(h=1)
# axis(side=2, at=(5:15)/10, labels=percent((5:15)/10), las=2, lwd=0, lwd.ticks=1)
# for (indiv in unique(lumtho_syr_rel$indiv)) {
#   subs = lumtho_syr_rel$indiv == indiv
#   points(lumtho_syr_rel$rel_syringe[subs], lumtho_syr_rel$rel_prp[subs], type='l', lwd=3, col='red')
# }
m = lm(rel_prp ~ rel_syringe, data=lumtho_syr_rel)
summary(m)


dev.off()

# possible different way to plot this

# info for Figure S4 legend
sporadic_ratio
sporadic_pval
genetic_ratio
genetic_pval
gsymp_ratio
gsymp_pval
anova_pval


max(csf$prp_ngml) / min(csf$prp_ngml)
sd(csf$prp_ngml) / mean(csf$prp_ngml)


if(sex_age_available) {
  age_rho
  age_pval
  sex_ks_p
  m = lm(prp_ngml ~ sample_source + prion_category + age, data=csf)
  summary(m)
  
}




cohort_smry = sqldf("
select   sample_source cohort, stdev(prp_ngml) sd, avg(prp_ngml) mean, count(*) n, max(prp_ngml) highest, min(prp_ngml) lowest
from     csf
where    prion_category not in ('symptomatic sporadic', 'symptomatic genetic')
group by 1
order by 1
;")

cohort_smry$cv = cohort_smry$sd / cohort_smry$mean
cohort_smry
mean(cohort_smry$cv)
mean(cohort_smry$highest / cohort_smry$lowest)


diag_smry = sqldf("
select   prion_category diagnosis, stdev(prp_ngml) sd, avg(prp_ngml) mean, count(*) n, max(prp_ngml) highest, min(prp_ngml) lowest
from     csf
where    sample_source in ('Zerr','Parchi')
group by 1
order by 1
;")
diag_smry

diag_smry$cv = diag_smry$sd / diag_smry$mean
diag_smry
mean(diag_smry$cv)
mean(diag_smry$highest / diag_smry$lowest)


sum(csf$sample_source=='Zerr')
sum(samples$sample_source=='Zerr')

# ----------------------------------------------------------------------------------------------------------------
# Figure S5. Additional evidence regarding loss of PrP to plastic adsorption
# A. Effect of aliquot size
# B. Effect of syringe volume
# C. Effect of post-hoc addition of CHAPS
# D. Effect of adding BSA


pdf('figures/manuscript/figure_s5.pdf',width=8,height=9)
#layout_matrix = matrix(c(1,2,2,1,2,2,3,3,3),nrow=3,byrow=T)
#layout_matrix = matrix(c(1,1,2,2,2,2,3,3,3,4,4,4),nrow=2,byrow=T)
#layout_matrix = matrix(c(1,2,3,4),nrow=2,byrow=T)
#layout_matrix = matrix(c(1,2,3,3,4,4),nrow=3,byrow=T)
#layout_matrix = matrix(c(1,3,3,2,3,3,4,4,4),nrow=3,byrow=T)
layout_matrix = matrix(c(1,1,3,3,3,3,2,2,3,3,3,3,4,4,4,4,5,5,4,4,4,4,5,5,6,6,6,6,6,6),nrow=5,byrow=T)
layout(layout_matrix, heights=c(1,1))
par(oma=c(1,2,2,1))

# 


# aliquot size
aqsize = subset(elisa_raw, plate==14 & grepl('aliquot',elisa_raw$detail))
aqsize$size = as.integer(gsub('.+-','',gsub('uL.+','',aqsize$detail)))
aqsize$number = as.integer(gsub('.+aliquot','',aqsize$detail))
aq = sqldf("
           select   size, number, avg(ngml) mean, stdev(ngml) sd, count(*) n
           from     aqsize
           group by 1, 2
           order by 1 desc, 2 asc
           ;")
aq$se = aq$sd / sqrt(aq$n)
aq$x = c(1:4,6:9)

# -- Begin Panel A
par(mar=c(4,4,4,2))
plot(NA, NA, xlim=c(0,10), ylim=c(0,425), xaxs='i', yaxs='i', axes=FALSE, ann=FALSE)
axis(side=1, at=aq$x, labels=aq$number, line=-1, cex.axis=.8, lwd.ticks=0, lwd=0)
axis(side=1, at=c(2.5,7.5), labels=c('40uL','10uL'), line=0, cex.axis=1, font=2, lwd.ticks=0, lwd=0)
abline(h=0, lwd=2)
abline(v=0, lwd=2)
abline(v=c(0,5,10), lwd=1, col='#777777')
axis(side=2, at=(0:4)*100, labels=(0:4)*100, lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=3, text="[PrP] (ng/mL)")
points(aq$x, aq$mean, pch=19, cex=1.2, col='black')
segments(x0=aq$x, x1=aq$x, y0=aq$mean-1.96*aq$se, y1=aq$mean+1.96*aq$se, lwd=2, col='black')
mtext(side=1, line=2, text='aliquots of same sample')

mtext('A', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel A



# --- begin panel B - syringe vol

syringevol = subset(elisa, plate==33 & grepl('pool',sample))
syringevol$u95 = syringevol$prp_ngml + 1.96 * syringevol$prp_se
syringevol$l95 = syringevol$prp_ngml - 1.96 * syringevol$prp_se
syringevol$vol = gsub('pool-','',syringevol$sample)
syringevol$vol = gsub('mL',' mL',gsub('std','STD',syringevol$vol))
syringevol$x = 1:4

plot(NA, NA, xlim=c(0.5, 4.5), ylim=c(0, 300), axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
abline(h=0, lwd=2)
abline(v=c(0.5,4.5), lwd=2)

axis(side=1, at=syringevol$x, labels=syringevol$vol, lwd=0, lwd.ticks=0, line=0)
axis(side=2, at=(0:3)*100, labels=(0:3)*100, lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=3.0, text='CSF [PrP] (ng/mL)')
mtext(side=1, line=2.5, text='syringe volume')
segments(x0=syringevol$x, x1=syringevol$x, y0=syringevol$l95, y1=syringevol$u95, lwd=2, col='#000000')
points(x=syringevol$x, y=syringevol$prp_ngml, pch=19, col='#000000')

mtext('B', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel B


# --- begin panel C -- post-hoc chaps addition
# post-hoc CHAPS addition
postchaps = elisa[elisa$plate==21 & grepl('-without|-chaps-post',elisa$sample),]
postchaps$condition = gsub('.+-','',postchaps$sample)
postchaps$id = gsub('-.*','',postchaps$sample)
parms = data.frame(id=rep(unique(postchaps$id),2),
                   condition=rep(unique(postchaps$condition),each=nrow(postchaps)/2),
                   col=rep(c('#000000','#0001CD'),each=nrow(postchaps)/2),
                   y=c((nrow(postchaps)/2):1 - .33, (nrow(postchaps)/2):1 - .67))
pdat = sqldf("
             select   c.id, c.condition, p.y, p.col, c.prp_ngml, c.prp_se
             from     postchaps c, parms p
             where    c.condition = p.condition and c.id = p.id
             ;")
pdat$l95 = pdat$prp_ngml - 1.96*pdat$prp_se
pdat$u95 = pdat$prp_ngml + 1.96*pdat$prp_se
ylims = c(0, nrow(postchaps)/2)

without = postchaps$prp_ngml[postchaps$condition=='without']
post = postchaps$prp_ngml[postchaps$condition=='post']

mean_increase = mean(post/without) - 1 
pval = t.test(without, post, paired=TRUE, alternative='two.sided')$p.value
message = paste('mean ',percent(mean_increase),' increase with CHAPS, P = ',formatC(pval,format='fg',digits=2),sep='')


par(mar=c(4,6,3,3))
plot(NA, NA, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i', xlim=c(0,400), ylim=ylims)
axis(side=1, at=(0:4)*100, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
axis(side=2, at=pdat$y[postchaps$condition=='without']-.17, labels=pdat$id[postchaps$condition=='without'], las=2, lwd=0, lwd.ticks=0)
abline(h=min(ylims):max(ylims), lwd=.5, col='#222222')
abline(h=ylims, lwd=2)
abline(v=0, lwd=2)
abline(v=10, col='red', lwd=3, lty=2)
mtext(side=1, at=10, col='red', text='LLQ')
mtext(side=2, at=parms$y, text=parms$disp, cex=.8, las=2, line=.5)
mtext(side=3, line=0, cex=.7, text=message)
segments(x0=pdat$l95, x1=pdat$u95, y0=pdat$y, y1=pdat$y, col=pdat$col, lwd=2)
points(pdat$prp_ngml, pdat$y, pch=19, col=pdat$col)
legend('bottomright',c('no CHAPS','0.03% CHAPS added'),col=c('#000000','#0001CD'),lwd=2,pch=19,bg='white',box.lwd=2,box.col='#222222')

mtext('C', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel B

## Begin panels D & E - replication experiment of Figure 1 A-B


# 1A. Storage & handling without CHAPS
dilution = 12 # 50 for plates 15 & 19 in orig expt, 12 for plates 35 & 36 in replication expt
qc = elisa[elisa$plate==35 & elisa$sample != 'v1205.6+C-IPC',c('sample','prp_ngml')] # 15 is original, 35 is replication
qc$condition = gsub('-.*','',qc$sample)
qc$aliquot = gsub('.+-','',qc$sample)
parms = data.frame(condition=c("control", "large80uL", "small10uL", "RT24h", "FT3x", "Tx1x", "Tx2x", "Tx3x", "tip5x"),
                   y=c(9, 4, 3, 2, 1, 8:6, 5),
                   disp=c("control 40 uL aliquot","larger (80 uL) aliquot","smaller (10 uL) aliquot","24 hours @ RT","freeze/thawed 3X","transferred 1X","transferred 2X","transferred 3X","mixed 10X each with 5 tips"))
qc$y = parms$y[match(qc$condition, parms$condition)]
ylims = c(.5,9.5)
line_breaks = c(8.5, 5.5, 4.5, 2.5, 1.5)

qcs = sqldf("
            select   condition, avg(prp_ngml) mean, stdev(prp_ngml) sd, count(*) n, avg(y) y
            from     qc
            group by 1
            ;")

qcs$se = qcs$sd / sqrt(qcs$n)
qcs$lower95 = qcs$mean - 1.96*qcs$se
qcs$upper95 = qcs$mean + 1.96*qcs$se

# summary stat - how much does 3x transfer reduce PrP?
percent(1 - (mean(qc$prp_ngml[qc$condition=='Tx3x']) / mean(qc$prp_ngml[qc$condition=='control'])))
t.test(qc$prp_ngml[qc$condition=='control'], qc$prp_ngml[qc$condition=='Tx3x'], alternative='two.sided', paired=FALSE)$p.value

par(mar=c(4,20,3,3))
plot(NA, NA, pch=19, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i', xlim=c(0,300), ylim=ylims)
axis(side=1, at=(0:3)*100, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
abline(h=ylims, lwd=2)
abline(h=line_breaks, lwd=.5, col='#777777')
abline(v=0, lwd=2)
abline(v=1*dilution, col='red', lwd=3, lty=2)
mtext(side=1, at=1*dilution, col='red', text='LLQ', line=0.3)
abline(v=20*dilution, col='red', lwd=3, lty=2)
mtext(side=1, at=20*dilution, col='red', text='ULQ', line=0.3)
mtext(side=2, at=parms$y, text=parms$disp, cex=.8, las=2, line=.5)
points(qcs$mean, qcs$y, pch=19)
segments(x0=qcs$lower95, x1=qcs$upper95, y0=qcs$y, y1=qcs$y, col='black', lwd=2)
#mtext(side=3, line=1, font=2, cex=1.2, text='effects of CSF storage and handling', adj=0)

mtext('D', side=3, cex=2, adj = 0.0, line = 0.3)


# QC - storage & handling experiments - with CHAPS
qc_chaps = elisa[elisa$plate==36 & elisa$sample != 'v1205.6+C-IPC',c('sample','prp_ngml')] # 19 is original, 36 is replication
qc_chaps$condition = gsub('-.*','',qc_chaps$sample)
qc_chaps$aliquot = gsub('.+-','',qc_chaps$sample)
parms = data.frame(condition=c("control", "large80uL", "small10uL", "RT24h", "FT3x", "Tx1x", "Tx2x", "Tx3x", "tip5x"),
                   y=c(9, 4, 3, 2, 1, 8:6, 5),
                   disp=c("control 40 uL aliquot","larger (80 uL) aliquot","smaller (10 uL) aliquot","24 hours @ RT","freeze/thawed 3X","transferred 1X","transferred 2X","transferred 3X","mixed 10X each with 5 tips"))
qc_chaps$y = parms$y[match(qc_chaps$condition, parms$condition)]


# summary stat - how much does 3x transfer reduce PrP?
percent(1 - (mean(qc_chaps$prp_ngml[qc$condition=='Tx3x']) / mean(qc_chaps$prp_ngml[qc$condition=='control'])))
t.test(qc_chaps$prp_ngml[qc$condition=='control'], qc_chaps$prp_ngml[qc$condition=='Tx3x'], alternative='two.sided', paired=FALSE)$p.value


qccs = sqldf("
             select   condition, avg(prp_ngml) mean, stdev(prp_ngml) sd, count(*) n, avg(y) y
             from     qc_chaps
             group by 1
             ;")

qccs$se = qccs$sd / sqrt(qccs$n)
qccs$lower95 = qccs$mean - 1.96*qccs$se
qccs$upper95 = qccs$mean + 1.96*qccs$se



par(mar=c(4,0,3,3))
plot(NA, NA, pch=19, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i', xlim=c(0,300), ylim=ylims)
axis(side=1, at=(0:3)*100, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
abline(h=ylims, lwd=2)
abline(h=line_breaks, lwd=.5, col='#777777')
abline(v=0, lwd=2)
abline(v=1*dilution, col='red', lwd=3, lty=2)
mtext(side=1, at=1*dilution, col='red', text='LLQ', line=0.3)
abline(v=20*dilution, col='red', lwd=3, lty=2)
mtext(side=1, at=20*dilution, col='red', text='ULQ', line=0.3)
points(qccs$mean, qccs$y, pch=19)
segments(x0=qccs$lower95, x1=qccs$upper95, y0=qccs$y, y1=qccs$y, col='black', lwd=2)


#mtext(side=2, at=parms$y, text=parms$disp, cex=.8, las=2, line=.5)
#mtext(side=3, line=1, font=2, cex=1.2, text='effects of CSF storage and handling', adj=0)

mtext('E', side=3, cex=2, adj = 0.0, line = 0.3)




# -- Begin Panel F - BSA expt

bsa_expt = elisa[elisa$plate==32 & grepl('165.2',elisa$sample),]
bsa_expt$additive = substr(bsa_expt$sample,7,9)
bsa_expt$transfers = as.integer(gsub('x-.*','',gsub('.*Tx','',bsa_expt$sample)))
bsa_expt$vial_replicate = as.integer(gsub('.*-','',bsa_expt$sample))

bsa_smry = sqldf("
select   additive, transfers, avg(prp_ngml) mean, stdev(prp_ngml) sd, count(*) n
from     bsa_expt
group by 1, 2
order by 1, 2
;")
bsa_smry$y = 4 - bsa_smry$transfers
bsa_smry$col = ''
bsa_smry$col[bsa_smry$additive=='PBS'] = '#000000'
bsa_smry$col[bsa_smry$additive=='BSA'] = '#0001CD'
bsa_smry$l95 = bsa_smry$mean - 1.96 * bsa_smry$sd / sqrt(bsa_smry$n)
bsa_smry$u95 = bsa_smry$mean + 1.96 * bsa_smry$sd / sqrt(bsa_smry$n)
ylims=c(0.5, 4.5)
dilution = 20

par(mar=c(4,8,3,3))
plot(NA, NA, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i', xlim=c(0,100), ylim=ylims)
axis(side=1, at=(0:5)*20, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='[PrP] (ng/mL)')
axis(side=2, at=bsa_smry$y[bsa_smry$additive=='BSA'], labels=paste("transferred ",bsa_smry$transfers[bsa_smry$additive=='BSA'],"x",sep=""), las=2, lwd=0, lwd.ticks=0)
abline(h=min(ylims):max(ylims), lwd=.5, col='#222222')
abline(h=ylims, lwd=2)
abline(v=0, lwd=2)
abline(v=dilution, col='red', lwd=3, lty=2)
mtext(side=1, at=dilution, col='red', text='LLQ')
# mtext(side=2, at=parms$y, text=parms$disp, cex=.8, las=2, line=.5)

segments(x0=bsa_smry$l95, x1=bsa_smry$u95, y0=bsa_smry$y, y1=bsa_smry$y, col=bsa_smry$col, lwd=2)
points(bsa_smry$mean, bsa_smry$y, pch=19, col=bsa_smry$col)
legend('bottomright',paste("+",unique(bsa_smry$additive),sep=""),col=unique(bsa_smry$col),text.col=unique(bsa_smry$col),pch=19,bg='white',box.lwd=2,box.col='#222222')

mtext('F', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel B



dev.off() # -- end Figure S5





# ----------------------------------------------------------------------------------------------------------------
# Figure S6. CSF hemoglobin variation across test-retest samples

pdf('figures/manuscript/figure_s6.pdf',width=8,height=5)


arnold = subset(test_retest, cohort=='Arnold')
arnold$id = arnold$sample_id
arnold$x = 1:nrow(arnold)
breaks = which(!duplicated(arnold$indiv_id)) - 1
breaks = c(breaks, nrow(arnold))
mids = sqldf("
             select   indiv_id, avg(x) midx
             from     arnold
             group by 1
             order by 1
             ;")
arnold$prp_ngml = csf$prp_ngml[match(arnold$id, csf$sample)]
arnold$prp_se = csf$prp_se[match(arnold$id, csf$sample)]
arnold$prp_l95 = csf$prp_l95[match(arnold$id, csf$sample)]
arnold$prp_u95 = csf$prp_u95[match(arnold$id, csf$sample)]
arnold$hb_ngml = csf$hb_ngml[match(arnold$id, csf$sample)]
prp_cv_step1 = sqldf("
                 select   indiv_id, avg(prp_ngml) mean_prp, stdev(prp_ngml) sd_prp
                 from     arnold
                 group by indiv_id
                 order by indiv_id
                 ;")
prp_mean_cv = mean(prp_cv_step1$sd_prp/prp_cv_step1$mean_prp)
prp_mean_cv
arnold_cv = prp_mean_cv # for figure S7
hb_cv_step1 = sqldf("
                 select   indiv_id, avg(hb_ngml) mean_hb, stdev(hb_ngml) sd_hb
                     from     arnold
                     group by indiv_id
                     order by indiv_id
                     ;")
hb_mean_cv = mean(hb_cv_step1$sd_hb/hb_cv_step1$mean_hb)
hb_mean_cv

par(mar=c(4,4,3,5), mfrow=c(1,1))
plot(NA, NA, xlim=c(0,18), ylim=c(0,300), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=(1:9)*2-1, labels=1:9, lwd.ticks=0, lwd=0, cex.axis=1, las=1)
mtext(side=3, line=0, font=1, cex=.8, paste('PrP mean CV = ',percent(prp_mean_cv),'  Hb mean CV = ',percent(hb_mean_cv),sep=''))
abline(h=0,lwd=2)
abline(v=0,lwd=2)
abline(v=(0:9)*2,lwd=.25)
abline(h=25, lty=2, lwd=2, col='#0001CD')
mtext(side=2, line=.5, at=25, col='#0001CD', text='LLQ', las=2)
points(as.integer(arnold$id)-.5, arnold$prp_ngml, pch=19, col='#0001CD')
axis(side=1, at=arnold$x-.5, labels=arnold$visit, lwd=0, lwd.ticks=0, line=-1, cex.axis=.6)
mtext(side=2, line=0, at=-10, cex=.6, text='visit', las=2)
axis(side=2, at=100*(0:5), lwd.ticks=1, lwd=0, las=2, col.axis='#0001CD', col.ticks='#0001CD')
mtext(side=2, line=3, text='CSF [PrP] (ng/mL)', col='#0001CD')
mtext(side=1, line=2.5, text='Individual ID')
par(new=TRUE)
plot(NA, NA, xlim=c(0,18), ylim=c(1,50000), log='y', ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
points(as.integer(arnold$id)-.5, arnold$hb_ngml, pch=19, col='#710001')
axis(side=4, at=10^(0:5), labels=formatC(10^(0:5),format='fg',digits=0,big.mark=','), lwd=1, lwd.ticks=1, col.axis = '#710001', col.ticks='#710001', las=2)
abline(h=c(hb_llq, hb_ulq), lty=2, lwd=2, col='#710001')
mtext(side=4, line=.5, at=c(hb_llq, hb_ulq), col='#710001', text=c('LLQ','ULQ'), las=2)
mtext(side=4, line=3, text='CSF [Hb] (ng/mL)', col='#710001')
dev.off()




# ----------------------------------------------------------------------------------------------------------------
# Figure S7. Test-retest reliability of CSF PrP in additional cohorts

pdf('figures/manuscript/figure_s7.pdf',width=8,height=3)
par(mfrow=c(1,5), oma=c(1,1,1,1))

# test retest arnold
arnold_rel = sqldf("
select   a0.indiv_id indiv, aother.visit visit, aother.prp_ngml / a0.prp_ngml rel_prp
from     arnold a0, arnold aother
where    a0.visit = 1
and      a0.indiv_id = aother.indiv_id
order by 1, 2
;")
par(mar=c(4,5,3,1))
plot(NA, NA, xlim=c(.5,2.5), ylim=c(-4,4), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=c(1,2), labels=c(1,2), lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='visit number')
abline(h=0)
axis(side=2, at=-4:4, labels=percent(2^(-4:4)), las=2, lwd=0, lwd.ticks=1)
mtext(side=2, line=4, text='CSF [PrP] relative to first visit')
mean_cv = arnold_cv

for (indiv in unique(arnold_rel$indiv)) {
  subs = arnold_rel$indiv == indiv
  points(arnold_rel$visit[subs], log2(arnold_rel$rel_prp[subs]), type='l', lwd=3, col='red')
}
mtext(side=3, line=0, text='metformin\ntrial', font=2, cex=0.7)
mtext(side=3, line=-1, text=paste('Mean CV = ',percent(mean_cv)), cex=0.7)
mtext('A', side=3, cex=2, adj = -0.2, line = 0.3) # -- end panel A

# smaller margins for rest
par(mar=c(4,1,3,1))

# test retest kuvan
kuvan = subset(test_retest, cohort=='Swoboda')
kuvan$id = kuvan$sample_id
kuvan$x = 1:nrow(kuvan)

breaks = which(!duplicated(kuvan$indiv_id)) - 1
breaks = c(breaks, nrow(kuvan))
mids = sqldf("
             select   indiv_id, avg(x) midx
             from     kuvan
             group by 1
             order by 1
             ;")
kuvan$prp_ngml = csf$prp_ngml[match(kuvan$id, csf$sample)]
kuvan$prp_se = csf$prp_se[match(kuvan$id, csf$sample)]
kuvan$prp_l95 = csf$prp_l95[match(kuvan$id, csf$sample)]
kuvan$prp_u95 = csf$prp_u95[match(kuvan$id, csf$sample)]

cv_step1 = sqldf("
                 select   indiv_id, avg(prp_ngml) mean_prp, stdev(prp_ngml) sd_prp
                 from     kuvan
                 group by indiv_id
                 order by indiv_id
                 ;")
mean_cv = mean(cv_step1$sd_prp/cv_step1$mean_prp)
mean_cv

kuvan_rel = sqldf("
select   k0.indiv_id, kother.visit, kother.prp_ngml / k0.prp_ngml rel_prp
from     kuvan k0, kuvan kother
where    k0.visit = 1
and      k0.indiv_id = kother.indiv_id
order by 1, 2
;")

plot(NA, NA, xlim=c(.5,3.5), ylim=c(-4,4), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=1:3, labels=1:3, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='visit number')
abline(h=0)

for (indiv in unique(kuvan_rel$indiv)) {
  subs = kuvan_rel$indiv == indiv
  points(kuvan_rel$visit[subs], log2(kuvan_rel$rel_prp[subs]), type='l', lwd=3, col='red')
}
mtext(side=3, line=0, text='sapropterin\ntrial', font=2, cex=0.7)
mtext(side=3, line=-1, text=paste('Mean CV = ',percent(mean_cv)), cex=0.7)


mtext('B', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel B




m2 = subset(csf, sample %in% samples$id[samples$sample_source=='MIND'] & !grepl('\\.',sample))
tr002_match = subset(test_retest, cohort=='MIND')
tr002_match$id = gsub("MTR","",tr002_match$sample_id)
m2$iid = tr002_match$indiv_id[match(m2$sample, tr002_match$id)]
m2$visit = tr002_match$visit[match(m2$sample, tr002_match$id)]
m2 = m2[order(m2$iid, m2$visit),]
m2$x = 1:nrow(m2)
breaks = c(which(!duplicated(m2$iid)) - 1, max(m2$x))
mids = sqldf("
             select   iid, avg(x) midx
             from     m2
             group by 1
             order by 1
             ;")
tr002_match$prp_ngml = csf$prp_ngml[match(tr002_match$id, csf$sample)]
tr002_match$prp_se = csf$prp_se[match(tr002_match$id, csf$sample)]
tr002_match$prp_l95 = csf$prp_l95[match(tr002_match$id, csf$sample)]
tr002_match$prp_u95 = csf$prp_u95[match(tr002_match$id, csf$sample)]
# assign a rank of visit number within each individual -- using idea from https://stackoverflow.com/a/31859590
tr002_match = transform(tr002_match,visit=ave(1:nrow(tr002_match),indiv_id,
                                      FUN=function(x) order(visit[x],decreasing=FALSE)))


cv_step1 = sqldf("
                 select   indiv_id, avg(prp_ngml) mean_prp, stdev(prp_ngml) sd_prp
                 from     tr002_match
                 group by indiv_id
                 order by indiv_id
                 ;")
mean_cv = mean(cv_step1$sd_prp/cv_step1$mean_prp)
mean_cv


tr002_rel = sqldf("
                  select   m0.indiv_id indiv, mother.visit, mother.prp_ngml / m0.prp_ngml rel_prp
                  from     tr002_match m0, tr002_match mother
                  where    m0.visit = 1
                  and      m0.indiv_id = mother.indiv_id
                  order by 1, 2
                  ;")


plot(NA, NA, xlim=c(.5,3.5), ylim=c(-4,4), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=1:3, labels=1:3, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='visit number')
abline(h=0)
#axis(side=2, at=-4:4, labels=percent(2^(-4:4)), las=2, lwd=0, lwd.ticks=1)
#mtext(side=2, line=4, text='CSF [PrP] relative to first visit')

for (indiv in unique(tr002_rel$indiv)) {
  subs = tr002_rel$indiv == indiv
  points(tr002_rel$visit[subs], log2(tr002_rel$rel_prp[subs]), type='l', lwd=3, col='red')
}
mtext(side=3, line=0, text='MIND NPH\nlumbar drains', font=2, cex=0.7)
mtext(side=3, line=-1, text=paste('Mean CV = ',percent(mean_cv)), cex=0.7)



mtext('C', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel C


# -- begin panels D & E


ucsf_tr = subset(elisa, plate %in% c(33,34,38,39) & !grepl('IPC',sample) & !grepl('pool',sample))
ucsf_tr$l95 = ucsf_tr$prp_ngml - 1.96*ucsf_tr$prp_se
ucsf_tr$u95 = ucsf_tr$prp_ngml + 1.96*ucsf_tr$prp_se
ucsf_tr$indiv = as.integer(gsub('[A-Z]+','',ucsf_tr$sample))


ucsf_match = subset(test_retest, cohort=='Geschwind')
ucsf_tr$visit = ucsf_match$visit[match(ucsf_tr$sample,ucsf_match$sample_id)]
ucsf_tr$indiv_id = ucsf_match$indiv_id[match(ucsf_tr$sample,ucsf_match$sample_id)]

# get symptomatic/presymptomatic status
ucsf_tr$status = samples$prion_category[match(ucsf_tr$sample, samples$id)]
ucsf_n = sqldf("select indiv_id, count(*) n from ucsf_tr group by 1 order by 1;")
ucsf_tr$n_samples = ucsf_n$n[match(ucsf_tr$indiv_id, ucsf_n$indiv_id)]

# for test-retest, only care about those with >1 sample
ucsf_tr = subset(ucsf_tr, n_samples > 1)

ucsf_tr = ucsf_tr[with(ucsf_tr, order(indiv_id, visit)),]

ucsf_presymp = subset(ucsf_tr, status=='presymptomatic genetic')
ucsf_symp = subset(ucsf_tr, status=='symptomatic genetic')

ucsf_presymp$x = 1:nrow(ucsf_presymp)

breaks = c(ucsf_presymp$x[which(!duplicated(ucsf_presymp$indiv))] - 0.5, max(ucsf_presymp$x) + 0.5)
mids = sqldf("
             select   indiv, avg(x) midx
             from     ucsf_presymp
             group by 1
             order by 1
             ;")

ucsf_presymp_step1 = sqldf("
                           select   indiv, avg(prp_ngml) mean_prp, stdev(prp_ngml) sd_prp
                           from     ucsf_presymp
                           group by indiv
                           having   count(*) > 1
                           order by indiv
                           ;")
mean_cv = mean(ucsf_presymp_step1$sd_prp/ucsf_presymp_step1$mean_prp)
mean_cv

ucsf_presymp$indiv_id = as.character(ucsf_presymp$indiv_id)

ucsf_presymp_rel = sqldf("
                  select   p0.indiv_id, pother.visit, pother.prp_ngml / p0.prp_ngml rel_prp
                  from     ucsf_presymp p0, ucsf_presymp pother
                  where    p0.visit = 1
                  and      p0.indiv_id = pother.indiv_id
                  order by 1, 2
                  ;")
colnames(ucsf_presymp_rel)[1] = 'indiv'


plot(NA, NA, xlim=c(.5,5.5), ylim=c(-4,4), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=1:5, labels=1:5, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='visit number')
abline(h=0)
#axis(side=2, at=-4:4, labels=percent(2^(-4:4)), las=2, lwd=0, lwd.ticks=1)
#mtext(side=2, line=4, text='CSF [PrP] relative to first visit')

for (indiv in unique(ucsf_presymp_rel$indiv)) {
  subs = ucsf_presymp_rel$indiv == indiv
  points(ucsf_presymp_rel$visit[subs], log2(ucsf_presymp_rel$rel_prp[subs]), type='l', lwd=3, col='red')
}
mtext(side=3, line=0, text='UCSF\npresymptomatic', font=2, cex=0.7)
mtext(side=3, line=-1, text=paste('Mean CV = ',percent(mean_cv)), cex=0.7)


mtext('D', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel D




# ----------begin panel E

ucsf_symp$x = 1:nrow(ucsf_symp)

breaks = c(ucsf_symp$x[which(!duplicated(ucsf_symp$indiv))] - 0.5, max(ucsf_symp$x) + 0.5)
mids = sqldf("
             select   indiv, avg(x) midx
             from     ucsf_symp
             group by 1
             order by 1
             ;")

ucsf_symp_step1 = sqldf("
                        select   indiv, avg(prp_ngml) mean_prp, stdev(prp_ngml) sd_prp
                        from     ucsf_symp
                        group by indiv
                        having   count(*) > 1
                        order by indiv
                        ;")
mean_cv = mean(ucsf_symp_step1$sd_prp/ucsf_symp_step1$mean_prp)
mean_cv



ucsf_symp$indiv_id = as.character(ucsf_symp$indiv_id)

ucsf_symp_rel = sqldf("
                      select   p0.indiv_id, pother.visit, pother.prp_ngml / p0.prp_ngml rel_prp
                      from     ucsf_symp p0, ucsf_symp pother
                      where    p0.visit = 1
                      and      p0.indiv_id = pother.indiv_id
                      order by 1, 2
                      ;")
colnames(ucsf_symp_rel)[1] = 'indiv'


plot(NA, NA, xlim=c(.5,5.5), ylim=c(-4,4), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=1:5, labels=1:5, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='visit number')
abline(h=0)

for (indiv in unique(ucsf_symp_rel$indiv)) {
  subs = ucsf_symp_rel$indiv == indiv
  points(ucsf_symp_rel$visit[subs], log2(ucsf_symp_rel$rel_prp[subs]), type='l', lwd=3, col='red')
}
mtext(side=3, line=0, text='UCSF\nsymptomatic', font=2, cex=0.7)
mtext(side=3, line=-1, text=paste('Mean CV = ',percent(mean_cv)), cex=0.7)


mtext('E', side=3, cex=2, adj = 0.0, line = 0.3) # -- end panel E


dev.off() # -- end Figure S7







