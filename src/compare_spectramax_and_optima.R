options(stringsAsFactors=FALSE)
library(reshape2)
library(sqldf)
library(plyr)
setwd('~/d/sci/src/csf_prp_quantification/')

spectramax = read.table('data/elisa/processed/06.tsv',sep='\t',header=TRUE)
optima = read.table('data/elisa/processed/6_alt.tsv',sep='\t',header=TRUE)

both = sqldf("
select   s.row, s.col, s.abs450_620 sa, o.abs450_620 oa, s.ngml sp, o.ngml op
from     spectramax s, optima o
where    s.row = o.row
and      s.col = o.col
order by 1,2
;")

png('figures/spectramax_optima_compare_absorbance.png',width=800,height=500,res=150)
par(mar=c(4,4,3,1))
plot(NA,NA,xlim=c(0,3.5),ylim=c(0,3.5),ann=FALSE,axes=FALSE,xaxs='i',yaxs='i')
axis(side=1, at=(0:7)/2, lwd=0, lwd.ticks=1)
axis(side=2, at=(0:7)/2, lwd=0, lwd.ticks=1, las=2)
abline(v=0, lwd=2)
abline(h=0, lwd=2)
abline(a=0, b=1, col='red', lty=2, lwd=2)
points(x=both$sa, y=both$oa, pch=20, col='#555578')
mtext(side=1, line=2.5, text='Spectramax A450-A620')
mtext(side=2, line=3, text='Fluostar Optima A450-A620')
mtext(side=3, line=1, text='Absorbance comparison between 2 platereaders', font=2, cex=1.2)
mtext(side=3, line=0, text='CSF ELISA 6', cex=.9)
dev.off()

png('figures/spectramax_optima_compare_protein_conc.png',width=800,height=500,res=150)
par(mar=c(4,4,3,1))
plot(NA,NA,xlim=c(0,300),ylim=c(0,300),ann=FALSE,axes=FALSE,xaxs='i',yaxs='i')
axis(side=1, at=(0:3)*100, lwd=0, lwd.ticks=1)
axis(side=2, at=(0:3)*100, lwd=0, lwd.ticks=1, las=2)
abline(v=0, lwd=2)
abline(h=0, lwd=2)
abline(a=0, b=1, col='red', lty=2, lwd=2)
points(x=both$sp, y=both$op, pch=20, col='#555578')
mtext(side=1, line=2.5, text='Spectramax calculated [PrP] (ng/mL)')
mtext(side=2, line=3, text='Fluostar Optima calculated [PrP] (ng/mL)')
mtext(side=3, line=1, text='Protein results between 2 platereaders', font=2, cex=1.2)
mtext(side=3, line=0, text='CSF ELISA 6', cex=.9)
dev.off()

meta = read.table('data/elisa/meta/meta_plate06.tsv',sep='\t',header=TRUE)

# calculate mean cv for all wells
replicates = sqldf("
select   m.detail, avg(b.sp) sp_mean, stdev(b.sp) sp_sd, avg(b.op) op_mean, stdev(b.op) op_sd
from     both b, meta m
where    b.row = m.row
and      b.col = m.col
group by 1
order by 1
;")
s_mean_cv = mean(replicates$sp_sd/replicates$sp_mean)
o_mean_cv = mean(replicates$op_sd/replicates$op_mean)
s_mean_cv
o_mean_cv

# calculate mean cv just for CSF
replicates = sqldf("
select   m.detail, avg(b.sp) sp_mean, stdev(b.sp) sp_sd, avg(b.op) op_mean, stdev(b.op) op_sd
from     both b, meta m
where    b.row = m.row
and      b.col = m.col
and      m.stype = 'CSF'
group by 1
order by 1
;")
s_mean_cv = mean(replicates$sp_sd/replicates$sp_mean)
o_mean_cv = mean(replicates$op_sd/replicates$op_mean)
s_mean_cv
o_mean_cv

replicates$s_cv = replicates$sp_sd/replicates$sp_mean
replicates$o_cv = replicates$op_sd/replicates$op_mean

# exploratory plots
# plot(replicates$sp_mean, replicates$s_cv)
# plot(replicates$op_mean, replicates$o_cv)
# nothing interesting there

png('figures/spectramax_optima_compare_cv.png',width=800,height=500,res=150)
par(mar=c(4,4,3,1))
plot(NA,NA,xlim=c(0,.3),ylim=c(0,.3),ann=FALSE,axes=FALSE,xaxs='i',yaxs='i')
axis(side=1, at=(0:3)/10, lwd=0, lwd.ticks=1)
axis(side=2, at=(0:3)/10, lwd=0, lwd.ticks=1, las=2)
abline(v=0, lwd=2)
abline(h=0, lwd=2)
abline(a=0, b=1, col='red', lty=2, lwd=2)
points(x=replicates$s_cv, y=replicates$o_cv, pch=20, col='#555578')
mtext(side=1, line=2.5, text='Spectramax CV')
mtext(side=2, line=3, text='Fluostar Optima CV')
mtext(side=3, line=1, text='CV between 2 platereaders', font=2, cex=1.2)
mtext(side=3, line=0, text='CSF ELISA 6', cex=.9)
dev.off()
