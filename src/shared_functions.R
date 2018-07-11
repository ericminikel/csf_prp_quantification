options(stringsAsFactors=FALSE)
library(reshape2)
library(sqldf)
library(plyr)
setwd('~/d/sci/src/csf_prp_quantification/')

expand_range = function(x, by=.5) {
  return ( c(min(x)-by,max(x)+by) )
}

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}


# take a path to a raw datafile (.txt) from the spectramax platereader, and return a data frame
process_spectramax_data = function(path) {
  file_lines = readLines(path)
  lines_to_skip = suppressWarnings(grep('Plate:',file_lines) + 1)
  rawdata = read.table(path,sep='\t',header=F,skip=lines_to_skip,fill=T)
  if (ncol(rawdata) == 16) {
    rawdata = as.matrix(rawdata[1:8,3:14])
    colnames(rawdata) = 1:12
    rownames(rawdata) = LETTERS[1:8]
    rawdata = data.frame(melt(rawdata))
    wavelength = strsplit(suppressWarnings(readLines(path))[lines_to_skip-1],"\t")[[1]][16]
    a_wavelength = paste('a',wavelength,sep='')
    colnames(rawdata) = c('row','col',a_wavelength)
    rawdata$row = as.character(rawdata$row)
  } else if (ncol(rawdata) == 29) {
    rawdata = rawdata[1:8,c(3:14,16:27)]
    abs450 = as.matrix(rawdata[1:8,1:12])
    abs620 = as.matrix(rawdata[1:8,13:24])
    abs450_620 = abs450 - abs620
    colnames(abs450_620) = 1:12
    rownames(abs450_620) = LETTERS[1:8]
    rawdata = data.frame(melt(abs450_620))
    colnames(rawdata) = c('row','col','a450_620')
    rawdata$row = as.character(rawdata$row)
  }
  return (rawdata)
}

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

castable_to_integer = function(x) {
  tryCatch ({
    as.integer(x)
    return (TRUE)
  }, warning = function(w) {
    return (FALSE)
  })
}

# 
# # test
# filename = 'totprot4-2016-11-17.CSV'
# filename = 'SV-2016-09-30-123319-CSF-tp.txt'
# full_path = paste('data/dc/raw/',filename,sep='')
# path = full_path
# 
# filename = 'csf-elisa-08-30-16-164815.txt'
# full_path = paste('data/elisa/raw/',filename,sep='')
# path = full_path
# 
