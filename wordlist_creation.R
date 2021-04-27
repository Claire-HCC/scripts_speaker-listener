rm(list=ls(all=TRUE))
library(tidytext)
library(dplyr)
library(plyr)
library(stringr)

expdir='Y://claire//speaker-listener//'
exp='sherlock'
text=read.delim(paste(expdir,'//',exp,'//sound//text.txt',sep=""),stringsAsFactors=FALSE)

words=unnest_tokens(text,word,event)


write.table(words,col.names=TRUE,row.names = FALSE,file=paste(expdir,'//',exp,'//sound//words.txt',sep=""))

# https://clarin.phonetik.uni-muenchen.de/BASWebServices/help