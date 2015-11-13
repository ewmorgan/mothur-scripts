#!/usr/bin/Rscript
ar <- commandArgs(TRUE)

help_msg <- c('Script to calculate number of OTUs and sequences','shared between two,three, four of six sample groups.','\tUsage: count_shared_seqs.R SHARED_FILE [NUM_GROUPS ["GROUP NAMES"]]')


shared<-read.table(ar[1], header=TRUE, fill=TRUE)
if(length(ar) == 1) {
	if(length(grep('--help', ar[1])) == 1) {
		writeLines(paste(help_msg,collapse="\n"))
	}

	ar[2] <- nrow(shared)
	ar[3] <- paste0(shared[,2],collapse=' ')
} else if(length(ar)<3) {
	writeLines(paste(help_msg,collapse="\n"))
	writeLines('Provide NUM_GROUPS and "GROUP NAMES" as second and third arguments.')
	quit()
}

if(!(ar[2] %in% c(2,3,4,6))) {
	writeLines(paste(help_msg,collapse="\n"))
	writeLines("Can not process this number of groups!")
	print(ar)
	quit()
}
ar[3] <- paste0(sort(strsplit(ar[3],' ')[[1]]),collapse=' ')
#print(ar)

getIntersectsOf6<-function(g, s)
{
#looking for unique OTUs only

	sh<-s[s$Group %in% g, 4:length(s[1, ])]
	o<-list()
#unique
	i<-paste(c('+','-','-','-','-','-'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]==0 & sh[3,]==0 & sh[4,]==0 & sh[5,]==0 & sh[6,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{	#this check is required if there are initially more than 6 groups 
		o[[i]]<-c(sum(ish), sum(sh[,ish]), mapply(c,colnames(sh[1,ish]),sh[1,ish])
		)
	}

	i<-paste(c('-','+','-','-','-','-'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]>0 & sh[3,]==0 & sh[4,]==0 & sh[5,]==0 & sh[6,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{	#this check is required if there are initially more than 6 groups 
		o[[i]]<-c(sum(ish), sum(sh[,ish]), mapply(c,colnames(sh[2,ish]),sh[2,ish])
		
		)
	}

	i<-paste(c('-','-','+','-','-','-'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]==0 & sh[3,]>0 & sh[4,]==0 & sh[5,]==0 & sh[6,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{	#this check is required if there are initially more than 6 groups 
		o[[i]]<-c(sum(ish), sum(sh[,ish]), mapply(c,colnames(sh[3,ish]),sh[3,ish])
		)
	}

	i<-paste(c('-','-','-','+','-','-'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]==0 & sh[3,]==0 & sh[4,]>0 & sh[5,]==0 & sh[6,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{	#this check is required if there are initially more than 6 groups 
		o[[i]]<-c(sum(ish), sum(sh[,ish]), mapply(c,colnames(sh[4,ish]),sh[4,ish]))
	}


	i<-paste(c('-','-','-','-','+','-'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]==0 & sh[3,]==0 & sh[4,]==0 & sh[5,]>0 & sh[6,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{	#this check is required if there are initially more than 6 groups 
		o[[i]]<-c(sum(ish), sum(sh[,ish]), mapply(c,colnames(sh[5,ish]),sh[5,ish])
		)
	}

	i<-paste(c('-','-','-','-','-','+'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]==0 & sh[3,]==0 & sh[4,]==0 & sh[5,]==0 & sh[6,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{	#this check is required if there are initially more than 6 groups 
		o[[i]]<-c(sum(ish), sum(sh[,ish]), mapply(c,colnames(sh[6,ish]),sh[6,ish])
		)
	}

	o
}

getIntersectsOf4<-function(g, s)
{
	sh<-s[s$Group %in% g, 4:length(s[1, ])]
	o<-list()

	#print(s[1:4,1:5])
	#print(sh[1:4,1:5])
	#quit()

#shared between four
	i<-paste(c('+','+','+','+'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]>0 & sh[3,]>0 & sh[4,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

#shared between three
	i<-paste(c('+','+','+','-'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]>0 & sh[3,]>0 & sh[4,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('+','+','-','+'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]>0 & sh[3,]==0 & sh[4,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('+','-','+','+'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]==0 & sh[3,]>0 & sh[4,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('-','+','+','+'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]>0 & sh[3,]>0 & sh[4,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

#shared between two
	i<-paste(c('+','+','-','-'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]>0 & sh[3,]==0 & sh[4,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('+','-','+','-'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]==0 & sh[3,]>0 & sh[4,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('+','-','-','+'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]==0 & sh[3,]==0 & sh[4,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('-','+','+','-'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]>0 & sh[3,]>0 & sh[4,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('-','+','-','+'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]>0 & sh[3,]==0 & sh[4,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('-','-','+','+'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]==0 & sh[3,]>0 & sh[4,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

#unique
	i<-paste(c('+','-','-','-'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]==0 & sh[3,]==0 & sh[4,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('-','+','-','-'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]>0 & sh[3,]==0 & sh[4,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('-','-','+','-'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]==0 & sh[3,]>0 & sh[4,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('-','-','-','+'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]==0 & sh[3,]==0 & sh[4,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	o
}

getIntersectsOf3<-function(g, s)
{
	sh<-s[s$Group %in% g, 4:length(s[1, ])]
	o<-list()

	#print(s[1:4,1:5])
	#print(sh[1:4,1:5])
	#quit()
	
	i<-paste(rep('+', 3),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]>0 & sh[3,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('+','+','-'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]>0 & sh[3,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('+','-','+'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]==0 & sh[3,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('-','+','+'),g,collapse='', sep='')
	ish<-sh[1,]==0 & sh[2,]>0 & sh[3,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('+','-','-'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]==0 & sh[3,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}
	i<-paste(c('-','+','-'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]>0 & sh[3,]==0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	i<-paste(c('-','-','+'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]==0 & sh[3,]>0
	if(sum(ish)==0){o[[i]]<-c(0,0)}else{
		o[[i]]<-c(sum(ish), sum(sh[,ish]))
	}

	o
}

getIntersectsOf2<-function(g, s)
{
	sh<-s[s$Group %in% g, 4:length(s[1, ])]
	o<-list()

	i<-paste(rep('+', 2),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]>0
	o[[i]]<-c(sum(ish), sum(sh[,ish]))

	i<-paste(c('+','-'),g,collapse='',sep='')
	ish<-sh[1,]>0 & sh[2,]==0
	o[[i]]<-c(sum(ish), sum(sh[,ish]))

	i<-paste(c('-','+'),g,collapse='',sep='')
	ish<-sh[1,]==0 & sh[2,]>0
	o[[i]]<-c(sum(ish), sum(sh[,ish]))

	o
}

if(ar[2]==2) {
	getIntersectsOf2(g=strsplit(ar[3]," ")[[1]], s=shared)
} else if (ar[2]==3) {
	getIntersectsOf3(g=strsplit(ar[3]," ")[[1]], s=shared)
} else if (ar[2]==4) {
	getIntersectsOf4(g=strsplit(ar[3]," ")[[1]], s=shared)
} else if (ar[2]==6) {
	getIntersectsOf6(g=strsplit(ar[3]," ")[[1]], s=shared)
} else {}

quit()