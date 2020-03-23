#Call ancestry in all subjects, using R
#V4 - Mar 01, 2017 - New African-American cutoff
#V3 - Jun 17, 2016 - Export predpc file with header
#V2 - May 31, 2016 - Made plots and outputs pretty
#V1 - May 5 2016 - Initial release 

   	#WR colorlist
	colorlist <- rep(rgb(0,0,0,1), 8)
	colorlist[1] <- rgb(242,24,3,255, maxColorValue = 255) #1 are africans
	colorlist[2] <- rgb(0,6,253,255, maxColorValue = 255) #2 europeans
	colorlist[3] <- rgb(142,0,117,255, maxColorValue = 255) #3 han chinese
	colorlist[4] <- rgb(178,178,55,255, maxColorValue = 255) #4 sw asians
	colorlist[5] <- rgb(68,26,80,255, maxColorValue = 255)   #5 native americans
	colorlist[6] <- rgb(0,153,229,255, maxColorValue = 255) #6 oceanians
	colorlist[7] <- rgb(255,0,249,255, maxColorValue = 255) #7 cs asian
	colorlist[8] <- rgb(122,124,166,255, maxColorValue = 255) #8 Admixed

  	 
  popassign2 <- function(x,panel_calibration)
	{	
		#Column names of populations (In order of SNPweights reference file)
		europe <- 1
		oceania <- 2
		africa <- 3
		americas <- 4
		cs_asia <- 5
		e_asia <- 6
		
		eur <-  (x[europe] >= 0.90 )
		csa <-  (x[cs_asia] >= 0.90 )
		#afr <-  (x[africa] >= 0.90)
		eas <-  (x[e_asia] >= 0.90 )	
		aam <-  (x[europe] < 0.90 & x[oceania] < 0.05 & x[africa] >= 0.05 &  x[americas] < 0.05  & x[cs_asia] < 0.05 & x[e_asia] < 0.05 ) #AAM must be <90% African to insure exclusivity # 
		lat <-  (x[europe] < 0.90 & x[oceania] < 0.05 & x[africa] < 0.07 & x[americas] >= 0.05 & x[cs_asia] < 0.05 & x[e_asia] < 0.05)
		nat <-  (x[americas] >= 0.60 & x[oceania] < 0.05 & x[africa] < 0.05 & x[cs_asia] < 0.15 & x[e_asia] < 0.20 )
		pue <-  (x[europe] >= 0.40 & x[oceania] < 0.05 & x[africa] >= 0.07 & x[americas] >= 0.05 & x[cs_asia] < 0.05 & x[e_asia] < 0.05)
		oce <-  (x[oceania] >= 0.2)
		fil <-  (x[europe] >= 0.4 & x[oceania] < 0.05 & x[africa] < 0.05 & x[americas] < 0.05  & x[cs_asia] < 0.05 & x[e_asia] >= 0.4 )
		afr2 <- (x[europe] < 0.5 & x[oceania] < 0.05 & x[africa] >= 0.5 &  x[americas] < 0.05  & (x[cs_asia] + x[e_asia]) < 0.10 ) #New African-American addition , more lenient
                
	        #extra adjustment for Affy UKBB data, eur tends to clsuter with CSA, this is calibrated to matching self-report
	        if(panel_calibration == "AffyBB")
		 {
	           eur <- (x[europe] >= 0.75 & x[oceania] < 0.05 & x[africa] < 0.05 &  x[americas] < 0.05  & x[cs_asia] < 0.2 & x[e_asia] < 0.05 )
                 }
	  
		if(eur)
		{
			popfit <- "eur"
		}
		if(csa)
		{
			popfit <- "csa"
		}
		if(aam) #defines african americans
		{
			popfit <- "aam"
		}  
		if(!aam & afr2) #defines african americans
		{
			popfit <- "aam"
		}  
            
		#if(afr) #defines africans
		#{
		#	popfit <- "afr"
		#}
		if(eas) #defines east asia
		{
			popfit <- "eas"
		}
		if(lat) #latinos
		{
			popfit <- "lat"
		}
		if(!lat & nat) #latinos
		{
			popfit <- "nat"
		}
		if(oce) #oceanians
		{
			popfit <- "oce"
		}
		if(fil) #filipinos
		{
			popfit <- "fil"
		}	
		if(pue) #puerto rico
		{
			popfit <- "pue"
		}	
	
		#Others are defined as anyone not fitting into a specified group
		other <- eur + max(afr2,aam) + eas +    csa + max(nat,lat) + oce +fil + pue # max(afr,aam)
			
		if(other == 0)
		{
			popfit <- "oth"
		}
		#error handling for if a subject fits into more than 1 group
		if (other > 1)
		{
			popfit <- "ERROR"
		}
		return (popfit)
	}
args <- commandArgs(trailingOnly = TRUE)

bfilepreds <- args[1]
clustercenterpreds <- args[2]
panel_calibration <- args[3]


    unlist_split <- function(x, ...)
	{
		toret <- unlist(strsplit(x, ...) )
		return(t(toret))
	}

#Read 1 week paper SNPweights predictions

	datam <- read.table(bfilepreds, header=F,stringsAsFactors=F,na.strings=c("#N/A"))
	names(datam) <- c("FID_IID","affection", "markers","PC1","PC2","PC3","PC4","PC5","europe", "oceania", "africa", "americas", "cs_asia", "e_asia")
	
       datam$FID <- NA
	datam$IID <- NA
	datam[,c("FID","IID")] <- t(sapply(datam$FID_IID,unlist_split,split=":"))
    print("assigning populations to test data")
       datam$bestpop_oneweek <- apply(datam[,9:14], 1, popassign2,panel_calibration=panel_calibration)
       write.table(datam, paste(bfilepreds,'.header',sep=''),row.names=F,quote=F) #Write predPC out to a file with a header
       
print("Assigning colors to test data")
 datam$color <- colorlist[8]

 datam$color <- ifelse(datam$bestpop_oneweek == "eur", colorlist[2], datam$color)
 datam$color <- ifelse(datam$bestpop_oneweek == "afr", colorlist[1], datam$color)
 datam$color <- ifelse(datam$bestpop_oneweek == "aam", colorlist[1], datam$color)
 datam$color <- ifelse(datam$bestpop_oneweek == "lat", colorlist[5], datam$color)
 datam$color <- ifelse(datam$bestpop_oneweek == "eas", colorlist[3], datam$color)
 datam$color <- ifelse(datam$bestpop_oneweek == "csa", colorlist[7], datam$color)
 datam$color <- ifelse(datam$bestpop_oneweek == "nat", "green", datam$color)
 datam$color <- ifelse(datam$bestpop_oneweek == "oce", colorlist[6], datam$color)
 datam$color <- ifelse(datam$bestpop_oneweek == "fil", "black", datam$color)
 datam$color <- ifelse(datam$bestpop_oneweek == "pue", "#CC6699", datam$color)

 pop_list <- c("aam", "afr", "csa", "eas", "eur","fil", "lat", "nat", "pue", "oce", "oth")
 pop_listlong <- c("AfrAm/African","Alaska/NatAm","CS Asian", "East Asian (EA)", "European (EU)", "\"Filipino\"/EU-EA", "Latino/NatAm",  "Pacific Isl./Oceania", "\"Puerto Rican\"/3 way", "Other")
 popcolors <- c(colorlist[1],"green",colorlist[7],colorlist[3],colorlist[2],"black",colorlist[5],colorlist[6],"#CC6699",colorlist[8])


#Read reference data
 refdat <- read.csv(clustercenterpreds,stringsAsFactors=F,header=T)
 print(refdat)
	refdat$color <- NA

	refdat[which(refdat[,1] == "africa"),]$color <- colorlist[1]
	refdat[which(refdat[,1] == "europe"),]$color <- colorlist[2]
	refdat[which(refdat[,1] == "e_asia"),]$color <- colorlist[3]
	refdat[which(refdat[,1] == "cs_asia"),]$color <- colorlist[7]
	refdat[which(refdat[,1] == "americas"),]$color <- colorlist[5]
	refdat[which(refdat[,1] == "oceania"),]$color <- colorlist[6]

 pdf(paste(bfilepreds,'_pcs.pdf',sep=''),7,7)
	#Plot PCs
	plot(-datam$PC1,datam$PC2, col=datam$color,  pch="x", cex=.75, xlim=c(-0.025,0.04), ylim=c(-0.03,0.025),xlab="PC1",ylab="PC2",cex.axis=1.25,cex.lab=1.45)
	 points(-refdat$PC1,refdat$PC2,col='black', bg=refdat$color, pch=21,cex=1.5)
	 points(-refdat$PC1,refdat$PC2,col='black', bg=refdat$color, pch="R",cex=.75)

       legend('bottomright',col=popcolors, legend=pop_listlong, pch=19,cex=.6)
	plot(datam$PC2,datam$PC3, col=datam$color, pch="x", cex=.75 , xlim=c(-0.035,0.025), ylim=c(-0.105,0.11),xlab="PC2",ylab="PC3",cex.axis=1.25,cex.lab=1.45)
	 points(refdat$PC2,refdat$PC3,col='black', bg=refdat$color, pch=21,cex=1.5)
	 points(refdat$PC2,refdat$PC3,col='black', bg=refdat$color, pch="R",cex=.75)
	legend('bottomright',col=popcolors, legend=pop_listlong, pch=19,cex=.6)
	plot(-datam$PC1,datam$PC4, col=datam$color, pch="x", cex=.75 , xlim=c(-0.025,0.04), ylim=c(-0.075,0.02),xlab="PC1",ylab="PC4",cex.axis=1.25,cex.lab=1.45)
	 points(-refdat$PC1,refdat$PC4,col='black', bg=refdat$color, pch=21,cex=1.5)
	 points(-refdat$PC1,refdat$PC4,col='black', bg=refdat$color, pch="R",cex=.75)
	legend('bottomright',col=popcolors, legend=pop_listlong, pch=19,cex=.6)
	plot(-datam$PC1,datam$PC5, col=datam$color, pch="x", cex=.75 , xlim=c(-0.025,0.04), ylim=c(-0.15,0.025),xlab="PC1",ylab="PC5",cex.axis=1.25,cex.lab=1.45)
	 points(-refdat$PC1,refdat$PC5,col='black', bg=refdat$color, pch=21,cex=1.5)
	 points(-refdat$PC1,refdat$PC5,col='black', bg=refdat$color, pch="R",cex=.75)
	legend('bottomright',col=popcolors, legend=pop_listlong, pch=19,cex=.6)

 dev.off()


 
    write.table(subset(datam,select=c(FID,IID,bestpop_oneweek)),paste(bfilepreds,'_ancestries.txt',sep=''),quote=F,row.names=F)
    write.table(table(datam$bestpop_oneweek), paste(bfilepreds,'_ancestries_samplesizes.txt',sep=''),quote=F,row.names=F)
