#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ACE_functions.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Obtain functions used in ACE.R
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1) Add function to plot summary files (ACE adjusted segments for each fit)
#
# History:
#  23-03-2022: File creation, write code
#  30-03-2022: Create function to fix segmentdata
#  31-03-2022: Implement both squaremodel and fixed ploidy ACE
#  08-04-2022: Create Add_normal_segments, edit concatenate lengths
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
suppressWarnings(library(ggplot2))
#-------------------------------------------------------------------------------
segmentstotemplate <- function(segmentdf, chrci=2, startci=3, endci=4, binsci=5, meanci=6, seci, sdci, log=FALSE) {
  if(is(segmentdf, "character")) {segmentdf <- try(read.table(segmentdf, header = TRUE, comment.char = "", sep = "\t"))}
  if (inherits(segmentdf, "try-error")) {print(paste0("could not find ", segmentdf))
  } else {
    if(log==TRUE){log<-exp(1)}
    if(log) {
      segmentdf[,meanci] <- log^segmentdf[,meanci]
    }
    segments <- rep(segmentdf[,meanci],segmentdf[,binsci])
    chr <- rep(segmentdf[,chrci],segmentdf[,binsci])
    chr <- gsub("chr","",chr,ignore.case = TRUE)
    if(!missing(sdci)) {
      copynumbers <- c()
      for (i in seq_along(segmentdf[,1])) {copynumbers <- append(copynumbers,rnorm(n=segmentdf[i,binsci],
                                                                                  mean=segmentdf[i,meanci],
                                                                                  sd=segmentdf[i,sdci]))}
    } else if(!missing(seci)) {
      copynumbers <- c()
      for (i in seq_along(segmentdf[,1])) {copynumbers <- append(copynumbers,rnorm(n=segmentdf[i,binsci],
                                                                                  mean=segmentdf[i,meanci],
                                                                                  sd=(segmentdf[i,seci]*sqrt(segmentdf[i,binsci]))))}
    } else {
      copynumbers <- segments
    }
    bin <- seq_along(segments)
    start <- c()
    end <- c()
    for (i in seq_along(segmentdf[,1])) {
      start <- append(start,round(seq(from=segmentdf[i,startci],by=(segmentdf[i,endci]+1-segmentdf[i,startci])/segmentdf[i,binsci],
                                length.out=segmentdf[i,binsci])))
      end <- append(end,round(seq(to=segmentdf[i,endci],by=(segmentdf[i,endci]+1-segmentdf[i,startci])/segmentdf[i,binsci],
                            length.out=segmentdf[i,binsci])))
    }
    template <- data.frame(bin,chr,start,end,copynumbers,segments)
    return(template)
  }
}
squaremodel <- function(template, QDNAseqobjectsample = FALSE, prows=100, ptop=5, pbottom=1, method = 'RMSE',
                        exclude = c("X","Y"), sgc = c(), penalty = 0, penploidy = 0, 
                        cellularities = seq(5,100), highlightminima = TRUE, standard,errorgraph_out) {
  if(QDNAseqobjectsample) {template <- objectsampletotemplate(template, QDNAseqobjectsample)}
  template <- template[!template$chr %in% exclude,]
  if(missing(standard) || !is(standard, "numeric")) { standard <- median(template$segments, na.rm = T) }
  templateh <- template[template$chr %in% sgc,]
  template <- template[!template$chr %in% sgc,]
  segmentdata <- rle(as.vector(na.exclude(template$segments)))
  segmentdatah <- rle(as.vector(na.exclude(templateh$segments)))
  fraction <- sort(cellularities)/100
  error <- c()
  errormatrix <- matrix(nrow=(prows+1),ncol=length(fraction))
  listofploidy <- c()
  listofcellularity <- c()
  listoferrors <- c()
  for (t in seq(0,prows)) {
    ploidy <- ptop-((ptop-pbottom)/prows)*t
    listofploidy <- append(listofploidy, rep(ploidy,length(fraction)))
    expected <- c()
    hexpected <- c()
    temp <- c()
    errorlist <- c()
    for (i in seq_along(fraction)) {
      for (p in seq(1,12)) {
        expected[p] <- standard*(p*fraction[i] + 2*(1-fraction[i]))/(fraction[i]*ploidy + 2*(1-fraction[i]))
        hexpected[p] <- standard*(p*fraction[i] + 1*(1-fraction[i]))/(fraction[i]*ploidy + 2*(1-fraction[i]))
      }
      for (j in seq_along(segmentdata$values)) {
        if(method=='RMSE') {temp[j] <- (min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty))^2}
        else if(method=='SMRE') {temp[j] <- sqrt(min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty))}
        else if(method=='MAE') {temp[j] <- min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty)}
        else {print("Not a valid method")}
      }
      for (k in seq_along(segmentdatah$values)) {
        if(method=='RMSE') {temp[j+k] <- (min(abs(segmentdatah$values[k]-hexpected),0.5)/(fraction[i]^penalty))^2}
        else if(method=='SMRE') {temp[j+k] <- sqrt(min(abs(segmentdatah$values[k]-hexpected),0.5)/(fraction[i]^penalty))}
        else if(method=='MAE') {temp[j+k] <- min(abs(segmentdatah$values[k]-hexpected),0.5)/(fraction[i]^penalty)}
        else {print("Not a valid method")}
      }
      lengths <- c(segmentdata$lengths, segmentdatah$lengths)
      if(method=='RMSE') {errorlist[i] <- sqrt(sum(temp*lengths)/sum(lengths))}
      else if(method=='SMRE') {errorlist[i] <- sum(temp*lengths)/sum(lengths)^2}
      else if(method=='MAE') {errorlist[i] <- sum(temp*lengths)/sum(lengths)}
    }
    listofcellularity <- append(listofcellularity, fraction)
    listoferrors <- append(listoferrors, errorlist)
    errormatrix[t+1,] <- errorlist
    
  }
  minimat <- matrix(nrow=(prows+1),ncol=length(fraction))
  for (i in seq(1,prows+1)) {
    for (j in seq(1,length(fraction))) {
      if (i==1|i==(prows+1)) {
        minimat[i,j] <- FALSE
      } else if (j==length(fraction)) {
        minimat[i,j] <- errormatrix[i,j]==min(errormatrix[seq(i-1,i+1),seq(j-1,j)])
      } else {
        minimat[i,j] <- errormatrix[i,j]==min(errormatrix[seq(i-1,i+1),seq(j-1,j+1)])
      }
    }
  }
  round(errormatrix, digits = 10)
  round(listoferrors, digits = 10)
  errordf <- data.frame(ploidy=listofploidy,
                        cellularity=listofcellularity,
                        error=listoferrors/max(listoferrors),
                        minimum=as.vector(t(minimat)))
  minimadf <- errordf[errordf$minimum==TRUE,]
  minimadf <- minimadf[order(minimadf$error,-minimadf$cellularity),]
  if(highlightminima==TRUE){
      tempplot <- ggplot2::ggplot() +
          geom_raster(data=errordf, aes(x=cellularity, y=ploidy, fill=1/error)) +
          geom_point(data=minimadf, aes(x=cellularity, y=ploidy, alpha=min(error)/error), shape=16) +
          scale_fill_gradient(low="green", high="red") +
          ggtitle("Matrix of errors") +
          theme(plot.title = element_text(hjust = 0.5))
  } else {
      tempplot <- ggplot2::ggplot() +
          geom_raster(data=errordf, aes(x=cellularity, y=ploidy, fill=1/error)) +
          scale_fill_gradient(low="green", high="red") +
          ggtitle("Matrix of errors") +
          theme(plot.title = element_text(hjust = 0.5))
  }
  svg(errorgraph_out)
  tempplot
  dev.off()

  
  return(list(method=method, 
              penalty=penalty, 
              penploidy=penploidy,
              standard=standard,
              errormatrix=errormatrix, 
              minimatrix = minimat, 
              errordf=errordf, 
              minimadf=minimadf, 
              matrixplot=tempplot))
  
}


# Function to fix segments
Concatenate_lengths <- function(segmentvalues,segmentdata){
    # if segmentdata is valid do nothing
    if(length(segmentdata$lengths) == length(segmentdata$values)){
        return(segmentdata)
    }else{
        while(length(segmentdata$lengths) != length(segmentdata$values)){
            # find index where segmentvalues are concatenated
            ix <- min(which(!segmentdata$values == segmentvalues))
            # concatenate lengths
            segmentdata$lengths[ix] <- segmentdata$lengths[ix-1] + segmentdata$lengths[ix]
            # Remove redundant length
            segmentdata$lengths <- segmentdata$lengths[-(ix-1)]
        }
        return(segmentdata)
    }
}

Add_normal_segments <- function(segments,chromlengths){
    # intialize empty dataframe
    New_Segments <- data.frame(Chromosome = integer(), Start = integer(), End = integer(), Segment_Mean = numeric())
    # iterate over chromosomes
    for(chr in seq(1,22)){
        Segments_chr <- segments[Segments$Chromosome == chr,c(2,3,4,6)]
                                        # iterate over rows
        for(i in 1:nrow(Segments_chr)){
            if(i == 1){
                start <- 1
                end <- Segments_chr$Start[i] -1
                New_Segments <- rbind(New_Segments, data.frame(Chromosome=chr,Start=start,End=end,Segment_Mean=0))
            }else{
                start <- Segments_chr$End[i-1] + 1
                end <- Segments_chr$Start[i] - 1
                New_Segments <- rbind(New_Segments, data.frame(Chromosome=chr,Start=start,End=end,Segment_Mean=0))
            }
            New_Segments <- rbind(New_Segments,Segments_chr[i,])
        }
        start <- Segments_chr$End[i] +1
        end <- chromlengths[chr,3]
        if(end > start){
            New_Segments <- rbind(New_Segments, data.frame(Chromosome=chr,Start=start,End=end,Segment_Mean=0))
        }
    }
    return(unique(New_Segments))
}

if(exists('ploidy')){
    if(ploidy == "squaremodel"){
        Plot_errors <- function(errordf, minimadf,sample, outfile) {
            p <- ggplot2::ggplot() +
                geom_raster(data=errordf, aes(x=cellularity, y=ploidy, fill=1/error)) +
                geom_point(data=minimadf, aes(x=cellularity, y=ploidy, alpha=min(error)/error), shape=16) +
                scale_fill_gradient(low="green", high="red") +
                ggtitle(sample) +
                theme(plot.title = element_text(hjust = 0.5))
            suppressMessages(ggsave(outfile, plot = p))
        }

    }else{
        Plot_errors <- function(tempdf, minimadf,outfile) {
            p <- ggplot() +
                scale_y_continuous(name = "relative error", limits = c(0,1.05), expand=c(0,0)) +
                scale_x_continuous(name = "cellularity (%)") +
                geom_vline(xintercept = seq(from = 10, to = 100, by = 10), color = "#666666", linetype = "dashed") +
                geom_point(aes(y=errorlist, x=cellularity), data=tempdf) +
                geom_point(aes(y=rerror, x=minima), data=minimadf, color = 'red') +
                theme_classic() + theme(
                                      axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
                ggtitle(paste0(sample, " - errorlist")) +
                theme(plot.title = element_text(hjust = 0.5))
            suppressMessages(ggsave(outfile, plot = p))
        }
    }
}
# Function to run ACE from segmentdata
Run_ACE <- function(segmentdata, Segments, ploidy, sample, method, penalty,penploidy, fit_out, errorgraph_out){
    # calculate standard
    standard <- median(rep(segmentdata$values,segmentdata$lengths))
    # create empty factors
    fraction <- c()
    expected <- c()
    temp <- c()
    errorlist <- c()
    q <- as.integer(gsub("N$","",ploidy))
    #iterate over purities
    for (i in seq(5,100)) {
        fraction[i-4] <- i/100
        for (p in seq(1,12)) {
            expected[p] <- standard*(p*fraction[i-4] + 2*(1-fraction[i-4]))/(fraction[i-4]*q + 2*(1-fraction[i-4]))
        }
        for (j in seq_along(segmentdata$values)) {
            if(method=='RMSE') {temp[j] <- (min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(q-2))^penploidy/(fraction[i-4]^penalty))^2}
            else if(method=='SMRE') {temp[j] <- sqrt(min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(q-2))^penploidy/(fraction[i-4]^penalty))}
            else if(method=='MAE') {temp[j] <- min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(q-2))^penploidy/(fraction[i-4]^penalty)}
            else {print("Not a valid method")}
        }
        if(method=='RMSE') {errorlist[i-4] <- sqrt(sum(temp*segmentdata$lengths)/sum(segmentdata$lengths))}
        else if(method=='SMRE') {errorlist[i-4] <- sum(temp*segmentdata$lengths)/sum(segmentdata$lengths)^2}
        else if(method=='MAE') {errorlist[i-4] <- sum(temp*segmentdata$lengths)/sum(segmentdata$lengths)}
        
    }
    minima <- c()
    rerror <- c()
    
    if (round(errorlist[1], digits = 10) < round(errorlist[2], digits = 10)) {
        lastminimum <- fraction[1]
        minima[1] <- fraction[1]
        rerror[1] <- errorlist[1]/max(errorlist)
    }
    
    for (l in seq(6,99)) {
        if (round(errorlist[l-4], digits = 10) < round(errorlist[l-5], digits = 10) & round(errorlist[l-4], digits = 10) < round(errorlist[l-3], digits = 10)) { 
            lastminimum <- fraction[l-4]
            minima <- append(minima,fraction[l-4])
            rerror <- append(rerror,(errorlist[l-4]/max(errorlist)))
        }
    }
    
    if (errorlist[100-4] <= errorlist[100-5]) {
        lastminimum <- fraction[100-4]
        minima <- append(minima, fraction[100-4])
        rerror <- append(rerror, errorlist[100-4]/max(errorlist))
    }
    chr <- Segments$Chromosome
    bin <- seq_along(chr)
    rlechr <- rle(chr)
    lastchr <- length(rlechr$values)
    binchrend <- c()
    currentbin <- 0
    binchrmdl <- c()
    for (i in seq_along(chr)){
        currentmiddle <- currentbin+bin[i]/2
        currentbin <- currentbin+bin[i]
        binchrend <- append(binchrend, currentbin)
        binchrmdl <- append(binchrmdl, currentmiddle)
    }
    cellularity <- seq(5,100)
    tempdf <- data.frame(cellularity,errorlist=errorlist/max(errorlist))
    minimadf <- data.frame(minima=minima*100,rerror)
    bfi <- tail(which(rerror==min(rerror)))
    fitpicker <- matrix(ncol = 17, nrow = 1)
    colnames(fitpicker) <- c("sample","likely_fit","ploidy","standard","rerror","fit_1","fit_2","fit_3","fit_4","fit_5","fit_6","fit_7","fit_8","fit_9","fit_10","fit_11","fit_12")
    fitpicker[1,1] <- sample
    fitpicker[1,2] <- minima[bfi]
    fitpicker[1,3] <- q
    fitpicker[1,5] <- min(rerror)
    for (m in seq_along(minima)) {
        fitpicker[1,m+5] <- minima[m]
        adjustedsegments <- q + ((2^Segments$Segment_Mean-standard)*(minima[m]*(q-2)+2))/(minima[m]*standard)
        df <- as.data.frame(adjustedsegments[seq_along(bin)])
        colnames(df) <- "segments"
        dfna <- na.exclude(df)
        segments <- dfna$segments
        cappedsegments <- dfna[dfna$segments > 12,]
        if(length(cappedsegments)>0) {cappedsegments <- 12-0.1}
        toppedsegments <- dfna[dfna$segments <= 0,]
        if(length(toppedsegments)>0) {toppedsegments <- 0+0.1}
        line1 <- paste0("Cellularity: ", minima[m])
        line2 <- paste0("Relative error: ", round(rerror[m], digits = 3))
    }
    # Plot errorgraph
    Plot_errors(tempdf, minimadf,errorgraph_out)
    # Write fittable
    write.table(fitpicker, file = fit_out, row.names = F, quote = F, sep = "\t")
}


# Function to run ACE (square_model) from segmentdata
Run_squaremodel <- function(segmentdata, Segments,
                    sample, method, penalty,penploidy, fit_out, errorgraph_out,
                    #prows=4, # include only integer ploidy
                    prows=100, # to include any ploidy
                    ptop =5, pbottom = 1, cellularities = seq(5,100)){
    # calculate standard
    standard <- median(rep(segmentdata$values,segmentdata$lengths))
    # Initialize variables
    fraction <- sort(cellularities)/100
    error <- c()
    errormatrix <- matrix(nrow=(prows+1),ncol=length(fraction))
    listofploidy <- c()
    listofcellularity <- c()
    listoferrors <- c()
    for (t in seq(0,prows)){
        ploidy <- ptop-((ptop-pbottom)/prows)*t
        listofploidy <- append(listofploidy, rep(ploidy,length(fraction)))
        expected <- c()
        hexpected <- c()
        temp <- c()
        errorlist <- c()
        for (i in seq_along(fraction)) {
            for (p in seq(1,12)) {
                expected[p] <- standard*(p*fraction[i] + 2*(1-fraction[i]))/(fraction[i]*ploidy + 2*(1-fraction[i]))
                hexpected[p] <- standard*(p*fraction[i] + 1*(1-fraction[i]))/(fraction[i]*ploidy + 2*(1-fraction[i]))
            }
            for (j in seq_along(segmentdata$values)) {
                if(method=='RMSE') {temp[j] <- (min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty))^2}
                else if(method=='SMRE') {temp[j] <- sqrt(min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty))}
                else if(method=='MAE') {temp[j] <- min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty)}
                else {print("Not a valid method")}
            }
            for (k in seq_along(segmentdata$values)) {
                if(method=='RMSE') {temp[j+k] <- (min(abs(segmentdata$values[k]-hexpected),0.5)/(fraction[i]^penalty))^2}
                else if(method=='SMRE') {temp[j+k] <- sqrt(min(abs(segmentdata$values[k]-hexpected),0.5)/(fraction[i]^penalty))}
                else if(method=='MAE') {temp[j+k] <- min(abs(segmentdata$values[k]-hexpected),0.5)/(fraction[i]^penalty)}
                else {print("Not a valid method")}
            }
            lengths <- c(segmentdata$lengths, segmentdata$lengths)
            if(method=='RMSE') {errorlist[i] <- sqrt(sum(temp*lengths)/sum(lengths))}
            else if(method=='SMRE') {errorlist[i] <- sum(temp*lengths)/sum(lengths)^2}
            else if(method=='MAE') {errorlist[i] <- sum(temp*lengths)/sum(lengths)}
        }
        listofcellularity <- append(listofcellularity, fraction)
        listoferrors <- append(listoferrors, errorlist)
        errormatrix[t+1,] <- errorlist
        
    }
    minimat <- matrix(nrow=(prows+1),ncol=length(fraction))
    for (i in seq(1,prows+1)) {
        for (j in seq(1,length(fraction))) {
            if (i==1|i==(prows+1)) {
                minimat[i,j] <- FALSE
            } else if (j==length(fraction)) {
                minimat[i,j] <- errormatrix[i,j]==min(errormatrix[seq(i-1,i+1),seq(j-1,j)])
            } else {
                minimat[i,j] <- errormatrix[i,j]==min(errormatrix[seq(i-1,i+1),seq(j-1,j+1)])
            }
        }
    }
    errordf <- data.frame(ploidy=listofploidy,
                          cellularity=listofcellularity,
                          error=listoferrors/max(listoferrors),
                          minimum=as.vector(t(minimat)))
    minimadf <- errordf[errordf$minimum==TRUE,]
    minimadf <- minimadf[order(minimadf$error,-minimadf$cellularity),]
    
    # Plot error matrix
    Plot_errors(errordf, minimadf,sample, errorgraph_out)
    # Write fits to table
    write.table(minimadf, file = fit_out, row.names = F, quote = F, sep = "\t")
}

read_ACE_fit <- function(fit){
    fits <- read.table(fit,header=T)
    if(grepl("squaremodel",fit)){
        cellularity <- data.frame(ploidy = fits[1,1], cellularity = fits[1,2], error = fits[1,3])
    }else{
        cellularity <- data.frame(ploidy = fits[1,3], cellularity = fits[1,2], error = fits[1,5])
    }
    return(cellularity)
}

# this function takes a template dataframe as specified in objectsampletotemplate
# you can use a QDNAseq object as "template", then specify QDNAseqobjectsample by the sample number!
# don't forget to specify cellularity, error, ploidy, and standard; you can find these in the model output
# chrsubset restricts the plot to the selected CONTIGUOUS chromosomes: the function takes the lowest and highest 
#  integer from the entered vector and plots the range of chromosomes
#  e.g. chrsubset = c(3,6) and chrsubset = 3:6 result in the same plot
#  subsetting chromosomes messes up the position for the legend, 
#  using chrsubset=1:22 is therefore a "hack" to remove the legend but still get the same plot
# you can scale your upper and lower y-axis limit to the data using cap="max" and bottom="min" respectively
singleplot <- function(template, QDNAseqobjectsample = FALSE, cellularity = 1, error, ploidy = 2, standard, title,
                       trncname = FALSE, cap = 12, bottom = 0, chrsubset, onlyautosomes = TRUE) {
  if(QDNAseqobjectsample) {
    if(missing(title)) {
      pd<-Biobase::pData(template)
      if(trncname==TRUE) {pd$name <- gsub("_.*","",pd$name)}
      if(trncname!=FALSE&&trncname!=TRUE) {pd$name <- gsub(trncname,"",pd$name)}
      title <- pd$name[QDNAseqobjectsample]
    }
    template <- objectsampletotemplate(template, QDNAseqobjectsample)
  }
  if(missing(title)) {title <- "Plot"}
  segmentdata <- rle(as.vector(na.exclude(template$segments)))
  if(missing(standard) || !is(standard, "numeric")) { standard <- median(rep(segmentdata$values,segmentdata$lengths)) }
  adjustedcopynumbers <- ploidy + ((template$copynumbers-standard)*(cellularity*(ploidy-2)+2))/(cellularity*standard)
	adjustedsegments <- ploidy + ((template$segments-standard)*(cellularity*(ploidy-2)+2))/(cellularity*standard)
	template$chr <- gsub("chr","",template$chr,ignore.case = TRUE)
	df <- data.frame(bin=template$bin,chr=template$chr,adjustedcopynumbers,adjustedsegments)
	
	if(onlyautosomes==TRUE) {
	  rlechr <- rle(as.vector(template$chr[template$chr %in% seq(1,22)]))
	}	else if (onlyautosomes==FALSE) {rlechr <- rle(as.vector(template$chr))
	} else {rlechr <- rle(as.vector(template$chr[template$chr %in% seq(1,onlyautosomes)]))}
	binchrend <- c()
	currentbin <- 0
	binchrmdl <- c()
	for (i in seq_along(rlechr$values)) {
		currentmiddle <- currentbin+rlechr$lengths[i]/2
		currentbin <- currentbin+rlechr$lengths[i]
		binchrend <- append(binchrend, currentbin)
		binchrmdl <- append(binchrmdl, currentmiddle)
	} 
	colnames(df)[3] <- "copynumbers"
	colnames(df)[4] <- "segments"
	
	dfna <- na.exclude(df)
	bin <- dfna$bin
	copynumbers <- dfna$copynumbers
	segments <- dfna$segments
	if(bottom=="min"){bottom<-floor(min(dfna$copynumbers))}
	if(cap=="max"){cap<-ceiling(max(dfna$copynumbers))}
	if(!missing(chrsubset)){
	  firstchr <- range(chrsubset)[1]
	  lastchr <- range(chrsubset)[2]
	  dfna <- dfna[dfna$chr %in% rlechr$values[seq(firstchr,lastchr)],]
	}
	cappedcopynumbers <- dfna[dfna$copynumbers > cap,]
	if(length(cappedcopynumbers$copynumbers)>0) {cappedcopynumbers$copynumbers <- cap-0.1}
	cappedsegments <- dfna[dfna$segments > cap,]
	if(length(cappedsegments$segments)>0) {cappedsegments$segments <- cap-0.1}
	toppedcopynumbers <- dfna[dfna$copynumbers <= bottom,]
	if(length(toppedcopynumbers$copynumbers)>0) {toppedcopynumbers$copynumbers <- bottom+0.1}
	toppedsegments <- dfna[dfna$segments <= bottom,]
	if(length(toppedsegments$segments)>0) {toppedsegments$segments <- bottom+0.1}
	

	
	line1 <- paste0("Cellularity: ", cellularity)
	line2 <- paste0("")
	if(!missing(error)) {
	  line2 <- paste0("Relative error: ", round(error, digits = 3))
	}
	if(missing(chrsubset)){
  	ggplot2::ggplot() +
  		scale_y_continuous(name = "copies", limits = c(bottom,cap), breaks = seq(bottom,cap), expand=c(0,0)) +
  		scale_x_continuous(name = "chromosome", limits = c(0,tail(binchrend,1)), breaks = binchrmdl, labels = rlechr$values, expand = c(0,0)) +
  		geom_hline(yintercept = seq(0,4), color = '#333333', size = 0.5) +
  		geom_hline(yintercept = seq(5,(cap-1)), color = 'lightgray', size = 0.5) +
  		geom_vline(xintercept = binchrend, color = "#666666", linetype = "dashed") +
  		geom_point(aes(x = bin,y = copynumbers),data=dfna[dfna$copynumbers>bottom&dfna$copynumbers<cap,], size = 0.1, color = 'black') +
  	  geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'black', shape = 24) +
  	  geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'black', shape = 25) +
  		geom_point(aes(x = bin,y = segments),data=dfna[dfna$segments>bottom&dfna$segments<cap,], size = 1, color = 'darkorange') +
  	  geom_point(aes(x = bin,y = segments),data=cappedsegments, size = 1, color = 'darkorange', shape = 24) +
  	  geom_point(aes(x = bin,y = segments),data=toppedsegments, size = 1, color = 'darkorange', shape = 25) +
	    theme_classic() + theme(
	      axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
	    ggtitle(title) +
  	  theme(plot.title = element_text(hjust = 0.5)) +
  		geom_rect(aes(xmin=binchrmdl[1]/10, xmax = binchrmdl[4], ymin = cap-0.9, ymax = cap), fill = 'white') +
  		annotate("text", x = binchrmdl[2], y = cap-0.3, label = line1) +
  		annotate("text", x = binchrmdl[2], y = cap-0.7, label = line2)
	} else {
	  
	  if(firstchr==1){firstbin<-0
	  } else {firstbin<-binchrend[firstchr-1]+1}
	  ggplot2::ggplot() +
	    scale_y_continuous(name = "copies", limits = c(bottom,cap), breaks = seq(bottom,cap), expand=c(0,0)) +
	    scale_x_continuous(name = "chromosome", limits = c(firstbin,binchrend[lastchr]), 
	                       breaks = binchrmdl[seq(firstchr,lastchr)], 
	                       labels = rlechr$values[seq(firstchr,lastchr)], expand = c(0,0)) +
	    geom_hline(yintercept = seq(0,4), color = '#333333', size = 0.5) +
	    geom_hline(yintercept = seq(5,(cap-1)), color = 'lightgray', size = 0.5) +
	    geom_vline(xintercept = binchrend[seq(firstchr,lastchr)], color = "#666666", linetype = "dashed") +
	    geom_point(aes(x = bin,y = copynumbers),data=dfna[dfna$copynumbers>bottom&dfna$copynumbers<cap,], size = 0.1, color = 'black') +
	    geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'black', shape = 24) +
	    geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'black', shape = 25) +
	    geom_point(aes(x = bin,y = segments),data=dfna[dfna$segments>bottom&dfna$segments<cap,], size = 1, color = 'darkorange') +
	    geom_point(aes(x = bin,y = segments),data=cappedsegments, size = 1, color = 'darkorange', shape = 24) +
	    geom_point(aes(x = bin,y = segments),data=toppedsegments, size = 1, color = 'darkorange', shape = 25) +
	    theme_classic() + theme(
	      axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
	    ggtitle(title) +
	    theme(plot.title = element_text(hjust = 0.5))
	}
}
