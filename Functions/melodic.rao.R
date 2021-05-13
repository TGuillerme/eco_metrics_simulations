##MEan DIssimilarity Components##
#samp: community matrix; sites in lines, species in columns
#dis: dissimilarity matrix

melodic.rao <- function(samp,dis,abundance.weighted=TRUE){
	N<-dim(samp)[1]
	rao<-richness<-mpd<-simpson<-numeric(N)
	for (i in 1:N){
		sppInSample<-names(samp[i,samp[i, ]>0])
	    #richness[i]<-length(sppInSample)
		if (length(sppInSample)>1){
      		sample.dis<-dis[sppInSample,sppInSample]
			abund.w<-numeric(length(sppInSample))
			if (abundance.weighted){
				abund.w <- samp[i , sppInSample] / sum(samp[i , sppInSample])
			} 	else{
				abund.w <- rep(1 , length(sppInSample)) / length(sppInSample)
			}	
        	sample.weights <- outer(as.numeric( abund.w) , as.numeric( abund.w))
		  	rao[i] <- sum(sample.weights * sample.dis)
		  	mpd[i] <- weighted.mean(sample.dis[lower.tri(sample.dis)] , sample.weights[lower.tri(sample.weights)])
		  	simpson[i]<-sum(2*sample.weights[lower.tri(sample.weights)])

		}	else {
			rao[i]<- mpd[i]  <- simpson[i]<- 0
		}
	}	
	out<-list(rao=rao, mpd=mpd, simpson=simpson)
	return(out)	
}
