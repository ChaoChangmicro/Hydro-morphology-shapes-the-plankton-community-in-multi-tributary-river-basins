rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
    id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
    net.Raw=netRaw #don't want change netRaw
    net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
    if (abundance.weighted){
        net.stength= net.Raw*sp.ra
    } else {
        net.stength= net.Raw
    }
    
    sp.meanInteration<-colMeans(net.stength)
    
    id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
    remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
    #for simplicity, I only consider the immediate effects of removing the
    #'id.rm' species; not consider the sequential effects of extinction of
    # the 'id.rm2' species.
    
    #you can write out the network pruned
    #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
    # write.csv( net.Raw,"network pruned.csv")
    
    remain.percent
}

rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
    t(sapply(rm.p.list,function(x){
        remains=sapply(1:nperm,function(i){
            rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
        })
        remain.mean=mean(remains)
        remain.sd=sd(remains)
        remain.se=sd(remains)/(nperm^0.5)
        result<-c(remain.mean,remain.sd,remain.se)
        names(result)<-c("remain.mean","remain.sd","remain.se")
        result
    }))
}