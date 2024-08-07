vpa.arg2 <- function(ARG_beta,ENV,GEO){
    ARG_dis <- vegdist(ARG_beta)
    
    aa <- match(rownames(ARG_beta),rownames(ENV))
    ENV <- ENV[aa,]
    ENV_dis <- vegdist(ENV,method = "euclidean")
    ENV_pcoa <- pcoa(ENV_dis, correction = "none", rn = NULL)
    
    rownames(GEO) <- rownames(ENV)
    aa <- match(rownames(ARG_beta),rownames(GEO))
    GEO <- GEO[aa,]
    distsite <- dist(GEO[,1:2])
    PCNM.site <- pcnm(distsite)
    
    vpa_ARG <- varpart(ARG_dis, ENV_pcoa$vectors[,1:4], 
                       PCNM.site$vectors[,1:4])
    return(vpa_ARG)
}