vpa.mrm.arg3 <- function(spe,ENV1,GEO1,hydro1){
    spe_dis <- vegdist(spe, method = 'bray')
    
    aa <- match(rownames(spe),rownames(ENV1))
    ENV1 <- ENV1[aa,]
	ENV_1 <- vegdist(ENV1['NO3.N'], method = 'euclidean')
	ENV_2 <- vegdist(ENV1['NO2.N'], method = 'euclidean')
	ENV_3 <- vegdist(ENV1['TDN'], method = 'euclidean')
	ENV_4 <- vegdist(ENV1['TP'], method = 'euclidean')
	ENV_5 <- vegdist(ENV1['TDP'], method = 'euclidean')
	ENV_6 <- vegdist(ENV1['NH4.N'], method = 'euclidean')
	ENV_7 <- vegdist(ENV1['DOC'], method = 'euclidean')
	ENV_8 <- vegdist(ENV1['HIX'], method = 'euclidean')
	ENV_9 <- vegdist(ENV1['DO'], method = 'euclidean')
	ENV_10 <- vegdist(ENV1['ORP'], method = 'euclidean')
	ENV_11 <- vegdist(ENV1['pH'], method = 'euclidean')
	ENV_12 <- vegdist(ENV1['MAT'], method = 'euclidean')
    
    
	
	bb <- match(rownames(spe),rownames(hydro1))
    hydro1 <- hydro1[bb,]
	hydro_1 <- vegdist(hydro1['Velocity'], method = 'euclidean')
	hydro_2 <- vegdist(hydro1['Alt'], method = 'euclidean')
	hydro_3 <- vegdist(hydro1['Temp'], method = 'euclidean')
	hydro_4 <- vegdist(hydro1['Discharge'], method = 'euclidean')
	hydro_5 <- vegdist(hydro1['Width'], method = 'euclidean')
	hydro_6 <- vegdist(hydro1['Agriculture'], method = 'euclidean')
	hydro_7 <- vegdist(hydro1['Forest'], method = 'euclidean')
	hydro_8 <- vegdist(hydro1['Shrub'], method = 'euclidean')
	hydro_9 <- vegdist(hydro1['Grass'], method = 'euclidean')
	hydro_10 <- vegdist(hydro1['Urbanization'], method = 'euclidean')
	
	
  
	
	
	
    rownames(GEO1) <- rownames(ENV1)
    aa <- match(rownames(spe),rownames(GEO1))
    GEO1 <- GEO1[aa,]
	geo_dis <- distm(GEO1[c('Lon', 'Lat')]) / 1000
    geo_dis <- log(geo_dis+1, exp(1))
    rownames(geo_dis) <- rownames(ENV1)
    colnames(geo_dis) <- rownames(ENV1)
    geo_dis <- as.dist(geo_dis)
	



    
    r2part <- varpart4(
    spe_dis,
    list(geo_dis),
    list(ENV_1, ENV_2, ENV_3, ENV_4, ENV_5, ENV_6, ENV_7, ENV_8, ENV_9, ENV_10, ENV_11, ENV_12),
    list(hydro_1, hydro_2, hydro_3, hydro_4, hydro_5, hydro_6, hydro_7, hydro_8, hydro_9, hydro_10)
)
    return(r2part)
}