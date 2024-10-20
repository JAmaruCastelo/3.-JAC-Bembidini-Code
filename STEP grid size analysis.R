# Packages necessary ------------------------------------------------------
library (terra)
library(biomapME)
library(macrounlp)
library(glue)
library (ggplot2)
library(dplyr)
# list_files and data -----------------------------------------------------
path_ocurrences="D:/00017-Bebidini_Articulo/0. Ocurrencias/Into Study Area/5.- Ocurrences that have points in our study area"
path_save="D:/00017-Bebidini_Articulo/3. Endemism analysis"
path_ocurrences<- list.files(path_ocurrences, full.names=T)

# Create the extension and grid size --------------------------------------
minx=-90
maxx=-30
miny=-60
maxy=20

# Create the matrix for the focal analysis around each point --------------
matriz<-matrix(c(0,1,0,1,1,1,0,1,0), nrow=3)

# Open the area of study --------------------------------------------------
study<-vect(glue("{path_save}/study_area_Bembidini.sqlite"))


# Calculate the proportion of null cells in the sampled area when  --------
secuencia_scale<-seq(0.1,5, 0.1)
values=c()
for (e in secuencia_scale){
  extension= rast(nrows=(maxy-miny)/e,
                  ncols=(maxx-minx)/e,
                  xmin=minx, xmax=maxx, ymin=miny, ymax=maxy)
  s<-sampled_area(path_ocurrences, extension) ## in order to obtained all the values in the sampled area
  s<-focal_w_p(s, matriz, rast=T)
  s<-crop(s, study)
  s<-mask(s, study)
  # caculation of the proportion of area
  selected_cells <- s == 0 # to cells without values
  area_cells <- cellSize(s, unit="km")
  null_area<- sum(area_cells[selected_cells], na.rm = TRUE)
  # calculation of the total area
  total<-mask(area_cells,study)
  total_area<-sum(values(total), na.rm=T)
  # calculation of the proportion
  val=null_area/total_area
  values=c(values, val)
}


# Create the graphic of the analysis  -------------------------------------
area=secuencia_scale*secuencia_scale
info<- data.frame(secu=secuencia_scale, val=values, area=area)

svg(paste(path_save,"/analysis_null.svg", sep=""), 6,5)
ggplot(data=info, aes(x=secu, y=values))+
  geom_point(color="black")+
  geom_smooth()+
  geom_vline(aes(xintercept=0.5))
dev.off()

val_y=info$val[which(info$secu >= 0.5)[1]]

# Montecarlo iteration to see if the grid sample is good ----------------------

# calculation of the presence absence matrix used in the analysis
extension= rast(nrows=(maxy-miny)/0.5,
                ncols=(maxx-minx)/0.5,
                xmin=minx, xmax=maxx, ymin=miny, ymax=maxy)
s<-sampled_area(path_ocurrences, extension)
presence <- Data_form (path_ocurrences, extension, matriz, s=s, map=F)

# make the matrix
matrix<-presence$names[, -ncol(presence$names)]
rownames(matrix)<-presence$names$code

# make the site coords
site_coor<-presence$cells
rownames(site_coor)<-site_coor$code
site_coor<-site_coor[, -ncol(site_coor)]

# Calculate the observed values
richness<-weighted.endemism(matrix,records="site", site.coords=site_coor,
                            deg.resolution=0.5, weight.type= "richness",
                            type="weighted", plot.raster=F)
valor_obs<-richness$WE # this is the we of each cell observed that we used to compare with the monte carlo randomization

# make the iterations with the randomized values
values_df<-data.frame(valor_obs, names=names(valor_obs))
for (e in seq(1, 200)){
  matrix_randomly<-randomize_matrix(matrix)
  rownames(matrix_randomly)<-rownames(site_coor)
  richness2<-weighted.endemism(data.frame(matrix_randomly),records="site", site.coords=site_coor,
                               deg.resolution=0.5, weight.type= "richness",
                               type="weighted", plot.raster=F)
  values<-data.frame(richness2$WE, names=names(richness2$WE))
  values_df<-left_join(values_df, values, by="names")
}
colnames(values_df)<- c("obs","name",seq(1, 200))

# Calculate the p values in montecarlo randomization
valores_nulos<-values_df[seq(3,202)] # the cells with the values od the Monte Carlo randomization

p_value <- valores_nulos >= valor_obs
p_value <- apply(p_value, 1, function(x){sum(x,na.rm=T)/200}) ### apply to different values needed for the analysis

names(p_value)<-names(valor_obs)
s_p_value=richness$WE_raster
s_p_value[as.numeric(names(p_value))] <- p_value
s_p_value<-s_p_value<0.05

writeRaster(s_p_value, glue("{path_save}/p_value_richness.tif"))
