setwd('C:/Users/chaiy/Documents/R/Project/tutorial_viirs/scripts')

#setwd('C:/Users/Administrator/Documents/R/Project/tutorial_viirs/scripts/')

#install.packages("ggmap", dependencies = TRUE)

library(doParallel)
library(foreach)
library(raster)
library(sp)
library(rgdal)
library(ggmap)
library(plotly)

##Set directory path
imagery = "images/SVDNB_npp_20180501-20180531_75N180W_vcmslcfg_v10_c201806061100"

##Obtain a list of TIF files, load in the first file in list
tifs = list.files(imagery,pattern = "\\.tif")
rast <- raster(paste(imagery,"/",tifs[1],sep=""))


##Specify WGS84 as the projection of the raster file
wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
projection(rast) <- CRS(wgs84)

shape_direct <- function(url, shp) {
  library(rgdal)
  temp = tempfile()
  download.file(url, temp) ##download the URL taret to the temp file
  unzip(temp,exdir=getwd()) ##unzip that file
  return(readOGR(paste(shp,".shp",sep=""),shp))
}
msa <- shape_direct(url="http://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_cbsa_20m.zip", 
                    shp= "cb_2014_us_cbsa_20m")
projection(msa) <- CRS(wgs84)

#msa_pop <- read.csv("http://www.census.gov/popest/data/metro/totals/2014/files/CBSA-EST2014-alldata.csv")
msa_pop <- read.csv("https://www2.census.gov/programs-surveys/popest/datasets/2010-2014/metro/totals/cbsa-est2014-alldata.csv")

msa_pop <- msa_pop[msa_pop$LSAD=="Metropolitan Statistical Area",]
msa_pop <- msa_pop[order(msa_pop$POPESTIMATE2014),]
msa_pop$NAME <- as.character(msa_pop$NAME) 

# cities <- c("New York, NY", "Los Angeles, CA","Chicago, IL", "Houston, TX",
#             "Philadelphia, PA", "Phoenix, AZ", "San Antonio, TX", "San Diego, CA",     
#             "Dallas, TX", "San Jose, CA", "Austin, TX", "Jacksonville, FL",
#             "San Francisco, CA", "Indianapolis, IN", "Columbus, OH", "Fort Worth, TX",
#             "Charlotte, NC", "Detroit, MI", "El Paso, TX", "Seattle, WA",
#             "Denver, CO","Washington, DC", "Memphis, TN", "Boston, MA",
#             "Nashville, TN", "Baltimore, MD", "Oklahoma City, OK", "Portland, OR",
#             "Las Vegas, NV", "Louisville, KY","Milwaukee, WI","Albuquerque, NM",
#             "Tucson, AZ","Fresno, CA","Sacramento, CA")

cities <- c("New York, NY", "Los Angeles, CA","Chicago, IL", "Houston, TX")

##Set graph layout
par(mai=c(0,0,0,0),mfrow = c(7,5),bg='#001a4d', bty='n')

##Loop through data
coords <- data.frame() ##place holder

for(i in 1:length(cities)){
  
  ##Coords
  temp_coord <- geocode(cities[i], source = "google")
  coords <- rbind(coords,temp_coord)
  
  e <- extent(temp_coord$lon - 1, temp_coord$lon + 1,
              temp_coord$lat - 0.25, temp_coord$lat + 0.25)
  rc <- crop(rast, e)    
  
  ##Rescale brackets
  sampled <- as.vector(rc)
  clusters <- 15
  clust <- kmeans(sampled,clusters)$cluster
  combined <- as.data.frame(cbind(sampled,clust))
  brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])
  
  #Plots
  plot(rc, breaks=brk, col=colorRampPalette(c("#001a4d","#0066FF", "yellow"))(clusters), 
       legend=F,yaxt='n',xaxt='n',frame = F, asp=1.5)
  text(temp_coord$lon ,temp_coord$lat + 0.15,
       substr(cities[i],1,regexpr(",",cities[i])-1), 
       col="white", cex=1.25)
  
  rm(combined)
}



masq <- function(shp,rast,i){
  
  #Extract one polygon based on index value i
  polygon <- shp[i,] #extract one polygon
  extent <- extent(polygon) #extract the polygon extent 
  
  #Raster extract
  outer <- crop(rast, extent) #extract raster by polygon extent
  inner <- mask(outer,polygon) #keeps values from raster extract that are within polygon
  
  #Convert cropped raster into a vector
  #Specify coordinates
  coords <- expand.grid(seq(extent@xmin,extent@xmax,(extent@xmax-extent@xmin)/(ncol(inner)-1)),
                        seq(extent@ymin,extent@ymax,(extent@ymax-extent@ymin)/(nrow(inner)-1)))
  #Convert raster into vector
  data <- as.vector(inner)
  
  #package data in neat dataframe
  data <- cbind(as.character(shp@data$CBSAFP[i]),coords, data) 
  colnames(data)<-c("GEOID","lon","lat","avg_rad") #note that 
  data <- data[!is.na(data$avg_rad),] #keep non-NA values only
  
  return(data)
}



##MSAs by GEOID
msa_list <- c(16180,19140,45820,42540,35620)

##Placeholder
radiances <- data.frame() 

##Loop MSA file
for(i in msa_list){
  
  print(i)
  
  #Extract MSA i polygon
  shp_temp <- msa[msa@data$GEOID==i,]
  
  #Extract MSA abbreviated name
  if(regexpr("-",as.character(shp_temp@data$NAME)[1])[1]==-1){
    loc = as.character(substr(as.character(shp_temp@data$NAME)[1],1,regexpr(",",as.character(shp_temp@data$NAME)[1])-1))
  } else{
    loc = as.character(substr(as.character(shp_temp@data$NAME)[1],1,regexpr("-",as.character(shp_temp@data$NAME)[1])-1))
  }
  
  #Extract the radiances, append to radiances placeholder
  rad <- masq(shp_temp,rast,1)$avg_rad 
  temp <- data.frame(loc = as.character(paste(loc,"(TNL = ",round(sum(rad),0),")",sep="")), avg_rad = rad) 
  radiances <- rbind(radiances,temp)
}

#Use ggplot to create histograms by MSA group. Preload.
ggplot(radiances, aes(x=log(avg_rad))) +
  geom_histogram(position="identity", alpha=0.4) +
  facet_grid(. ~ loc)

#Remove all axes labels for style
x <- list(
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)
y <- list(
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
) 

#Initiate a plotly graph without axes
ggplotly()  %>% layout(xaxis=x, yaxis=y)




##Set up comparisons
registerDoParallel(cores=2)
extract <- foreach(i=1:nrow(msa@data), .combine=rbind) %dopar% {
  data <- masq(msa,rast,i)
  data.frame(GEOID = data$GEOID[1],sum = sum(data$avg_rad))
}
extract$GEOID <- as.numeric(as.character(extract$GEOID))

##Join in data
joined<-merge(extract, msa_pop[,c("CBSA","NAME","POPESTIMATE2014")],by.x="GEOID",by.y="CBSA")

colnames(joined) <- c("GEOID","TNL","MSA","Population")

plot_ly(joined, 
        x = log(TNL), 
        y = log(Population), 
        text = paste("MSA: ", MSA),
        mode = "markers", 
        color = TNL,colors="PuOr")  %>% 
  layout(title="Total Nighttime Light vs. Population", showlegend = F)

