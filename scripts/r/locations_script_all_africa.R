setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/maps")


#------------------#
# Loading packages
#------------------#
library("maps")
library("mapdata")
library("maptools")  #for shapefiles
library("scales")  #for transparency
library("gplots")	# for table


#-----------------#
# loading Dataset 
#-----------------#

mydata <-read.table("Locations_all_africa.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")


#---------#
# Figures
#---------#


#postscript(file="sampling_map_tanzania.eps",paper='special',horizontal=FALSE, width=15, height=10)
svg(file="sampling_map_all_africa_publi.svg")

layout(matrix(c(1,2,3,4,5,3,6,6,3),3,3), heights= c(0.8,0.8,0.4), widths=c(0.8,0.8,0.8))
#layout(matrix(c(1,2,3,4,5,3,6,7,3),3,3), heights= c(0.8,0.8,0.4), widths=c(0.8,0.8,0.8))
#layout(matrix(c(1,2,3,4,4,3,4,4,3,5,6,3),3,4), heights= c(0.8,0.8,0.4), widths=c(0.8,0.4,0.4,0.8))

ocean <- readShapePoly("Ocean_and_lake_shapes/ne_10m_ocean.shp")
ocean_110 <- readShapePoly("Ocean_and_lake_shapes/ne_110m_ocean.shp")
coasts <- readShapeLines("Ocean_and_lake_shapes/ne_10m_coastline.shp")
coasts_110 <- readShapeLines("Ocean_and_lake_shapes/ne_110m_coastline.shp")
lands <- readShapePoly("Ocean_and_lake_shapes/ne_10m_land_scale_rank.shp")
lands <- readShapePoly("Ocean_and_lake_shapes/ne_10m_land.shp")
lands_110 <- readShapePoly("Ocean_and_lake_shapes/ne_110m_land.shp")
lakes <- readShapePoly("Ocean_and_lake_shapes/ne_10m_lakes.shp")
lakes_110 <- readShapePoly("Ocean_and_lake_shapes/ne_110m_lakes.shp")
pluv_lake <- readShapePoly("Ocean_and_lake_shapes/ne_10m_lakes_pluvial.shp")
playas <- readShapePoly("Ocean_and_lake_shapes/ne_10m_playas.shp")				
rivers <- readShapeLines("Ocean_and_lake_shapes/ne_10m_rivers_lake_centerlines.shp")
frontiers <- readShapeLines("Ocean_and_lake_shapes/ne_10m_admin_0_boundary_lines_land.shp")

#---------------------
# Display Africa map
#---------------------

par(mar=c(2,2,2,2))

mylongitude <- c(-42,43)
mylatitude <- c(-20,20)

plot(ocean_110, xlim=mylongitude, ylim=mylatitude, col=alpha("deepskyblue")) 
plot(coasts_110, add=TRUE, xlim=mylongitude, ylim=mylatitude, col=alpha("black"))
plot(lands_110, add=TRUE, xlim=mylongitude, ylim=mylatitude, col=alpha("grey", 0.2))
#plot(rivers, add=TRUE, xlim=mylongitude, ylim=mylatitude, col=alpha("deepskyblue")) 
plot(frontiers, add=TRUE, xlim=mylongitude, ylim=mylatitude, col=alpha("black")) 
plot(lakes_110, add=TRUE, xlim=mylongitude, ylim=mylatitude, col=alpha("lightskyblue")) 
#plot(pluv_lake, add=TRUE, xlim=mylongitude, ylim=mylatitude, col=alpha("lightskyblue")) 
#plot(playas, add=TRUE, xlim=mylongitude, ylim=mylatitude, col=alpha("lightskyblue")) 

#brazil <- map("worldHires", regions="Brazil", plot=FALSE, fill=TRUE)
#senegal <- map("worldHires", regions="Senegal", plot=FALSE, fill=TRUE)
#niger <- map("worldHires", regions="Niger", plot=FALSE, fill=TRUE)
#tanzania <- map("worldHires", regions="Tanzania", plot=FALSE, fill=TRUE)

#par(mar=c(2,2,2,2))
#map("worldHires", xlim=c(-42,-43), ylim=c(-8,10))
#map(brazil, col="brown", fill=TRUE, add=TRUE)
#map(senegal, col="brown", fill=TRUE, add=TRUE)
#map(niger, col="brown", fill=TRUE, add=TRUE)
#map(tanzania, col="brown", fill=TRUE, add=TRUE)

# Add the location on the map
mycol <- rep("chocolate3", nrow(mydata))
mycollist <- mydata[,5] == "Sh"
mycol[mycollist] <- "darkmagenta" 

mycolpch <- rep(19, nrow(mydata))
mycollistpch <- mydata[,5] == "Sh"
mycolpch[mycollistpch] <- 17 

points(mydata$Long, mydata$Lat, pch=mycolpch, col=mycol, cex=1)


# Add the scale
#my.x <- mylongitude[1]
#my.y <- par("usr")[3]+(abs(par("usr")[4])-abs(par("usr")[3]))*0.2

map.axes()
#axis(1, col="white")
map.scale(-40,-28,cex=1)
box(lwd=2)

text(-41.5,-9,"BRAZIL",cex=1,col="red",font=2)
text(-41.5,-11,"(BR)",cex=1,col="red",font=2)
text(-41.50,-18.6,"1",cex=1,col="coral4",font=2)
text(-18,14.5,"SENEGAL",cex=1,col="red",font=2)
#text(-11,20,"MAURITANIA",cex=0.8,col="brown",font=2)
text(10,18,"NIGER",cex=1,col="red",font=2)
#text(13,3,"CAMEROON",cex=0.6,col="brown",font=2)
text(-20,0,"Atlantic sea",cex=1.5,col="white",font=4)
text(36.2,-7.3,"TANZANIA",cex=1,col="red",font=2)
#text(40,-5.5,"Zanzibar archipelago",cex=0.8,col="orange",font=2)
#text(43,-7,"Indian",cex=1.5,col="white",font=4)
#text(43,-9,"Ocean",cex=1.5,col="white",font=4)

#---------------
# zoom on Niger 
#---------------

par(mar=c(2,2,2,2))

nilongitude <- c(1.1,2.5)
nilatitude <- c(13,15)

plot(ocean, xlim=nilongitude, ylim=nilatitude, col=alpha("deepskyblue"))
plot(coasts, add=TRUE, xlim=nilongitude, ylim=nilatitude, col=alpha("black"))
plot(lands, add=TRUE, xlim=nilongitude, ylim=nilatitude, col=alpha("grey", 0.2))
plot(rivers, add=TRUE, xlim=nilongitude, ylim=nilatitude, col=alpha("deepskyblue"))  
plot(frontiers, add=TRUE, xlim=nilongitude,ylim=nilatitude, col=alpha("black"))
plot(lakes, add=TRUE, xlim=nilongitude,ylim=nilatitude, col=alpha("lightskyblue"))
plot(pluv_lake, add=TRUE, xlim=nilongitude,ylim=nilatitude, col=alpha("lightskyblue"))
plot(playas, add=TRUE, xlim=nilongitude,ylim=nilatitude, col=alpha("lightskyblue"))

# Add the location on the map
mycol <- rep("chocolate3", nrow(mydata))
mycollist <- mydata[,5] == "Sh"
mycol[mycollist] <- "darkmagenta" 

mycolpch <- rep(19, nrow(mydata))
mycollistpch <- mydata[,5] == "Sh"
mycolpch[mycollistpch] <- 17 

points(mydata$Long, mydata$Lat, pch=mycolpch, col=mycol, cex=1.8)

mypos <- c(1,1,1,3,3,3,1,1,1,1,3,3,1,1,1,1,1,4,1,1,4,1)
			
map.axes()
map.scale(2.2, 14.9, cex=1)
text(2,14.5,"NIGER (NE)",cex=2,col="brown",font=2)

mycol_nb <- rep("coral4", nrow(mydata))
mycollist_nb <- mydata[,5] == "Sh"
mycol_nb[mycollist_nb] <- "darkorchid4" 

mycex <- 1
text(mydata$Long, mydata$Lat, col=mycol_nb, offset=mycex*0.5, cex=mycex, pos=mypos, font=2)

#----------#
# Legends
#----------#

mycol <- rep("chocolate3", nrow(mydata))
mycollist <- mydata[,5] == "Sh"
mycol[mycollist] <- "darkmagenta" 

myloc <- matrix(ncol=9, nrow=5)
for (i in 1:nrow(mydata)){
	myloc.tmp <- paste(i, ": ", mydata[i,2], " (", mydata[i,6], ")", sep="")
	myloc[[i]]<- myloc.tmp
}

myloc.lg <- ncol(myloc)*nrow(myloc)

if (length(mycol) != myloc.lg) {
	mycol[(length(mycol)+1):myloc.lg] <- NA
}

mycol <- matrix(mycol, ncol=ncol(myloc), nrow=nrow(myloc))

# Layout margins calibration
par(mar=c(0,0,0,0))

textplot(myloc, halign="center", hadj=0, cex=1.2, show.rownames=FALSE, show.colnames=FALSE, col.data=mycol)
Sm <- expression("Sampling site"~ italic("S. mansoni"))
Sh <- expression("Sampling site"~ italic("S. haematobium"))
legend(0.35, 1.05, Sm, pch=19, bty="n",cex=2, col="chocolate3")
legend(0.55, 1.05, Sh, pch=17, bty="n",cex=2, col="darkmagenta")

#hadj=0

#------------------
# zoom on Senegal
#------------------

par(mar=c(2,2,2,2))

selongitude <- c(-17,-15.5)
selatitude <- c(16,16.5)

plot(ocean, xlim=selongitude, ylim=selatitude, col=alpha("deepskyblue"))
plot(coasts, add=TRUE, xlim=selongitude, ylim=selatitude, col=alpha("black"))
plot(lands, add=TRUE, xlim=selongitude, ylim=selatitude, col=alpha("grey", 0.2))
plot(frontiers, add=TRUE, xlim=selongitude,ylim=selatitude, col=alpha("black"))
plot(rivers, add=TRUE, xlim=selongitude, ylim=selatitude, col=alpha("deepskyblue")) 
plot(lakes, add=TRUE, xlim=selongitude,ylim=selatitude, col=alpha("lightskyblue"))
plot(pluv_lake, add=TRUE, xlim=selongitude,ylim=selatitude, col=alpha("lightskyblue"))
plot(playas, add=TRUE, xlim=selongitude,ylim=selatitude, col=alpha("lightskyblue"))

# Add the location on the map
points(mydata$Long, mydata$Lat, pch=19, col="chocolate3", cex=1.6)

map.axes()
map.scale(-16,15.75, cex=1)
text(-16,16,"SENEGAL (SN)",cex=2,col="brown",font=2)
text(-16,16.6,"MAURITANIA",cex=2,col="brown",font=2)
text(-16.8,16.4,"Atlantic",cex=1.5,col="white",font=4)
text(-16.8,16.3,"sea",cex=1.5,col="white",font=4)

mycex <- 1
text(mydata$Long, mydata$Lat, col="coral4", offset=mycex*0.5, cex=mycex, pos=c(1,1), font=2)


#-----------------------
# zoom on Victoria Lake
#-----------------------

par(mar=c(2,2,2,2))

plot(ocean, xlim=c(31.5,34), ylim=c(-2.8,-2.5), col=alpha("deepskyblue"))
plot(coasts, add=TRUE, xlim=c(31.5,34), ylim=c(-2.8,-2.5), col=alpha("black"))
plot(lands, add=TRUE, xlim=c(31.5,34), ylim=c(-2.8,-2.5), col=alpha("grey", 0.2))
plot(rivers, add=TRUE, xlim=c(31.5,34), ylim=c(-2.8,-2.5), col=alpha("deepskyblue"))  
plot(frontiers, add=TRUE, xlim=c(31.5,34),ylim=c(-2.8,-2.5), col=alpha("black"))
plot(lakes, add=TRUE, xlim=c(31.5,34),ylim=c(-2.8,-2.5), col=alpha("lightskyblue"))

# Add the location on the map
points(mydata$Long, mydata$Lat, pch=19, col="chocolate3", cex=1.8)

map.axes()
map.scale(cex=1)
#text(33,1,"Uganda",cex=0.8,col="brown",font=2)
text(33,-3.5,"TANZANIA (TZ)",cex=2,col="brown",font=2)
#text(32.75,-3.1,"TANZANIA",cex=1.5,col="brown",font=2)
text(32.6,-1.7,"Lake Victoria",cex=2,col="white",font=4)
#text(32.75,-1.8,"Lake Victoria",cex=1.5,col="white",font=4)

mypos <- c(1,1,1,3,3,3,1,1,1,1,3,3,1,1,1,1,1,4,1,1,4,1)
mycex <- 1
text(mydata$Long, mydata$Lat, col="coral4", offset=mycex*0.5, cex=mycex, pos=mypos, font=2)

#-------------------------------
# zoom on Zanzibar 
#-------------------------------

par(mar=c(2,2,2,2))

cmlongitude <- c(39.1,39.9)
cmlatitude <- c(-6.5,-4.85)

plot(ocean, xlim=cmlongitude, ylim=cmlatitude, col=alpha("deepskyblue"))
plot(coasts, add=TRUE, xlim=cmlongitude, ylim=cmlatitude, col=alpha("black"))
plot(lands, add=TRUE, xlim=cmlongitude, ylim=cmlatitude, col=alpha("grey", 0.2))
plot(rivers, add=TRUE, xlim=cmlongitude, ylim=cmlatitude, col=alpha("deepskyblue"))
plot(frontiers, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("black"))
plot(lakes, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("lightskyblue"))
plot(pluv_lake, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("lightskyblue"))
plot(playas, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("lightskyblue"))

# Add the location on the map
points(mydata$Long, mydata$Lat, pch=17, col="darkmagenta", cex=1.8)

mypos <- c(1,1,1,3,3,3,1,1,1,1,3,3,1,1,1,1,1,4,1,1,4,1,
			4,1,2,3,1,3,2,2,3,1,4,4,2,2,2,3,1,4,4,1,3,1) #Zanzibar Sh positions
			
map.axes()
map.scale(39.7,-6.5, cex=1)
#text(38.3,-5.5,"TANZANIA",cex=2,col="brown",font=2)
text(39.02,-5,"(TZ)",cex=2,col="brown",font=2)
text(39.4,-5.6,"Zanzibar Archipelago",cex=2,col="coral4",font=3)
text(39.5,-5.2,"Pemba island",cex=1.5,col="brown",font=3)
text(39.5,-6,"Unguja island",cex=1.5,col="brown",font=3)
text(39.8,-5.7,"Indian",cex=1.5,col="white",font=4)
text(39.8,-5.75,"Ocean",cex=1.5,col="white",font=4)
mycex <- 1

text(mydata$Long, mydata$Lat, col="darkorchid", offset=mycex*0.5, cex=mycex, pos=mypos, font=2)

##-------------------------------
## zoom on Zanzibar Pemba island
##-------------------------------

#par(mar=c(2,2,2,2))

#cmlongitude <- c(39.6,39.9)
#cmlatitude <- c(-5.6,-4.85)

#plot(ocean, xlim=cmlongitude, ylim=cmlatitude, col=alpha("deepskyblue"))
#plot(coasts, add=TRUE, xlim=cmlongitude, ylim=cmlatitude, col=alpha("black"))
#plot(lands, add=TRUE, xlim=cmlongitude, ylim=cmlatitude, col=alpha("grey", 0.2))
#plot(rivers, add=TRUE, xlim=cmlongitude, ylim=cmlatitude, col=alpha("deepskyblue"))
#plot(frontiers, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("black"))
#plot(lakes, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("lightskyblue"))
#plot(pluv_lake, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("lightskyblue"))
#plot(playas, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("lightskyblue"))

## Add the location on the map
#points(mydata$Long, mydata$Lat, pch=17, col="darkmagenta", cex=1.8)

#mypos <- c(1,1,1,3,3,3,1,1,1,1,3,3,1,1,1,1,1,4,1,1,4,1,
#			4,1,2,3,1,3,2,2,3,1,4,4,2,2,2,3,1,4,4,1,3,1) #Zanzibar Sh positions
			
#map.axes()
#map.scale(40,-5.55, cex=1)
##text(38.3,-5.5,"TANZANIA",cex=2,col="brown",font=2)
##text(38.5,-5.7,"TZ",cex=2,col="brown",font=2)
#text(39.7,-5.6,"Zanzibar Archipelago",cex=2,col="orange",font=3)
#text(39.5,-5.2,"Pemba island",cex=1.5,col="brown",font=3)
##text(39.5,-6,"Unguja island",cex=1.5,col="white",font=3)
#text(40.1,-5.1,"Indian",cex=1.5,col="white",font=4)
#text(40.1,-5.15,"Ocean",cex=1.5,col="white",font=4)
#mycex <- 1

#text(mydata$Long, mydata$Lat, col="darkorchid", offset=mycex*0.5, cex=mycex, pos=mypos, font=2)

##---------------------------------
## zoom on Zanzibar Unguja island
##---------------------------------

#par(mar=c(2,2,2,2))

#cmlongitude <- c(39.1,39.6)
#cmlatitude <- c(-6.5,-5.6)

#plot(ocean, xlim=cmlongitude, ylim=cmlatitude, col=alpha("deepskyblue"))
#plot(coasts, add=TRUE, xlim=cmlongitude, ylim=cmlatitude, col=alpha("black"))
#plot(lands, add=TRUE, xlim=cmlongitude, ylim=cmlatitude, col=alpha("grey", 0.2))
#plot(rivers, add=TRUE, xlim=cmlongitude, ylim=cmlatitude, col=alpha("deepskyblue"))
#plot(frontiers, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("black"))
#plot(lakes, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("lightskyblue"))
#plot(pluv_lake, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("lightskyblue"))
#plot(playas, add=TRUE, xlim=cmlongitude,ylim=cmlatitude, col=alpha("lightskyblue"))

## Add the location on the map
#points(mydata$Long, mydata$Lat, pch=17, col="darkmagenta", cex=1.8)

#mypos <- c(1,1,1,3,3,3,1,1,1,1,3,3,1,1,1,1,1,4,1,1,4,1,
#			4,1,2,3,1,3,2,2,3,1,4,4,2,2,2,3,1,4,4,1,3,1) #Zanzibar Sh positions
			
#map.axes()
#map.scale(39.6,-5.65,cex=1)
##text(38.3,-5.5,"TANZANIA",cex=2,col="brown",font=2)
#text(38.78,-5.7,"(TZ)",cex=2,col="brown",font=2)
#text(39.3,-5.65,"Zanzibar Archipelago",cex=2,col="orange",font=3)
##text(40.2,-5.2,"Pemba island",cex=1,col="white",font=3)
#text(39,-6.1,"Unguja island",cex=1.5,col="brown",font=3)
#text(39.8,-6,"Indian",cex=1.5,col="white",font=4)
#text(39.8,-6.05,"Ocean",cex=1.5,col="white",font=4)
#mycex <- 1

#text(mydata$Long, mydata$Lat, col="darkorchid", offset=mycex*0.5, cex=mycex, pos=mypos, font=2)
#----------------#
# Saving graphs
#----------------#


dev.off()
