install.packages(c("ggplot2", "devtools", "dplyr", "stringr"))
install.packages(c("maps", "mapdata"))
install.packages("ggmap")

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/maps")

samples<-read.csv("map.csv", header = TRUE, sep = ",")


df<-data.frame(lat=samples$LAT, lon=samples$LON, size=samples$SIZE)

colors<-c(rep("deeppink", 10),
         rep("blue", 13),
         rep("aquamarine", 9),
         "darkviolet",
         "tan",
         "brown",
         "green",
         "chocolate1"
         )

pchs=c(rep(19, 10+13+9), rep(17, 5))

#S. intercalatum = yellow 
#curs - 0.2327267767725841,6.591667471026312
#guin - 14.715459989581978,-17.4691970898632
#marg -17.847388905922173,31.03320066454785
#bovis - -7.768059, 35.686072
#haem 28.307603438234516,30.712324786537692


#map africa
af_lat<-c(-6, 13)
af_lon<-c(1.5, 39.8)

af_map<-get_map(location=c(lon = mean(af_lon), lat = mean(af_lat)), 
                color="color",
                source="google", 
                zoom=3,
                maptype="terrain")
af<-ggmap(af_map)
af<-af + geom_point(data=df, aes(lon, lat), col=colors, cex=4, pch=pchs)
af

#map NE
ne_lat<-c(14.3, 13)
ne_lon<-c(1.5, 2.45)
ne_map <- get_map(location = c(lon = mean(ne_lon), lat = mean(ne_lat)),
                  color="color", 
                  source = "google", 
                  zoom=9, 
                  maptype = "terrain")

ne<-ggmap(ne_map)
ne<-ne + geom_point(data=df, aes(lon, lat), col=colors, cex=4)
ne

#map tz
tz_lat<-c(-6.5, -4.95)
tz_lon<-c(39.25, 39.8)

tz_map <- get_map(location = c(lon = mean(tz_lon), lat = mean(tz_lat)),
                  color="color", 
                  source = "google", 
                  zoom=9, 
                  maptype = "terrain")

tz<-ggmap(tz_map)
tz<-tz + geom_point(data=df, aes(lon, lat), col=colors, cex=4)
tz
