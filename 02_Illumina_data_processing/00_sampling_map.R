# Developed here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Manuscripts\Mya_genome\Figures\intro\mya_map.Rmd

# LOAD PACKAGES
# install.packages("mapdata")
# install.packages("maptools")
# install.packages("rnaturalearth")
# install.packages("raster")
# install.packages("rnaturalearthhires")
# library(raster)
# library(rnaturalearth)
library(mapdata)
library(maptools)
library(tidyverse)

# DATA FOR BLANK MAP
states <- map_data("state")
usa <- map_data("usa")
canada <- map_data("worldHires", "Canada")

# DATA TO PLOT
samples <- read.delim("C:/Users/shart/Metzger Lab Dropbox/Sam_data/Manuscripts/Mya_genome/Figures/intro/sequenced_samples.txt") %>% print()
btn <- read.delim("C:/Users/shart/Metzger Lab Dropbox/Sam_data/Manuscripts/Mya_genome/Figures/intro/Mya_Sampling_Locations.txt") %>% 
  filter(BTN == "Yes") %>% print()
nobtn <- read.delim("C:/Users/shart/Metzger Lab Dropbox/Sam_data/Manuscripts/Mya_genome/Figures/intro/Mya_Sampling_Locations.txt") %>% 
  filter(BTN == "No") %>% print()
PEI2 <- filter(samples, Sample == "PEI-DF488")
NY <- filter(samples, Sample == "NYTC-C9")
MELC <- filter(samples, Sample == "MELC-2E11")
FFM <- filter(samples, Sample == "FFM-19G1")
FFM2 <- filter(samples, Sample == "FFM-22F10")


# CREATE SCALE BAR
## I have taken this more or less directly from:
#http://egallic.fr/scale-bar-and-north-arrow-on-a-ggplot2-map/

#
# Result #
#--------#
# Return a list whose elements are :
#   - rectangle : a data.frame containing the coordinates to draw the first rectangle ;
#   - rectangle2 : a data.frame containing the coordinates to draw the second rectangle ;
#   - legend : a data.frame containing the coordinates of the legend texts, and the texts as well.
#
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distanceLon : length of each rectangle ;
# distanceLat : width of each rectangle ;
# distanceLegend : distance between rectangles and legend texts ;
# dist.units : units of distance "km" (kilometers) (default), "nm" (nautical miles), "mi" (statute miles).
createScaleBar <- function(lon,lat,distanceLon,distanceLat,
                           distanceLegend, dist.units = "km"){
    # First rectangle
    bottomRight <- gcDestination(lon = lon, lat = lat, bearing = 90, 
                                 dist = distanceLon, dist.units = dist.units,
                                 model = "WGS84")
    
    topLeft <- gcDestination(lon = lon, lat = lat, bearing = 0, 
                             dist = distanceLat, dist.units = dist.units, 
                             model = "WGS84")
    rectangle <- cbind(lon=c(lon, lon, bottomRight[1,"long"],
                             bottomRight[1,"long"], lon),
                       lat = c(lat, topLeft[1,"lat"], topLeft[1,"lat"],
                               lat, lat))
    rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)
    
    # Second rectangle t right of the first rectangle
    bottomRight2 <- gcDestination(lon = lon, lat = lat, bearing = 90, 
                                  dist = distanceLon*2, dist.units = dist.units,
                                  model = "WGS84")
    rectangle2 <- cbind(lon = c(bottomRight[1,"long"], bottomRight[1,"long"],
                                bottomRight2[1,"long"], bottomRight2[1,"long"],
                                bottomRight[1,"long"]),
                        lat=c(lat, topLeft[1,"lat"], topLeft[1,"lat"], 
                              lat, lat))
    rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)
    
    # Now let's deal with the text
    onTop <- gcDestination(lon = lon, lat = lat, bearing = 0, 
                           dist = distanceLegend, dist.units = dist.units, 
                           model = "WGS84")
    onTop2 <- onTop3 <- onTop
    onTop2[1,"long"] <- bottomRight[1,"long"]
    onTop3[1,"long"] <- bottomRight2[1,"long"]
    
    legend <- rbind(onTop, onTop2, onTop3)
    legend <- data.frame(cbind(legend, text = c(0, distanceLon, distanceLon*2)),
                         stringsAsFactors = FALSE, row.names = NULL)
    return(list(rectangle = rectangle, rectangle2 = rectangle2, 
                legend = legend))
}


#
# Result #
#--------#
# This function enables to draw a scale bar on a ggplot object, and optionally an orientation arrow #
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distanceLon : length of each rectangle ;
# distanceLat : width of each rectangle ;
# distanceLegend : distance between rectangles and legend texts ;
# dist.units : units of distance "km" (kilometers) (by default), "nm" (nautical miles), "mi" (statute miles) ;
# rec.fill, rec2.fill : filling colour of the rectangles (default to white, and black, resp.);
# rec.colour, rec2.colour : colour of the rectangles (default to black for both);
# legend.colour : legend colour (default to black);
# legend.size : legend size (default to 3);
# orientation : (boolean) if TRUE (default), adds an orientation arrow to the plot ;
# arrow.length : length of the arrow (default to 500 km) ;
# arrow.distance : distance between the scale bar and the bottom of the arrow (default to 300 km) ;
# arrow.North.size : size of the "N" letter (default to 6).
scaleBar <- function(lon, lat, distanceLon, distanceLat, 
                     distanceLegend, dist.unit = "km", rec.fill = "white",
                     rec.colour = "black", rec2.fill = "black", 
                     rec2.colour = "black", legend.colour = "black", 
                     legend.size = 3, orientation = TRUE, arrow.length = 500,
                     arrow.distance = 300, arrow.North.size = 6){
    laScaleBar <- createScaleBar(lon = lon, lat = lat, 
                                 distanceLon = distanceLon, 
                                 distanceLat = distanceLat, 
                                 distanceLegend = distanceLegend, 
                                 dist.unit = dist.unit)
    # First rectangle
    rectangle1 <- geom_polygon(data = laScaleBar$rectangle, 
                               aes(x = lon, y = lat), fill = rec.fill, 
                               colour = rec.colour)
    
    # Second rectangle
    rectangle2 <- geom_polygon(data = laScaleBar$rectangle2, 
                               aes(x = lon, y = lat), fill = rec2.fill, 
                               colour = rec2.colour)
    
    # Legend
    scaleBarLegend <- annotate("text", label = paste(laScaleBar$legend[,"text"],
                                                     dist.unit, sep=""), 
                               x = laScaleBar$legend[,"long"], 
                               y = laScaleBar$legend[,"lat"], 
                               size = legend.size, 
                               colour = legend.colour, fontface="bold")
    
    res <- list(rectangle1, rectangle2, scaleBarLegend)
    
    if(orientation){# Add an arrow pointing North
        coordsArrow <- createOrientationArrow(scaleBar = laScaleBar, 
                                              length = arrow.length, 
                                              distance = arrow.distance,
                                              dist.unit = dist.unit)
        arrow <- list(geom_segment(data = coordsArrow$res, 
                                   aes(x = x, y = y, xend = xend, yend = yend)),
                      annotate("text", label = "N", 
                               x = coordsArrow$coordsN[1,"x"], 
                               y = coordsArrow$coordsN[1,"y"], 
                               size = arrow.North.size, colour = "black"))
        res <- c(res, arrow)
    }
    return(res)
}


# COMPILE AND PRINT MAP
MAP <- ggplot() +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), 
                 fill = "gray80", color="black") + 
    geom_polygon(data = states, aes(x=long, y = lat, group = group), 
                fill = "gray90", color="black") +     #"white" gray90 antiquewhite
# RACHAELS OBS DATA
    geom_point(data=btn, aes(x=long, y=lat), 
               fill="black", color = "black", shape=4, size=1, stroke = 2) + 
    geom_point(data=samples, aes(x=long, y=lat), 
               fill="black", color = "black", shape=4, size=1, stroke = 2) + 
    # geom_point(data=nobtn, aes(x=long, y=lat), 
    #            fill="white", color = "black", shape=21, size=4) + 
# SEQUENCED SAMPLES
# NY
    geom_segment(data=NY, aes(x = long, y = lat-0.75, xend = long, yend = lat),
                 color = "black", size = 1.5) + 
    geom_point(data=NY, aes(x=long, y=lat-0.75), 
               fill="blue", color = "black", shape=21, size=7) +
    geom_text(data=NY, aes(x=long, y=lat-0.75, label = "5"), color = "white", size = 4,fontface = "bold") +
# FFM
    geom_segment(data=FFM, aes(x = long, y = lat-0.75, xend = long, yend = lat),
                 color = "black", size = 1.5) + 
    geom_point(data=FFM, aes(x=long, y=lat-0.75), 
               fill="blue", color = "black", shape=21, size=7) + 
    geom_text(data=FFM, aes(x=long, y=lat-0.75, label = "1"), color = "white", size = 4,fontface = "bold") +
    geom_segment(data=FFM, aes(x = long+0.7, y = lat-0.75, xend = long, yend = lat),
                 color = "black", size = 1.5) + 
    geom_point(data=FFM, aes(x=long+0.7, y=lat-0.75), 
               fill="blue", color = "black", shape=21, size=7) + 
    geom_text(data=FFM, aes(x=long+0.7, y=lat-0.75, label = "2"), color = "white", size = 4,fontface = "bold") +
# FFM2
    geom_segment(data=FFM2, aes(x = long+0.7, y = lat-0.75, xend = long, yend = lat),
                 color = "black", size = 1.5) + 
    geom_point(data=FFM2, aes(x=long+0.7, y=lat-0.75), 
               fill="blue", color = "black", shape=21, size=7) + 
    geom_text(data=FFM2, aes(x=long+0.7, y=lat-0.75, label = "3"), color = "white", size = 4,fontface = "bold") +
# MELC
    geom_segment(data=MELC, aes(x = long-0.7, y = lat-0.75, xend = long, yend = lat),
                 color = "black", size = 1.5) + #, arrow=arrow(length=unit(0.10,"cm"), type = "closed"), size = 1)+
    geom_point(data=MELC, aes(x=long-0.7, y=lat-0.75), 
               fill="blue", color = "black", shape=21, size=7) +
    geom_text(data=MELC, aes(x=long-0.7, y=lat-0.75, label = "4"), color = "white", size = 4,fontface = "bold") +
    geom_segment(data=MELC, aes(x = long, y = lat-0.75, xend = long, yend = lat),
                 color = "black", size = 1.5) +
    geom_point(data=MELC, aes(x=long, y=lat-0.75), 
               fill="black", color = "black", shape=21, size=7) +
    geom_segment(data=MELC, aes(x = long+0.7, y = lat-0.75, xend = long, yend = lat),
                 color = "black", size = 1.5) +
    geom_point(data=MELC, aes(x=long+0.7, y=lat-0.75), 
               fill="black", color = "black", shape=21, size=7) +
    geom_point(data=MELC, aes(x=long+0.7, y=lat-0.75), 
               fill="white", color = "white", shape=8, size=2, stroke = 2) +
  # PEI
    geom_segment(data=PEI2, aes(x = long+1.05, y = lat-0.75, xend = long, yend = lat),
                 color = "black", size = 1.5) +
    geom_point(data=PEI2, aes(x=long+1.05, y=lat-0.75), 
               fill="red", color = "black", shape=21, size=7) + 
    geom_text(data=PEI2, aes(x=long+1.05, y=lat-0.75, label = "3"), color = "white", size = 4,fontface = "bold") +
    geom_segment(data=PEI2, aes(x = long+0.35, y = lat-0.75, xend = long, yend = lat),
                 color = "black", size = 1.5) +
    geom_point(data=PEI2, aes(x=long+0.35, y=lat-0.75), 
               fill="red", color = "black", shape=21, size=7) +
    geom_text(data=PEI2, aes(x=long+0.35, y=lat-0.75, label = "2"), color = "white", size = 4,fontface = "bold") +
    geom_segment(data=PEI2, aes(x = long-0.35, y = lat-0.75, xend = long, yend = lat),
                 color = "black", size = 1.5) +
    geom_point(data=PEI2, aes(x=long-0.35, y=lat-0.75), 
               fill="red", color = "black", shape=21, size=7) + 
    geom_text(data=PEI2, aes(x=long-0.35, y=lat-0.75, label = "1"), color = "white", size = 4,fontface = "bold") +
    geom_segment(data=PEI2, aes(x = long-1.05, y = lat-0.75, xend = long, yend = lat),
                 color = "black", size = 1.5) +
    geom_point(data=PEI2, aes(x=long-1.05, y=lat-0.75), 
               fill="black", color = "black", shape=21, size=7) + 
# LABELS
    geom_text(aes(-63.4, 48, label = "Prince"), size = 4,fontface = "bold") +
    geom_text(aes(-63.4, 47.6, label = "Edward"), size = 4,fontface = "bold") +
    geom_text(aes(-63.4, 47.2, label = "Island (PEI)"), size = 4, fontface = "bold") +
    geom_text(aes(-69, 45.5, label = "Maine"), size = 4, fontface = "bold") +
    geom_text(aes(-75, 43, label = "New York"), size = 4, fontface = "bold") +
    geom_text(aes(-74.5, 44.5, label = "USA"), size = 6, fontface = "bold") +
    geom_text(aes(-74.5, 45.5, label = "Canada"), size = 6, fontface = "bold") +
# KEY
    #geom_text(aes(-65, 40.8, label = "Observed"), size = 4, fontface = "bold", hjust="left") +
    geom_text(aes(-69, 39.0, label = "Disseminated neoplasia observed"), size = 4, fontface = "bold", hjust="left") +
    geom_point(aes(-69.5, 39.0), fill="black", color = "black", shape=4, size=1, stroke = 2) +
    #geom_text(aes(-65, 39.4, label = "Sequenced"), size = 4, fontface = "bold", hjust="left") +
    geom_text(aes(-69.975, 38.5, label = "PEI"), size = 4, color = "red", fontface = "bold") +
    geom_text(aes(-69.025, 38.5, label = "USA"), size = 4, color = "blue", fontface = "bold") +
    geom_text(aes(-68.625, 38, label = "Cancer sequenced"), size = 4, fontface = "bold", hjust="left") +
    geom_point(aes(-69.125, 38), fill="blue", color = "black", shape=21, size=7) +
    geom_point(aes(-69.875, 38), fill="red", color = "black", shape=21, size=7) +
    geom_text(aes(-69, 40.4, label = "Healthy clam sequenced"), size = 4, fontface = "bold", hjust="left") +
    geom_point(aes(-69.5, 40.4), fill="black", color = "black", shape=21, size=7) +
    #geom_point(aes(-70.25, 40.4), fill="white", color = "red", shape=21, size=7) +
    geom_text(aes(-69, 39.7, label = "Healthy clam reference genome"), size = 4, fontface = "bold", hjust="left") +
    #geom_point(aes(-69.5, 39.7), fill="white", color = "black", shape=21, size=7) +
    geom_point(aes(-69.5, 39.7), fill="black", color = "black", shape=21, size=7) +
    geom_point(aes(-69.5, 39.7), fill="white", color = "white", shape=8, size=2, stroke = 2) +
# SCALE BAR
    coord_fixed(xlim = c(-77, -62),  ylim = c(38, 48), ratio = 1.2) +
    theme(line = element_blank(),
              text = element_blank(), 
              panel.background = element_rect(fill = "aliceblue")) + # steelblue
  scaleBar(lon = -77, lat = 47, distanceLon = 200, # lon = -72, lat = 38
            distanceLat = 50, distanceLegend = 100, 
            dist.unit = "km", legend.size = 4, 
            orientation = FALSE)
MAP
setwd("C:/Users/shart/Metzger Lab Dropbox/Sam_data/Manuscripts/Mya_genome/Figures_new")
pdf("mya_map_allsamples_final.pdf", height=6,width=6)
MAP
dev.off()