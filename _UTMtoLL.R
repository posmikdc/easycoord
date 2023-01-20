################################################################################
################################################################################
################################################################################
###
### Title: UTM to Lat/Long (LL) Coordinate Conversion
###
### Author: Daniel Posmik (posmikdc@uchicago.edu)
###
### Last Updated: January 20th, 2023
###
################################################################################
################################################################################
################################################################################

# Load packages ----------------------------------------------------------------
library(tidyverse)

################################################################################
###  Define Ellipsoids and Their Parameters
################################################################################

# Aggregate equitorial and polar radii for 12 common ellipsoids ----------------
ellipsoid_df <- data.frame(
  ellipsoid_name = c("Airy", "Australian", "Bessel_1841", 
                     "Clarke_1866", "Clarke_1880", "Everest", 
                     "Fischer_1960", "International_1924", "South_American_1969", 
                     "WGS_72", "GRS_80", "WGS_84"), 
  a = c(6377563, 6378160, 6377397,                               #Semi-Major Axis "Equitorial Radius" [m]
        6378206, 6378249, 6377276, 
        6378155, 6378388, 6378160, 
        6378135, 6378137, 6378137),
  b = c(6356257, 6356775, 6356079,                               #Semi-Minor Axis "Polar Radius" [m]
        6356584, 6356516, 6356075, 
        6356773, 6356912, 6356774,
        6356751, 6356752, 6356752)
)

#Data retrieved from an NPS Distributed Learning Module (Coordinates and Maps)
#https://www.oc.nps.edu/oc2902w/c_mtutor/shape/shape5.htm
#Department of Oceanography, Naval Postgraduate School
#Contact: clynch@nps.navy.mil
#Accessed on: December 27th, 2022

# Calculate ellipsoidal parameters ---------------------------------------------
ellipsoid_df %<>%
  mutate(f = (a - b) / a) %>%                                    #Calculate flattening
  mutate(ecc_sq = 2*f - f^2) %>%                                 #Calculate eccentricity squared
  mutate(e1 = (1-sqrt(1-ecc_sq))/(1+sqrt(1-ecc_sq)))             #USGS Professional Paper 1395 (1987), p.47 (3-24)

#Formulas retrieved from Peter H. Dana's website 
#https://foote.geography.uconn.edu/gcraft/notes/coordsys/coordsys_f.html
#The Geographer’s Craft Project, Department of Geography, University of Connecticut
#Contact: pdana@pdana.com; ken.foote@uconn.edu
#Accessed on: December 27th, 2022

#Formula retrieved from the U.S. Geological Survey (USGS) Professional Paper 1395 (1987) (p.62, Figure 11)
#https://pubs.usgs.gov/pp/1395/report.pdf
#Accessed on: January 10th, 2023

# Define additional ellipsoidal parameters -------------------------------------
k0 = 0.9996                                                      #Central meridian scale factor

#Value retrieved from the U.S. Geological Survey (USGS) Professional Paper 1395 (1987) (p.57, ¶3)
#https://pubs.usgs.gov/pp/1395/report.pdf
#Accessed on: December 27th, 2022

################################################################################
###  Define Ellipsoids and Their Parameters
################################################################################

UTMtoLL <- function (ellipsoid_name, UTMNorthing, UTMEasting, UTMZone, UTMBand){ 
  
  # Retrieve ellipsoidal parameters by ellipsoid -------------------------------
  a = ellipsoid_df$a[which(ellipsoid_df$ellipsoid_name           #Retrieve equatorial radius by ellipsoid
                           == ellipsoid_name)]
  ecc_sq = ellipsoid_df$ecc_sq[which(ellipsoid_df$ellipsoid_name #Retrieve eccentricity squared by ellipsoid  
                                     == ellipsoid_name)]
  e1 = (1-sqrt(1-ecc_sq))/(1+sqrt(1-ecc_sq))                     #USGS Professional Paper 1395 (1987), p.47 (3-24)
  
  #Formulas retrieved from the U.S. Geological Survey (USGS) Professional Paper 1395 (1987) (p.62, Figure 11)
  #https://pubs.usgs.gov/pp/1395/report.pdf
  #Accessed on: January 10th, 2023
  
  # Calculate x and y coordinates ----------------------------------------------
  # 500,000m offset for longitude
  x = UTMEasting - 500000.0  
  
  # Define y for both northern and southern hemipshere
  for(i in 1:length(UTMBand)){ 
      if(UTMBand[i] >= 'N'){
        y = UTMNorthing                                          #No adjustment in northern hemipshere
      } else {
        y = UTMNorthing - 10000000                               #10,000,000m adjustment in southern hemisphere
      }
  }
  
  # Calculate origin longitude -------------------------------------------------
  long_origin = (UTMZone - 1)*6 - 180 + 3
  
  #Formulas retrieved from the U.S. Geological Survey (USGS) Professional Paper 1395 (1987) (p.61)
  #https://pubs.usgs.gov/pp/1395/report.pdf
  #Accessed on: January 20th, 2023
  
  #Variables -------------------------------------------------------------------
  eccprime_sq = (ecc_sq^2)/(1 - ecc_sq^2)                        #USGS Professional Paper 1395 (1987), p.61 (8-12)
  M = y / k0                                                     #USGS Professional Paper 1395 (1987), p.63 (8-20)
  mu = M/(a*(1-ecc_sq/4-3*ecc_sq*ecc_sq/                         #USGS Professional Paper 1395 (1987), p.47 (7-19)
               64-5*ecc_sq*ecc_sq*ecc_sq/256))                   
  phi1Rad = (mu + (3*e1/2-27*e1*e1*e1/32)*sin(2*mu)              #USGS Professional Paper 1395 (1987), p.46 (3-26)
             + (21*e1*e1/16-55*e1*e1*e1*e1/32)*sin(4*mu)
             +(151*e1*e1*e1/96)*sin(6*mu))
  phi1 = phi1Rad*(180/pi)
  N1 = a/sqrt(1-ecc_sq*sin(phi1Rad)*sin(phi1Rad))                #USGS Professional Paper 1395 (1987), p.61 (4-20)
  T1 = tan(phi1Rad)*tan(phi1Rad)                                 #USGS Professional Paper 1395 (1987), p.61 (8-13)
  C1 = eccprime_sq*cos(phi1Rad)*cos(phi1Rad)                     #USGS Professional Paper 1395 (1987), p.61 (8-14)
  R1 = a*(1-ecc_sq)/(1-ecc_sq*sin(phi1Rad)*sin(phi1Rad))^1.5     #USGS Professional Paper 1395 (1987), p.64 (8-24)
  D = x/(N1*k0)                                                  #USGS Professional Paper 1395 (1987), p.64 (8-25)
  
  #Formulas retrieved from the U.S. Geological Survey (USGS) Professional Paper 1395 (1987) (p.61)
  #https://pubs.usgs.gov/pp/1395/report.pdf
  #Accessed on: January 20th, 2023
  
  # Calculate lat/long ---------------------------------------------------------
  lat = phi1Rad - (N1*tan(phi1Rad)/R1)*                          #USGS Professional Paper 1395 (1987), p.63 (8-17)
                  (D*D/2-(5+3*T1+10*C1-4*C1*C1-9*eccprime_sq)*D*D*D*D/24 +
                  (61+90*T1+298*C1+45*T1*T1-252*eccprime_sq-3*C1*C1)*
                  D*D*D*D*D*D/720)
  lat = lat*(180/pi)                                             
  
  long = (D-(1+2*T1+C1)*D*D*D/6                                  #USGS Professional Paper 1395 (1987), p.63 (8-18)
          +(5-2*C1+28*T1-3*C1*C1+8*eccprime_sq+24*T1*T1)
          *D*D*D*D*D/120)/cos(phi1Rad)
  long = long_origin + long*(180/pi)                
  
  # Return results -------------------------------------------------------------
  View(cbind(UTMZone, UTMBand, UTMEasting, UTMNorthing, lat, long))
  
}


################################################################################
###  Testing
################################################################################

# Example ----------------------------------------------------------------------
ellipsoid_name <- "WGS_84"
UTMNorthing <- c(4344683, 5569493) 
UTMEasting <- c(706832, 95458)
UTMZone <- c(10, 51)
UTMBand <- c("S","H")

# Execute Function -------------------------------------------------------------
UTMtoLL("WGS_84", UTMNorthing, UTMEasting, UTMZone, UTMBand)
