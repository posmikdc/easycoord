# Easy Coordinate Conversion (EasyCoord) Package

The `EasyCoord` package facilitates the conversion of coordinate systems in R. Previously, the conversion of coordinate systems was only possible in GIS software or through subscription-based third party entities. This package enables an easy-to-use and intuitive alternative for R users. 

The 'EasyCoord' package is a work in progress and submission to CRAN is forthcoming. All code shared in the repository is preliminary. 

This package enables coordinate conversion between the following popular coordinate systems:
* Decimal Degrees (LatLong)
* Universal Transverse Mercator (UTM)
* United States National Grid (USNG)/ Military Grid Reference System (MGRS)
* Global Area Reference System (GARS)
* Geographic Reference System (GEOREF)
* Plus Code

For the most common coordinate systems, there exist direct conversion functions (e.g., Lat/Long <> UTM). However, less popular coordinate systems may require conversion via a more popular intermediary system (e.g., USNG <> UTM <> Lat/Long). More coordinate systems may be added to the package depending on demand and feasibility. 

This code is written as lightweight and computationally efficient. The 'EasyCoord' package seeks to address the need to convert a large amount of data. Moreover, the code is structured intuitively. Inputs are segmented into idiosyncratic input parameters. This may require the user to conduct pre-emptive wrangling to shape the data into the correct format. 

Each common conversion pair comes with a function. For instance, the `_LLtoUTM.R` file contains the function converting Decimal Degrees to UTM and the `_UTMtoLL.R` file converts in the opposite order. 

Let us take a look at this function:

## Example 1: The Function Shell

The [`EasyCoord` package](https://github.com/posmikdc/easycoord) seeks to make dealing with intricate spatial data easy. To avoid inaccuracy, misspecification, or errors, each function is designed to receive each component of the coordinate system separately.

```{r}
library(EasyCoord)

# The shell of the LL to UTM function
LLtoUTM <- function (ellipsoid_name, lat, long) { 
                    ...
                    }
# The shell of the UTM to LL function
UTMtoLL <- function (ellipsoid_name, UTMNorthing, UTMEasting, UTMZone, UTMBand){ 
                    ...
                    }
```

Note how both function require different input parameters due to the structure of the respective input coordinate system. Aside from the coordinate-specific inputs, all functions require the specification of the ellipsoid. While almost all modern U.S.-based coordinate systems rely on the WGS_84 datum, this may vary for older data and for different regions. The WGS_84 datum can be thought of as the basis for the three-dimensional projection of any coordinate system.  

The function body is mainly derived from the U.S. Geological Survey, such as [Professional Paper 1395](https://pubs.er.usgs.gov/publication/pp1395). In that sense, most of the code is derived from the equations set forth in these key pieces of literature. 

## Example 2: The Function Body

As a next example, the mathematical mechanics of the conversion process are detailed. The following code is taken from the `_LLtoUTM.R` file. 

```{r}
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
  mutate(inv_f = 1 / f) %>%                                      #Calculate inverse flattening
  mutate(ecc_sq = 2*f - f^2)                                     #Calculate eccentricity squared

#Formulas retrieved from Peter H. Dana's website 
#https://foote.geography.uconn.edu/gcraft/notes/coordsys/coordsys_f.html
#The Geographer’s Craft Project, Department of Geography, University of Connecticut
#Contact: pdana@pdana.com; ken.foote@uconn.edu
#Accessed on: December 27th, 2022

# Define additional ellipsoidal parameters -------------------------------------
k0 = 0.9996                                                      #Central meridian scale factor

#Value retrieved from the U.S. Geological Survey (USGS) Professional Paper 1395 (1987) (p.57, ¶3)
#https://pubs.usgs.gov/pp/1395/report.pdf
#Accessed on: December 27th, 2022

################################################################################
###  Define Ellipsoids and Their Parameters
################################################################################
  
LLtoUTM <- function (ellipsoid_name, lat, long){ 
  
  # Retrieve ellipsoidal parameters by ellipsoid -------------------------------
  a = ellipsoid_df$a[which(ellipsoid_df$ellipsoid_name           #Retrieve equatorial radius by ellipsoid
                           == ellipsoid_name)]
  ecc_sq = ellipsoid_df$ecc_sq[which(ellipsoid_df$ellipsoid_name #Retrieve eccentricity squared by ellipsoid  
                                     == ellipsoid_name)]
  
  # Convert to Radians ---------------------------------------------------------
  lat_rad <- lat*pi/180 
  long_rad <- long*pi/180 
  
  # Calculate Longitudinal Projection Zones ("Zones") --------------------------
  UTMZone<- NULL                                                 #Empty vector shell
  
  for(i in 1:length(long)){ 
    if(lat[i] >= 72.0 && lat[i] < 84.0){                         #Defintion of Svalbard zones
        #There are four irregular zones near Svalbard
        if(long[i] >= 0.0 && long[i] < 9.0){                     #Zone 31
          UTMZone[i]<- "31"                    
        } else if(long[i] >= 9.0 && long[i] < 21.0){             #Zone 33
          UTMZone[i]<- "33" 
        } else if(long[i] >= 21.0 && long[i] < 33.0){            #Zone 35
          UTMZone[i]<- "35" 
        } else if(long[i] >= 33.0 && long[i] < 42.0){            #Zone 37
          UTMZone[i]<- "37" 
        }
    }
        #There is one irregular zone for Norway
    else if(lat[i] >= 56.0 && lat[i] < 64.0 &&                   #Definition of Norway Zone
            long[i] >= 3.0 && long[i] < 12.0){
          UTMZone[i]<- "32"                                      #Zone 32
        }
        #Calculation of regular zones
    else{
          UTMZone<- floor(((long + 180)/6) %% 60) + 1            #Calculation and standardization
    }
  }
  
  #Convert results to numeric strings
  UTMZone <- as.numeric(UTMZone)
  
  #Data retrieved from the U.S. Geological Survey (USGS) Professional Paper 1395 (1987) (p.62, Figure 11)
  #https://pubs.usgs.gov/pp/1395/report.pdf
  #Accessed on: January 10th, 2023
  
  
  
  # Define Formulas for Calculation --------------------------------------------
  #Calculate the central meridian of a 6 degree wide UTM zone
  long_origin = (UTMZone - 1)*6 - 180 + 3                        
  long_origin_rad = long_origin*(pi/180)
  
  #Variables
  eccprime_sq = (ecc_sq^2)/(1 - ecc_sq^2)                        #USGS Professional Paper 1395 (1987), p.61 (8-12)
  N = a/sqrt(1-ecc_sq*sin(lat_rad)*sin(lat_rad))                 #USGS Professional Paper 1395 (1987), p.61 (4-20)
  T = tan(lat_rad)*tan(lat_rad)                                  #USGS Professional Paper 1395 (1987), p.61 (8-13)
  C = eccprime_sq*cos(lat_rad)*cos(lat_rad)                      #USGS Professional Paper 1395 (1987), p.61 (8-14)
  A = cos(lat_rad)*(long_rad - long_origin_rad)                  #USGS Professional Paper 1395 (1987), p.61 (8-15)
  
  #Calculate the true distance along the central meridian from the Equator to lat_rad
  M = a*((1                                                      #USGS Professional Paper 1395 (1987), p.61 (3-21)
          - ecc_sq/4
          - 3*ecc_sq*ecc_sq/64
          - 5*ecc_sq*ecc_sq*ecc_sq/256)*lat_rad 
          - (3*ecc_sq/8
          + 3*ecc_sq*ecc_sq/32
          + 45*ecc_sq*ecc_sq*ecc_sq/1024)*sin(2*lat_rad)
          + (15*ecc_sq*ecc_sq/256 + 45*ecc_sq*ecc_sq*ecc_sq/1024)*sin(4*lat_rad) 
          - (35*ecc_sq*ecc_sq*ecc_sq/3072)*sin(6*lat_rad))
  
  #Formulas retrieved from the U.S. Geological Survey (USGS) Professional Paper 1395 (1987) (p.61)
  #https://pubs.usgs.gov/pp/1395/report.pdf
  #Accessed on: January 10th, 2023
  
  # Calculate UTM Coordinates --------------------------------------------------
  #Calculate Easting Coordinates
  UTMEasting = (k0*N*(A+(1-T+C)*A*A*A/6                          #USGS Professional Paper 1395 (1987), p.61 (8-9)
                      + (5-18*T+T*T+72*C-58*eccprime_sq)*A*A*A*A*A/120)
                      + 500000.0)
  #Calculate Northing Coordinates
  UTMNorthing = (k0*(M+N*tan(lat_rad)*                           #USGS Professional Paper 1395 (1987), p.61 (8-10)
                      (A*A/2+(5-T+9*C+4*C*C)*A*A*A*A/24
                      + (61-58*T+T*T+600*C
                      - 330*eccprime_sq)*A*A*A*A*A*A/720)))
      
  #Adjust 1e7 meters if the coordinates are in the Southern hemisphere
  if (lat[i] < 0) {
    UTMNorthing = UTMNorthing + 1e7                              #USGS Professional Paper 1395 (1987), p.58, ¶2 
  } else {
    UTMNorthing = UTMNorthing
  }

  #Formulas retrieved from the U.S. Geological Survey (USGS) Professional Paper 1395 (1987) (p.61)
  #https://pubs.usgs.gov/pp/1395/report.pdf
  #Accessed on: January 10th, 2023
  
  # Assign Latitudinal Letter Designator ---------------------------------------
  UTMBand<-NULL                                                  #Empty vector shell
  
  for(i in 1:length(lat)){ 
    if(84 >= lat[i] && lat[i] >= 72){                            #X band
      UTMBand[i] <- "X"
    } else if(72 > lat[i] && lat[i] >= 64){                      #W band
      UTMBand[i] <- "W"
    } else if(64 > lat[i] && lat[i] >= 56){                      #V band
      UTMBand[i] <- "V"
    } else if(56 > lat[i] && lat[i] >= 48){                      #U band
      UTMBand[i] <- "U"
    } else if(48 > lat[i] && lat[i] >= 40){                      #T band
      UTMBand[i] <- "T"
    } else if(40 > lat[i] && lat[i] >= 32){                      #S band
      UTMBand[i] <- "S"
    } else if(32 > lat[i] && lat[i] >= 24){                      #R band
      UTMBand[i] <- "R"
    } else if(24 > lat[i] && lat[i] >= 16){                      #Q band
      UTMBand[i] <- "Q"
    } else if(16 > lat[i] && lat[i] >= 8){                       #P band
      UTMBand[i] <- "P"
    } else if(8 > lat[i] && lat[i] >= 0){                        #N band
      UTMBand[i] <- "N" 
    } else if(0 > lat[i] && lat[i] >= -8){                       #M band
      UTMBand[i] <- "M"
    } else if(-8 > lat[i] && lat[i] >= -16){                     #L band
      UTMBand[i] <- "L"
    } else if(-16 > lat[i] && lat[i] >= -24){                    #K band
      UTMBand[i] <- "K"
    } else if(-24 > lat[i] && lat[i] >= -32){                    #J band
      UTMBand[i] <- "J"
    } else if(-32 > lat[i] && lat[i] >= -40){                    #H band
      UTMBand[i] <- "H"
    } else if(-40 > lat[i] && lat[i] >= -48){                    #G band
      UTMBand[i] <- "G"
    } else if(-48 > lat[i] && lat[i] >= -56){                    #F band
      UTMBand[i] <- "F"
    } else if(-56 > lat[i] && lat[i] >= -64){                    #E band
      UTMBand[i] <- "E"
    } else if(-64 > lat[i] && lat[i] >= -72){                    #D band
      UTMBand[i] <- "D"
    } else if(-72 > lat[i] && lat[i] >= -80){                    #C band
      UTMBand[i] <- "C"
    } else{
      UTMBand[i] <- "Z"                                          #Z band
    }
  }
```
Each function begins with defining a set of ellipsoids ("datum") used for the three-dimensional projection of the coordinate system. Depending on the choice of datum, a set of ellipsoidal parameters are calculated. Those can be thought of as the parameters thats et forth our closest approximation to the shape of the earth (The choice of different ellipsoids is motivated by people fundamentally disagreeing on this set of parameters). 

The `_LLtoUTM.R` file then proceeds by
* Calculating the UTM Zone (while accounting for exceptions in the grid such as Svalbard)
* Obtaining the Northing and Easting coordinates through a sequence of geodetic formulas
* Assigning a UTM Band

The result is the conversion of decimal coordinates to UTM grid: 

![Chart1](https://earth-info.nga.mil/img/utm_fig1.jpg)

For suggestions, comments, or concerns, please do not hesitate to reach out to me at posmikdc[at]uchicago[dot]edu.
