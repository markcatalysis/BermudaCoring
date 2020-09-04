###################################################################################################################################
##’ @title Scripts to analyze coral calcification rate data from Bermuda coral cores with respect to the North Atlantic Oscillation
##’
##’ @author Travis A Courtney
##’ @contact traviscourtney@gmail.com
##’ @date 2020-09–04
###################################################################################################################################

#Load required libraries####################
library(readr)
library(dplyr)
library(nlme)
library(reshape2)
library(emmeans)
library(ggplot2)
library(patchwork)
library(dplR)
library(DescTools)
library(car)
library(stringr)

#Compile core data from scan 1####################

#Disable scientific notation
options(scipen=999)

#Establish importCoralXDS function for importing data
importCoralXDS = function(x){
  read_csv(x, col_names = FALSE, locale = locale(), skip = 9)
}

#navigate to core standards for scan 1
setwd("Scan_1/Coral_Standards")

#import scan 1 density standards
Coral_1_4_5_std1B <- importCoralXDS('Coral_1_4_5_std1B.dat')
Coral_1_4_4_std1C <- importCoralXDS('Coral_1_4_4_std1C.dat')
Coral_1_4_3_std1A <- importCoralXDS('Coral_1_4_3_std1A.dat')
Coral_1_4_2_std1A <- importCoralXDS('Coral_1_4_2_std1A.dat')
Coral_1_3_4_std1B <- importCoralXDS('Coral_1_3_4_std1B.dat')
Coral_1_3_3_std1A <- importCoralXDS('Coral_1_3_3_std1A.dat')
Coral_1_3_2_std2A <- importCoralXDS('Coral_1_3_2_std2A.dat')
Coral_1_2_4_std1B <- importCoralXDS('Coral_1_2_4_BONE_std1B.dat')
Coral_1_2_3_std3A <- importCoralXDS('Coral_1_2_3_BONE_std3A.dat')
Coral_1_2_2_std1C <- importCoralXDS('Coral_1_2_2_BONE_std1C.dat')
Coral_1_1_4_std2A <- importCoralXDS('Coral_1_1_4_std2A.dat')
Coral_1_1_3_std1A <- importCoralXDS('Coral_1_1_3_std1A.dat')
Coral_1_1_2_std1D <- importCoralXDS('Coral_1_1_2_std1D.dat')

#pair scan1 density standards intensity data with density measurements
ssid1B_1 = c(mean(Coral_1_4_5_std1B$X2),1.509)
ssid1C_1 = c(mean(Coral_1_4_4_std1C$X2),1.457)
pstr1A_1 = c(mean(Coral_1_4_3_std1A$X2),1.204)
ssid1A_1 = c(mean(Coral_1_4_2_std1A$X2),1.561)
carb1B_1 = c(mean(Coral_1_3_4_std1B$X2),1.840)
carb1A_1 = c(mean(Coral_1_3_3_std1A$X2),2.316)
ssid2A_1 = c(mean(Coral_1_3_2_std2A$X2),1.613)
mcav1B_1 = c(mean(Coral_1_2_4_std1B$X2),1.651)
ssid3A_1 = c(mean(Coral_1_2_3_std3A$X2),1.756)
mcav1C_1 = c(mean(Coral_1_2_2_std1C$X2),1.362)
pstr2A_1 = c(mean(Coral_1_1_4_std2A$X2),1.345)
unk1A_1 = c(mean(Coral_1_1_3_std1A$X2),0.846)
ssid1D_1 = c(mean(Coral_1_1_2_std1D$X2),1.503)

#compile dataframe of scan1 core standard intensity and density
standards_scan_1 = data.frame(rbind(ssid1A_1,ssid1B_1,ssid1C_1,ssid1D_1,ssid2A_1,ssid3A_1,pstr1A_1,pstr2A_1,mcav1B_1,mcav1C_1,carb1A_1,carb1B_1,unk1A_1))
colnames(standards_scan_1)=c("intensity","density")

#construct scan1 standards curve relating CT intensity to measured density
scan1_standards_curve = lm(density~intensity,data=standards_scan_1)
summary(scan1_standards_curve)

#construct intensity range data for subsequent plotting of standard calibration curve
intensityrange = as.data.frame(seq(130, 220, by = 1))
colnames(intensityrange)=c("intensity")

#construct core standards curve for importing data
corestandards1 = function(x){
  colnames(x)=c("intensity")
  intensity = as.data.frame(x)
  coredensity = as.data.frame(predict(scan1_standards_curve, newdata = x))
  conf.int <- as.data.frame(predict(scan1_standards_curve, newdata = x, interval = "confidence", level = 0.95))
  return(cbind(intensity,conf.int))
}

#fit corestandards data to intensity range for plotting
modelfit1 = corestandards1(intensityrange)

#navigate to scan 1 coral samples data
setwd('..')
setwd('..')
setwd("Scan_1/Coral_Samples")

#construct scan1 corals data import
coretransect_scan1=function(a,b,c,d,e,f,g,h,i,j,k,l){
  #import the data for each transect
  t1 <- read_csv(a, skip = 7,col_types=cols())
  ext_high1=subset(t1, select=extension, luminance>=0)
  lum_high1=subset(t1, select=luminance, luminance>=0)
  ext_low1=subset(t1, select=extension_1, luminance_1>=0)
  lum_low1=subset(t1, select=luminance_1, luminance_1>=0)
  
  t2 <- read_csv(b, skip = 7,col_types=cols())
  ext_high2=subset(t2, select=extension, luminance>=0)
  lum_high2=subset(t2, select=luminance, luminance>=0)
  ext_low2=subset(t2, select=extension_1, luminance_1>=0)
  lum_low2=subset(t2, select=luminance_1, luminance_1>=0)
  
  t3 <- read_csv(c, skip = 7,col_types=cols())
  ext_high3=subset(t3, select=extension, luminance>=0)
  lum_high3=subset(t3, select=luminance, luminance>=0)
  ext_low3=subset(t3, select=extension_1, luminance_1>=0)
  lum_low3=subset(t3, select=luminance_1, luminance_1>=0)
  
  t4 <- read_csv(d, skip = 7,col_types=cols())
  ext_high4=subset(t4, select=extension, luminance>=0)
  lum_high4=subset(t4, select=luminance, luminance>=0)
  ext_low4=subset(t4, select=extension_1, luminance_1>=0)
  lum_low4=subset(t4, select=luminance_1, luminance_1>=0)
  
  t5 <- read_csv(e, skip = 7,col_types=cols())
  ext_high5=subset(t5, select=extension, luminance>=0)
  lum_high5=subset(t5, select=luminance, luminance>=0)
  ext_low5=subset(t5, select=extension_1, luminance_1>=0)
  lum_low5=subset(t5, select=luminance_1, luminance_1>=0)
  
  t6 <- read_csv(f, skip = 7,col_types=cols())
  ext_high6=subset(t6, select=extension, luminance>=0)
  lum_high6=subset(t6, select=luminance, luminance>=0)
  ext_low6=subset(t6, select=extension_1, luminance_1>=0)
  lum_low6=subset(t6, select=luminance_1, luminance_1>=0)
  
  t7 <- read_csv(g, skip = 7,col_types=cols())
  ext_high7=subset(t7, select=extension, luminance>=0)
  lum_high7=subset(t7, select=luminance, luminance>=0)
  ext_low7=subset(t7, select=extension_1, luminance_1>=0)
  lum_low7=subset(t7, select=luminance_1, luminance_1>=0)
  
  t8 <- read_csv(h, skip = 7,col_types=cols())
  ext_high8=subset(t8, select=extension, luminance>=0)
  lum_high8=subset(t8, select=luminance, luminance>=0)
  ext_low8=subset(t8, select=extension_1, luminance_1>=0)
  lum_low8=subset(t8, select=luminance_1, luminance_1>=0)
  
  t9 <- read_csv(i, skip = 7,col_types=cols())
  ext_high9=subset(t9, select=extension, luminance>=0)
  lum_high9=subset(t9, select=luminance, luminance>=0)
  ext_low9=subset(t9, select=extension_1, luminance_1>=0)
  lum_low9=subset(t9, select=luminance_1, luminance_1>=0)
  
  t10 <- read_csv(j, skip = 7,col_types=cols())
  ext_high10=subset(t10, select=extension, luminance>=0)
  lum_high10=subset(t10, select=luminance, luminance>=0)
  ext_low10=subset(t10, select=extension_1, luminance_1>=0)
  lum_low10=subset(t10, select=luminance_1, luminance_1>=0)
  
  t11 <- read_csv(k, skip = 7,col_types=cols())
  ext_high11=subset(t11, select=extension, luminance>=0)
  lum_high11=subset(t11, select=luminance, luminance>=0)
  ext_low11=subset(t11, select=extension_1, luminance_1>=0)
  lum_low11=subset(t11, select=luminance_1, luminance_1>=0)
  
  t12 <- read_csv(l, skip = 7,col_types=cols())
  ext_high12=subset(t12, select=extension, luminance>=0)
  lum_high12=subset(t12, select=luminance, luminance>=0)
  ext_low12=subset(t12, select=extension_1, luminance_1>=0)
  lum_low12=subset(t12, select=luminance_1, luminance_1>=0)
  
  #bind all transect high and low intensity band data
  high=cbind(rbind(ext_high1,ext_high2,ext_high3,ext_high4,ext_high5,ext_high6,ext_high7,ext_high8,ext_high9,ext_high10,ext_high11,ext_high12),rbind(lum_high1,lum_high2,lum_high3,lum_high4,lum_high5,lum_high6,lum_high7,lum_high8,lum_high9,lum_high10,lum_high11,lum_high12))
  high = tibble::rownames_to_column(high, "year")
  low=cbind(rbind(ext_low1,ext_low2,ext_low3,ext_low4,ext_low5,ext_low6,ext_low7,ext_low8,ext_low9,ext_low10,ext_low11,ext_low12),rbind(lum_low1,lum_low2,lum_low3,lum_low4,lum_low5,lum_low6,lum_low7,lum_low8,lum_low9,lum_low10,lum_low11,lum_low12))
  low = tibble::rownames_to_column(low, "year")
  
  #join ext and intensity data for each year and calculate annual values
  core_data = as.data.frame(sapply(full_join(high,low,by='year'), function(x) as.numeric(x)))
  colnames(core_data)=c('year','ext_high','lum_high','ext_low','lum_low')
  ext_annual = with(core_data,ext_high+ext_low) #annual extension rates from paired bands
  lum_annual = with(core_data,(ext_high*lum_high+ext_low*lum_low)/(ext_high+ext_low)) #annual luminance intensities from paired bands
  
  #calculate high, low, and annual densities from from intensity data
  dens_high = as.data.frame(corestandards1(as.data.frame(core_data$lum_high))$fit) #calculate high density bands from high luminance data
  colnames(dens_high)=c("dens_high")
  dens_low = as.data.frame(corestandards1(as.data.frame(core_data$lum_low))$fit) #calculate low density bands from low luminance data
  colnames(dens_low)=c("dens_low")
  dens_annual = as.data.frame(corestandards1(as.data.frame(lum_annual))$fit) #calculate annual density bands from annual luminance data
  colnames(dens_annual)=c("dens_annual")
  
  #compile core extension and density data frame
  core_data_density = data.frame(core_data$year,core_data$ext_high,dens_high,core_data$ext_low,dens_low,ext_annual,dens_annual)
  colnames(core_data_density) = c('year','ext_high','dens_high','ext_low','dens_low','ext_annual','dens_annual')
  
  #calculate high, low, and annual calcification rate data from extension and density data
  calc_high = with(core_data_density,ext_high*dens_high)
  calc_low = with(core_data_density,ext_low*dens_low)
  calc_annual = with(core_data_density,ext_annual*dens_annual)
  
  #compile data frame of core physiological measurements
  core_data_calc = data.frame(core_data$year,core_data$ext_high,dens_high,calc_high,core_data$ext_low,dens_low,calc_low,ext_annual,dens_annual,calc_annual)
  colnames(core_data_calc) = c('year','ext_high','dens_high','calc_high','ext_low','dens_low','calc_low','ext_annual','dens_annual','calc_annual')
  
  #return core data frame
  return(core_data_calc)
}

OFRA4 = coretransect_scan1("OFRA4_Coral_1_1_6_t1.dat",
                           "OFRA4_Coral_1_1_6_t2.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

OFRA5 = coretransect_scan1("OFRA5_Coral_1_4_6_t1.dat",
                           "OFRA5_Coral_1_4_6_t1.5.dat",
                           "OFRA5_Coral_1_4_6_t2.dat",
                           "OFRA5_Coral_1_4_6_t3.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

OFRA6 = coretransect_scan1("OFRA6_Coral_1_4_7_t1.dat",
                           "OFRA6_Coral_1_4_7_t2.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

OFRA9 = coretransect_scan1("OFRA9_Coral_1_2_6_BONE_t1.dat",
                           "OFRA9_Coral_1_2_6_BONE_t2.dat",
                           "OFRA9_Coral_1_2_6_BONE_t3.dat",
                           "OFRA9_Coral_1_2_6_BONE_t4.dat",
                           "OFRA9_Coral_1_2_6_BONE_t5.dat",
                           "OFRA9_Coral_1_2_6_BONE_t6.dat",
                           "OFRA9_Coral_1_2_6_BONE_t7.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

OFRA10 = coretransect_scan1("OFRA10_Coral_1_3_6_t1.dat",
                            "OFRA10_Coral_1_3_6_t2.dat",
                            "OFRA10_Coral_1_3_6_t3.dat",
                            "OFRA10_Coral_1_3_6_t4.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

OFRA11 = coretransect_scan1("OFRA11_Coral_1_3_5_t1.dat",
                            "OFRA11_Coral_1_3_5_t2.dat",
                            "OFRA11_Coral_1_3_5_t3.dat",
                            "OFRA11_Coral_1_3_5_t4.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

OFRA14 = coretransect_scan1("OFRA14_Coral_1_1_5_t1.dat",
                            "OFRA14_Coral_1_1_5_t2.dat",
                            "OFRA14_Coral_1_1_5_t3.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

PSTR3 = coretransect_scan1("PSTR3_Coral_1_4_8_t1.dat",
                           "PSTR3_Coral_1_4_8_t2.dat",
                           "PSTR3_Coral_1_4_8_t3.dat",
                           "PSTR3_Coral_1_4_8_t4.dat",
                           "PSTR3_Coral_1_4_8_t5.dat",
                           "PSTR3_Coral_1_4_8_t6.dat",
                           "PSTR3_Coral_1_4_8_t7.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

PSTR4 = coretransect_scan1("PSTR4_Coral_1_1_7_t1.dat",
                           "PSTR4_Coral_1_1_7_t2.dat",
                           "PSTR4_Coral_1_1_7_t3.dat",
                           "PSTR4_Coral_1_1_7_t4.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

PSTR5 = coretransect_scan1("PSTR5_Coral_1_2_7_BONE_t1.dat",
                           "PSTR5_Coral_1_2_7_BONE_t2.dat",
                           "PSTR5_Coral_1_2_7_BONE_t3.dat",
                           "PSTR5_Coral_1_2_7_BONE_t4.dat",
                           "PSTR5_Coral_1_2_7_BONE_t5.dat",
                           "PSTR5_Coral_1_2_7_BONE_t6.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

PSTR6 = coretransect_scan1("PSTR6_Coral_1_3_t1.dat",
                           "PSTR6_Coral_1_3_t2.dat",
                           "PSTR6_Coral_1_3_t3.dat",
                           "PSTR6_Coral_1_3_t4.dat",
                           "NoData.dat", #"PSTR6_Coral_1_3_t5.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

PSTR7 = coretransect_scan1("PSTR7_Coral_1_2_5_BONE_t1.dat",
                           "PSTR7_Coral_1_2_5_BONE_t2.dat",
                           "PSTR7_Coral_1_2_5_BONE_t3.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

PSTR9 = coretransect_scan1("PSTR9_Coral_1_4_1_t1.dat",
                           "PSTR9_Coral_1_4_1_t2.dat",
                           "PSTR9_Coral_1_4_1_t3.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

PSTR10 = coretransect_scan1("PSTR10_Coral_1_1_1_t1.dat",
                            "PSTR10_Coral_1_1_1_t2.dat",
                            "PSTR10_Coral_1_1_1_t3.dat",
                            "PSTR10_Coral_1_1_1_t4.dat",
                            "PSTR10_Coral_1_1_1_t5.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

PSTR11 = coretransect_scan1("PSTR11_Coral_1_3_1_t1.dat",
                            "PSTR11_Coral_1_3_1_t2.dat",
                            "PSTR11_Coral_1_3_1_t3.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

PSTR12 = coretransect_scan1("PSTR12_Coral_1_2_1_BONE_t1.dat",
                            "PSTR12_Coral_1_2_1_BONE_t2.dat",
                            "PSTR12_Coral_1_2_1_BONE_t3.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

#navigate back to project working directory
setwd('..')
setwd('..')

#Compile core data from scan 2####################
setwd("Scan_2/Coral_Standards")

#import scan2 coral standard data
Coral_2_1_2_std1C <- importCoralXDS('Coral_2_1_2_std1C.dat')
Coral_2_1_3_std3A <- importCoralXDS('Coral_2_1_3_std3A.dat')
Coral_2_1_4_std1A <- importCoralXDS('Coral_2_1_4_std1A.dat')
Coral_2_1_5_std1D <- importCoralXDS('Coral_2_1_5_std1D.dat')
Coral_2_1_10_std1A <- importCoralXDS('Coral_2_1_10_std1A.dat')
Coral_2_1_11_std1A <- importCoralXDS('Coral_2_1_11_std1A.dat')
Coral_2_2_1_std1B <- importCoralXDS('Coral_2_2_1_std1B.dat')
Coral_2_2_7_std1C <- importCoralXDS('Coral_2_2_7_std1C.dat')
Coral_2_2_8_std1A <- importCoralXDS('Coral_2_2_8_std1A.dat')
Coral_2_2_9_std2A <- importCoralXDS('Coral_2_2_9_std2A.dat')
Coral_2_3_2_std1B <- importCoralXDS('Coral_2_3_2_std1B.dat')
Coral_2_3_3_std1B <- importCoralXDS('Coral_2_3_3_std1B.dat')
Coral_2_3_4_std2A <- importCoralXDS('Coral_2_3_4_std2A.dat')

#compile scan2 CT intensity data with measured densities
ssid1B_2 = c(mean(Coral_2_3_2_std1B$X2),1.509)
ssid1C_2 = c(mean(Coral_2_2_7_std1C$X2),1.457)
pstr1A_2 = c(mean(Coral_2_1_4_std1A$X2),1.204)
ssid1A_2 = c(mean(Coral_2_2_8_std1A$X2),1.561)
carb1B_2 = c(mean(Coral_2_2_1_std1B$X2),1.840)
carb1A_2 = c(mean(Coral_2_1_11_std1A$X2),2.316)
ssid2A_2 = c(mean(Coral_2_2_9_std2A$X2),1.613)
mcav1B_2 = c(mean(Coral_2_3_3_std1B$X2),1.651)
ssid3A_2 = c(mean(Coral_2_1_3_std3A$X2),1.756)
mcav1C_2 = c(mean(Coral_2_1_2_std1C$X2),1.362)
pstr2A_2 = c(mean(Coral_2_3_4_std2A$X2),1.345)
unk1A_2 = c(mean(Coral_2_1_10_std1A$X2),0.846)
ssid1D_2 = c(mean(Coral_2_1_5_std1D$X2),1.503)

#compile scan 2 core standards data frame
standards_scan_2 = data.frame(rbind(ssid1A_2,ssid1B_2,ssid1C_2,ssid1D_2,ssid2A_2,ssid3A_2,pstr1A_2,pstr2A_2,mcav1B_2,mcav1C_2,carb1A_2,carb1B_2,unk1A_2))
colnames(standards_scan_2)=c("intensity","density")

#evaluate standard calibration curve between measured density and CT intensity
scan2_standards_curve = lm(density~intensity,data=standards_scan_2)
summary(scan2_standards_curve)

#compile range to fit scan2 standards curve for plots
intensityrange = as.data.frame(seq(130, 220, by = 1))
colnames(intensityrange)=c("intensity")

#construct scan2 standards curve function for data imports of scan2 cores
corestandards2 = function(x){
  colnames(x)=c("intensity")
  intensity = as.data.frame(x)
  coredensity = as.data.frame(predict(scan2_standards_curve, newdata = x))
  conf.int <- as.data.frame(predict(scan2_standards_curve, newdata = x, interval = "confidence", level = 0.95))
  return(cbind(intensity,conf.int))
}

#fit predicted line for scan2 core standards
modelfit2 = corestandards2(intensityrange)

#navigate to scan 2 core data
setwd('..')
setwd('..')
setwd("Scan_2/Coral_Samples")

#construct scan2 core import function
coretransect_scan2=function(a,b,c,d,e,f,g,h,i,j,k,l){
  #imports data for each transect
  t1 <- read_csv(a, skip = 7,col_types=cols())
  ext_high1=subset(t1, select=extension, luminance>=0)
  lum_high1=subset(t1, select=luminance, luminance>=0)
  ext_low1=subset(t1, select=extension_1, luminance_1>=0)
  lum_low1=subset(t1, select=luminance_1, luminance_1>=0)
  
  t2 <- read_csv(b, skip = 7,col_types=cols())
  ext_high2=subset(t2, select=extension, luminance>=0)
  lum_high2=subset(t2, select=luminance, luminance>=0)
  ext_low2=subset(t2, select=extension_1, luminance_1>=0)
  lum_low2=subset(t2, select=luminance_1, luminance_1>=0)
  
  t3 <- read_csv(c, skip = 7,col_types=cols())
  ext_high3=subset(t3, select=extension, luminance>=0)
  lum_high3=subset(t3, select=luminance, luminance>=0)
  ext_low3=subset(t3, select=extension_1, luminance_1>=0)
  lum_low3=subset(t3, select=luminance_1, luminance_1>=0)
  
  t4 <- read_csv(d, skip = 7,col_types=cols())
  ext_high4=subset(t4, select=extension, luminance>=0)
  lum_high4=subset(t4, select=luminance, luminance>=0)
  ext_low4=subset(t4, select=extension_1, luminance_1>=0)
  lum_low4=subset(t4, select=luminance_1, luminance_1>=0)
  
  t5 <- read_csv(e, skip = 7,col_types=cols())
  ext_high5=subset(t5, select=extension, luminance>=0)
  lum_high5=subset(t5, select=luminance, luminance>=0)
  ext_low5=subset(t5, select=extension_1, luminance_1>=0)
  lum_low5=subset(t5, select=luminance_1, luminance_1>=0)
  
  t6 <- read_csv(f, skip = 7,col_types=cols())
  ext_high6=subset(t6, select=extension, luminance>=0)
  lum_high6=subset(t6, select=luminance, luminance>=0)
  ext_low6=subset(t6, select=extension_1, luminance_1>=0)
  lum_low6=subset(t6, select=luminance_1, luminance_1>=0)
  
  t7 <- read_csv(g, skip = 7,col_types=cols())
  ext_high7=subset(t7, select=extension, luminance>=0)
  lum_high7=subset(t7, select=luminance, luminance>=0)
  ext_low7=subset(t7, select=extension_1, luminance_1>=0)
  lum_low7=subset(t7, select=luminance_1, luminance_1>=0)
  
  t8 <- read_csv(h, skip = 7,col_types=cols())
  ext_high8=subset(t8, select=extension, luminance>=0)
  lum_high8=subset(t8, select=luminance, luminance>=0)
  ext_low8=subset(t8, select=extension_1, luminance_1>=0)
  lum_low8=subset(t8, select=luminance_1, luminance_1>=0)
  
  t9 <- read_csv(i, skip = 7,col_types=cols())
  ext_high9=subset(t9, select=extension, luminance>=0)
  lum_high9=subset(t9, select=luminance, luminance>=0)
  ext_low9=subset(t9, select=extension_1, luminance_1>=0)
  lum_low9=subset(t9, select=luminance_1, luminance_1>=0)
  
  t10 <- read_csv(j, skip = 7,col_types=cols())
  ext_high10=subset(t10, select=extension, luminance>=0)
  lum_high10=subset(t10, select=luminance, luminance>=0)
  ext_low10=subset(t10, select=extension_1, luminance_1>=0)
  lum_low10=subset(t10, select=luminance_1, luminance_1>=0)
  
  t11 <- read_csv(k, skip = 7,col_types=cols())
  ext_high11=subset(t11, select=extension, luminance>=0)
  lum_high11=subset(t11, select=luminance, luminance>=0)
  ext_low11=subset(t11, select=extension_1, luminance_1>=0)
  lum_low11=subset(t11, select=luminance_1, luminance_1>=0)
  
  t12 <- read_csv(l, skip = 7,col_types=cols())
  ext_high12=subset(t12, select=extension, luminance>=0)
  lum_high12=subset(t12, select=luminance, luminance>=0)
  ext_low12=subset(t12, select=extension_1, luminance_1>=0)
  lum_low12=subset(t12, select=luminance_1, luminance_1>=0)
  
  #compile all transect data into dataframe
  high=cbind(rbind(ext_high1,ext_high2,ext_high3,ext_high4,ext_high5,ext_high6,ext_high7,ext_high8,ext_high9,ext_high10,ext_high11,ext_high12),rbind(lum_high1,lum_high2,lum_high3,lum_high4,lum_high5,lum_high6,lum_high7,lum_high8,lum_high9,lum_high10,lum_high11,lum_high12))
  high = tibble::rownames_to_column(high, "year")
  low=cbind(rbind(ext_low1,ext_low2,ext_low3,ext_low4,ext_low5,ext_low6,ext_low7,ext_low8,ext_low9,ext_low10,ext_low11,ext_low12),rbind(lum_low1,lum_low2,lum_low3,lum_low4,lum_low5,lum_low6,lum_low7,lum_low8,lum_low9,lum_low10,lum_low11,lum_low12))
  low = tibble::rownames_to_column(low, "year")
  
  #join transect to construct annual extension and luminance
  core_data = as.data.frame(sapply(full_join(high,low,by='year'), function(x) as.numeric(x)))
  colnames(core_data)=c('year','ext_high','lum_high','ext_low','lum_low')
  ext_annual = with(core_data,ext_high+ext_low) #construct annual extension rates
  lum_annual = with(core_data,(ext_high*lum_high+ext_low*lum_low)/(ext_high+ext_low)) #construct annual luminance values
  
  #compile densities based on corestandards2 relationship
  dens_high = as.data.frame(corestandards2(as.data.frame(core_data$lum_high))$fit) #calculates high density bands from luminance
  colnames(dens_high)=c("dens_high")
  dens_low = as.data.frame(corestandards2(as.data.frame(core_data$lum_low))$fit) #calculates low density bands from luminance
  colnames(dens_low)=c("dens_low")
  dens_annual = as.data.frame(corestandards2(as.data.frame(lum_annual))$fit) #calculates annual density bands from luminance
  colnames(dens_annual)=c("dens_annual")
  
  #compile core extension and density data
  core_data_density = data.frame(core_data$year,core_data$ext_high,dens_high,core_data$ext_low,dens_low,ext_annual,dens_annual)
  colnames(core_data_density) = c('year','ext_high','dens_high','ext_low','dens_low','ext_annual','dens_annual')
  
  #calculate calcification rates based on extension and density data
  calc_high = with(core_data_density,ext_high*dens_high)
  calc_low = with(core_data_density,ext_low*dens_low)
  calc_annual = with(core_data_density,ext_annual*dens_annual)
  
  #compile dataframe of all core physiological data
  core_data_calc = data.frame(core_data$year,core_data$ext_high,dens_high,calc_high,core_data$ext_low,dens_low,calc_low,ext_annual,dens_annual,calc_annual)
  colnames(core_data_calc) = c('year','ext_high','dens_high','calc_high','ext_low','dens_low','calc_low','ext_annual','dens_annual','calc_annual')
  
  #return core data frame
  return(core_data_calc)
}

PSTR1 = coretransect_scan2("Coral_2_1_1_PSTR1_t1.dat",
                           "Coral_2_1_1_PSTR1_t2.dat",
                           "NoData.dat", #"Coral_2_1_1_PSTR1_t3.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

DLAB10 = coretransect_scan2("Coral_2_1_6_DLAB10_t1.dat",
                            "Coral_2_1_6_DLAB10_t2.dat",
                            "Coral_2_1_6_DLAB10_t3.dat",
                            "NoData.dat", #"Coral_2_1_6_DLAB10_t4.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                            "NoData.dat", #"Coral_2_1_6_DLAB10_t5.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

OFRA3 = coretransect_scan2("Coral_2_1_7_OFRA3_t1.dat",
                           "Coral_2_1_7_OFRA3_t2.dat",
                           "NoData.dat", # "Coral_2_1_7_OFRA3_t3.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

PSTR8 = coretransect_scan2("Coral_2_1_8_PSTR8_1_t1.dat",
                           "Coral_2_1_8_PSTR8_1_t2.dat",
                           "Coral_2_1_8_PSTR8_1_t3.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

OFRA1 = coretransect_scan2("Coral_2_2_2_OFRA1_t1.dat",
                           "Coral_2_2_2_OFRA1_t2.dat",
                           "Coral_2_2_2_OFRA1_t3.dat",
                           "Coral_2_2_2_OFRA1_t4.dat",
                           "Coral_2_2_2_OFRA1_t5.dat",
                           "Coral_2_2_2_OFRA1_t6.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

DLAB1 = coretransect_scan2("Coral_2_2_3_DLAB1_t1.dat",
                           "Coral_2_2_3_DLAB1_t2.dat",
                           "Coral_2_2_3_DLAB1_t3.dat",
                           "NoData.dat", #"Coral_2_2_3_DLAB1_t4.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

OFRA13 = coretransect_scan2("Coral_2_2_4_OFRA13_t1.dat",
                            "Coral_2_2_4_OFRA13_t2.dat",
                            "Coral_2_2_4_OFRA13_t3.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

PSTR2 = coretransect_scan2("Coral_2_2_5_PSTR2_1_t1.dat",
                           "Coral_2_2_5_PSTR2_1_t2.dat",
                           "Coral_2_2_5_PSTR2_1_t3.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

OFRA2 = coretransect_scan2("Coral_2_2_6_OFRA2_1_t1.dat",
                           "Coral_2_2_6_OFRA2_1_t2.dat",
                           "Coral_2_2_6_OFRA2_1_t3.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

OFRA12 = coretransect_scan2("Coral_2_2_10_OFRA12_t1.dat",
                            "Coral_2_2_10_OFRA12_t2.dat",
                            "Coral_2_2_10_OFRA12_t3.dat",
                            "Coral_2_2_10_OFRA12_t4.dat",
                            "Coral_2_2_10_OFRA12_t5.dat",
                            "Coral_2_2_10_OFRA12_t6.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

OFRA8 = coretransect_scan2("Coral_2_2_11_OFRA8_t1.dat",
                           "Coral_2_2_11_OFRA8_t2.dat",
                           "Coral_2_2_11_OFRA8_t3.dat",
                           "Coral_2_2_11_OFRA8_t4.dat",
                           "Coral_2_2_11_OFRA8_t5.dat",
                           "Coral_2_2_11_OFRA8_t6.dat",
                           "Coral_2_2_11_OFRA8_t7.dat",
                           "Coral_2_2_11_OFRA8_t8.dat",
                           "Coral_2_2_11_OFRA8_t9.dat",
                           "Coral_2_2_11_OFRA8_t10.dat",
                           "NoData.dat",
                           "NoData.dat")

OFRA7 = coretransect_scan2("Coral_2_3_1_OFRA7_t1.dat",
                           "Coral_2_3_1_OFRA7_t2.dat",
                           "Coral_2_3_1_OFRA7_t3.dat",
                           "Coral_2_3_1_OFRA7_t4.dat",
                           "Coral_2_3_1_OFRA7_t5.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

OFRA15 = coretransect_scan2("Coral_2_3_5_OFRA15_img2_t1.dat",
                            "Coral_2_3_5_OFRA15_img2_t2.dat",
                            "Coral_2_3_5_OFRA15_img2_t3.dat",
                            "Coral_2_3_5_OFRA15_img2_t4.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

DLAB6 = coretransect_scan2("Coral_2_3_6_DLAB6_1_t1.dat",
                           "Coral_2_3_6_DLAB6_1_t2.dat",
                           "Coral_2_3_6_DLAB6_1_t3.dat",
                           "Coral_2_3_6_DLAB6_1_t4.dat",
                           "Coral_2_3_6_DLAB6_1_t5.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")
setwd('..')
setwd('..')
#Compile core data from scan 3####################
setwd("Scan_3/Coral_Standards")

#import standards for Scan3
Coral_2s2_1_2_std1C <- importCoralXDS('Coral_2s2_1_2_std1C.dat')
Coral_2s2_1_3_std3A <- importCoralXDS('Coral_2s2_1_3_std3A.dat')
Coral_2s2_1_4_std1A <- importCoralXDS('Coral_2s2_1_4_std1A.dat')
Coral_2s2_1_5_std1B <- importCoralXDS('Coral_2s2_1_5_std1B.dat')
Coral_2s2_1_12_std1A <- importCoralXDS('Coral_2s2_1_12_std1A.dat')
Coral_2s2_1_11_std1A <- importCoralXDS('Coral_2s2_1_11_std1A.dat')
Coral_2s2_1_13_std2A_2 <- importCoralXDS('Coral_2s2_1_13_std2A_2.dat')
Coral_2s2_2_4_std2A <- importCoralXDS('Coral_2s2_2_4_std2A.dat')
Coral_2s2_2_5_std1A <- importCoralXDS('Coral_2s2_2_5_std1A.dat')
Coral_2s2_2_6_std1B <- importCoralXDS('Coral_2s2_2_6_std1B.dat')
Coral_2s2_3_2_std1B <- importCoralXDS('Coral_2s2_3_2_std1B.dat')
Coral_2s2_3_4_std1D <- importCoralXDS('Coral_2s2_3_4_std1D.dat')
Coral_2s2_4_2_aka_2s2_3_3_std1C <- importCoralXDS('Coral_2s2_4_2_aka_2s2_3_3_std1C.dat')

#pair Scan3 standards with the respective measured densities
ssid1B_3 = c(mean(Coral_2s2_3_2_std1B$X2),1.509)
ssid1C_3 = c(mean(Coral_2s2_4_2_aka_2s2_3_3_std1C$X2),1.457)
pstr1A_3 = c(mean(Coral_2s2_1_4_std1A$X2),1.204)
ssid1A_3 = c(mean(Coral_2s2_2_5_std1A$X2),1.561)
carb1B_3 = c(mean(Coral_2s2_2_6_std1B$X2),1.840)
carb1A_3 = c(mean(Coral_2s2_1_12_std1A$X2),2.316)
ssid2A_3 = c(mean(Coral_2s2_1_13_std2A_2$X2),1.613)
mcav1B_3 = c(mean(Coral_2s2_1_5_std1B$X2),1.651)
ssid3A_3 = c(mean(Coral_2s2_1_3_std3A$X2),1.756)
mcav1C_3 = c(mean(Coral_2s2_1_2_std1C$X2),1.362)
pstr2A_3 = c(mean(Coral_2s2_2_4_std2A$X2),1.345)
unk1A_3 = c(mean(Coral_2s2_1_11_std1A$X2),0.846)
ssid1D_3 = c(mean(Coral_2s2_3_4_std1D$X2),1.503)

#compile density standards data frame
standards_scan_3 = data.frame(rbind(ssid1A_3,ssid1B_3,ssid1C_3,ssid1D_3,ssid2A_3,ssid3A_3,pstr1A_3,pstr2A_3,mcav1B_3,mcav1C_3,carb1A_3,carb1B_3,unk1A_3))
colnames(standards_scan_3)=c("intensity","density")

#linear model of core standard density as a function of CT intensity
scan3_standards_curve = lm(density~intensity,data=standards_scan_3)
summary(scan3_standards_curve)

#construct intensity range to predict scan3 calibration curve
intensityrange = as.data.frame(seq(130, 220, by = 1))
colnames(intensityrange)=c("intensity")

#apply scan3 density curve to a function for all cores from scan 3
corestandards3 = function(x){
  colnames(x)=c("intensity")
  intensity = as.data.frame(x)
  coredensity = as.data.frame(predict(scan3_standards_curve, newdata = x))
  conf.int <- as.data.frame(predict(scan3_standards_curve, newdata = x, interval = "confidence", level = 0.95))
  return(cbind(intensity,conf.int))
}

#fit core standards data to intensity range for subsequent plot
modelfit3 = corestandards3(intensityrange)

#go back to scan 3 coral samples
setwd('..')
setwd('..')
setwd("Scan_3/Coral_Samples")

#function to import coral sample data for scan 3
coretransect_scan3=function(a,b,c,d,e,f,g,h,i,j,k,l){
  #imports data for each transect
  t1 <- read_csv(a, skip = 7,col_types=cols()) 
  ext_high1=subset(t1, select=extension, luminance>=0)
  lum_high1=subset(t1, select=luminance, luminance>=0)
  ext_low1=subset(t1, select=extension_1, luminance_1>=0)
  lum_low1=subset(t1, select=luminance_1, luminance_1>=0)
  
  t2 <- read_csv(b, skip = 7,col_types=cols())
  ext_high2=subset(t2, select=extension, luminance>=0)
  lum_high2=subset(t2, select=luminance, luminance>=0)
  ext_low2=subset(t2, select=extension_1, luminance_1>=0)
  lum_low2=subset(t2, select=luminance_1, luminance_1>=0)
  
  t3 <- read_csv(c, skip = 7,col_types=cols())
  ext_high3=subset(t3, select=extension, luminance>=0)
  lum_high3=subset(t3, select=luminance, luminance>=0)
  ext_low3=subset(t3, select=extension_1, luminance_1>=0)
  lum_low3=subset(t3, select=luminance_1, luminance_1>=0)
  
  t4 <- read_csv(d, skip = 7,col_types=cols())
  ext_high4=subset(t4, select=extension, luminance>=0)
  lum_high4=subset(t4, select=luminance, luminance>=0)
  ext_low4=subset(t4, select=extension_1, luminance_1>=0)
  lum_low4=subset(t4, select=luminance_1, luminance_1>=0)
  
  t5 <- read_csv(e, skip = 7,col_types=cols())
  ext_high5=subset(t5, select=extension, luminance>=0)
  lum_high5=subset(t5, select=luminance, luminance>=0)
  ext_low5=subset(t5, select=extension_1, luminance_1>=0)
  lum_low5=subset(t5, select=luminance_1, luminance_1>=0)
  
  t6 <- read_csv(f, skip = 7,col_types=cols())
  ext_high6=subset(t6, select=extension, luminance>=0)
  lum_high6=subset(t6, select=luminance, luminance>=0)
  ext_low6=subset(t6, select=extension_1, luminance_1>=0)
  lum_low6=subset(t6, select=luminance_1, luminance_1>=0)
  
  t7 <- read_csv(g, skip = 7,col_types=cols())
  ext_high7=subset(t7, select=extension, luminance>=0)
  lum_high7=subset(t7, select=luminance, luminance>=0)
  ext_low7=subset(t7, select=extension_1, luminance_1>=0)
  lum_low7=subset(t7, select=luminance_1, luminance_1>=0)
  
  t8 <- read_csv(h, skip = 7,col_types=cols())
  ext_high8=subset(t8, select=extension, luminance>=0)
  lum_high8=subset(t8, select=luminance, luminance>=0)
  ext_low8=subset(t8, select=extension_1, luminance_1>=0)
  lum_low8=subset(t8, select=luminance_1, luminance_1>=0)
  
  t9 <- read_csv(i, skip = 7,col_types=cols())
  ext_high9=subset(t9, select=extension, luminance>=0)
  lum_high9=subset(t9, select=luminance, luminance>=0)
  ext_low9=subset(t9, select=extension_1, luminance_1>=0)
  lum_low9=subset(t9, select=luminance_1, luminance_1>=0)
  
  t10 <- read_csv(j, skip = 7,col_types=cols())
  ext_high10=subset(t10, select=extension, luminance>=0)
  lum_high10=subset(t10, select=luminance, luminance>=0)
  ext_low10=subset(t10, select=extension_1, luminance_1>=0)
  lum_low10=subset(t10, select=luminance_1, luminance_1>=0)
  
  t11 <- read_csv(k, skip = 7,col_types=cols())
  ext_high11=subset(t11, select=extension, luminance>=0)
  lum_high11=subset(t11, select=luminance, luminance>=0)
  ext_low11=subset(t11, select=extension_1, luminance_1>=0)
  lum_low11=subset(t11, select=luminance_1, luminance_1>=0)
  
  t12 <- read_csv(l, skip = 7,col_types=cols())
  ext_high12=subset(t12, select=extension, luminance>=0)
  lum_high12=subset(t12, select=luminance, luminance>=0)
  ext_low12=subset(t12, select=extension_1, luminance_1>=0)
  lum_low12=subset(t12, select=luminance_1, luminance_1>=0)
  
  #compile data from each transect
  high=cbind(rbind(ext_high1,ext_high2,ext_high3,ext_high4,ext_high5,ext_high6,ext_high7,ext_high8,ext_high9,ext_high10,ext_high11,ext_high12),rbind(lum_high1,lum_high2,lum_high3,lum_high4,lum_high5,lum_high6,lum_high7,lum_high8,lum_high9,lum_high10,lum_high11,lum_high12))
  high = tibble::rownames_to_column(high, "year")
  low=cbind(rbind(ext_low1,ext_low2,ext_low3,ext_low4,ext_low5,ext_low6,ext_low7,ext_low8,ext_low9,ext_low10,ext_low11,ext_low12),rbind(lum_low1,lum_low2,lum_low3,lum_low4,lum_low5,lum_low6,lum_low7,lum_low8,lum_low9,lum_low10,lum_low11,lum_low12))
  low = tibble::rownames_to_column(low, "year")
  
  #convert to core data to numbers by year
  core_data = as.data.frame(sapply(full_join(high,low,by='year'), function(x) as.numeric(x)))
  colnames(core_data)=c('year','ext_high','lum_high','ext_low','lum_low')
  ext_annual = with(core_data,ext_high+ext_low) #compute annual extension
  lum_annual = with(core_data,(ext_high*lum_high+ext_low*lum_low)/(ext_high+ext_low)) #compute annual luminance
  
  dens_high = as.data.frame(corestandards3(as.data.frame(core_data$lum_high))$fit) #convert luminance of high density bands to actual density
  colnames(dens_high)=c("dens_high")
  dens_low = as.data.frame(corestandards3(as.data.frame(core_data$lum_low))$fit) #convert luminance of low density bands to actual density
  colnames(dens_low)=c("dens_low")
  dens_annual = as.data.frame(corestandards3(as.data.frame(lum_annual))$fit) #convert luminance of annual density bands to actual density
  colnames(dens_annual)=c("dens_annual")
  
  #compile all core extension and density data
  core_data_density = data.frame(core_data$year,core_data$ext_high,dens_high,core_data$ext_low,dens_low,ext_annual,dens_annual)
  colnames(core_data_density) = c('year','ext_high','dens_high','ext_low','dens_low','ext_annual','dens_annual')
  
  #comput seasonal and annual calcification rates
  calc_high = with(core_data_density,ext_high*dens_high)
  calc_low = with(core_data_density,ext_low*dens_low)
  calc_annual = with(core_data_density,ext_annual*dens_annual)
  
  #compile chronlogy for the respective core
  core_data_calc = data.frame(core_data$year,core_data$ext_high,dens_high,calc_high,core_data$ext_low,dens_low,calc_low,ext_annual,dens_annual,calc_annual)
  colnames(core_data_calc) = c('year','ext_high','dens_high','calc_high','ext_low','dens_low','calc_low','ext_annual','dens_annual','calc_annual')
  
  #output coral chronology
  return(core_data_calc)
}

DLAB9 = coretransect_scan3("Coral_2s2_1_1_DLAB9_t1.dat",
                           "Coral_2s2_1_1_DLAB9_t2.dat",
                           "Coral_2s2_1_1_DLAB9_t3.dat",
                           "Coral_2s2_1_1_DLAB9_t4.dat",
                           "Coral_2s2_1_1_DLAB9_t5.dat",
                           "Coral_2s2_1_1_DLAB9_t6.dat",
                           "Coral_2s2_1_1_DLAB9_t7.dat",
                           "Coral_2s2_1_1_DLAB9_t8.dat",
                           "Coral_2s2_1_1_DLAB9_t9.dat",
                           "Coral_2s2_1_1_DLAB9_t10.dat",
                           "Coral_2s2_1_1_DLAB9_t11.dat",
                           "Coral_2s2_1_1_DLAB9_t12.dat")

DLAB5 = coretransect_scan3("Coral_2s2_1_6_DLAB5-1_t1.dat",
                           "Coral_2s2_1_6_DLAB5-1_t2.dat",
                           "NoData.dat",#"Coral_2s2_1_6_DLAB5-1_t3.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

DLAB12 = coretransect_scan3("Coral_2s2_1_8_DLAB12_t1.dat",
                            "Coral_2s2_1_8_DLAB12_t2.dat",
                            "Coral_2s2_1_8_DLAB12_t3.dat",
                            "Coral_2s2_1_8_DLAB12_t4.dat",
                            "NoData.dat", #"Coral_2s2_1_8_DLAB12_t5.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

DLAB4 = coretransect_scan3("Coral_2s2_1_9_DLAB4-1_t1.dat",
                           "Coral_2s2_1_9_DLAB4-1_t2.dat",
                           "Coral_2s2_1_9_DLAB4-1_t3.dat",
                           "NoData.dat",#"Coral_2s2_1_9_DLAB4-1_t4.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

DLAB3 = coretransect_scan3("Coral_2s2_1_14_DLAB3-1_t1.dat",
                           "Coral_2s2_2_1_DLAB3-2_t1.dat",
                           "Coral_2s2_2_1_DLAB3-2_t2.dat",
                           "Coral_2s2_2_1_DLAB3-2_t3.dat",
                           "Coral_2s2_2_1_DLAB3-2_t4.dat",
                           "Coral_2s2_2_1_DLAB3-2_t5.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

DLAB2 = coretransect_scan3("Coral_2s2_2_2_DLAB2-1_t1.dat",
                           "Coral_2s2_2_2_DLAB2-1_t2.dat",
                           "Coral_2s2_2_2_DLAB2-1_t3.dat",
                           "Coral_2s2_2_2_DLAB2-1_t4.dat",
                           "Coral_2s2_2_2_DLAB2-1_t5.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

DLAB15 = coretransect_scan3("Coral_2s2_2_3_DLAB15_t1.dat",
                            "Coral_2s2_2_3_DLAB15_t2.dat",
                            "Coral_2s2_2_3_DLAB15_t3.dat",
                            "Coral_2s2_2_3_DLAB15_t4.dat",
                            "Coral_2s2_2_3_DLAB15_t5.dat",
                            "Coral_2s2_2_3_DLAB15_t6.dat",
                            "NoData.dat",#"Coral_2s2_2_3_DLAB15_t7.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                            "NoData.dat",#"Coral_2s2_2_3_DLAB15_t8.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                            "NoData.dat",#"Coral_2s2_2_3_DLAB15_t9.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

DLAB14 = coretransect_scan3("Coral_2s2_2_7_DLAB14_t1.dat",
                            "Coral_2s2_2_7_DLAB14_t2.dat",
                            "NoData.dat",#"Coral_2s2_2_7_DLAB14_t3.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

DLAB11 = coretransect_scan3("Coral_2s2_2_8_DLAB11_t1.dat",
                            "Coral_2s2_2_8_DLAB11_t2.dat",
                            "Coral_2s2_2_8_DLAB11_t3.dat",
                            "Coral_2s2_2_8_DLAB11_t4.dat",
                            "Coral_2s2_2_8_DLAB11_t5.dat",
                            "Coral_2s2_2_8_DLAB11_t6.dat",
                            "NoData.dat",#"Coral_2s2_2_8_DLAB11_t7.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                            "NoData.dat",#"Coral_2s2_2_8_DLAB11_t8.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

DLAB7 = coretransect_scan3("Coral_2s2_3_1_DLAB7_t1.dat",
                           "Coral_2s2_3_1_DLAB7_t2.dat",
                           "Coral_2s2_3_1_DLAB7_t3.dat",
                           "Coral_2s2_3_1_DLAB7_t4.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat",
                           "NoData.dat")

DLAB13 = coretransect_scan3("Coral_2s2_3_5_DLAB13_t1.dat",
                            "Coral_2s2_3_5_DLAB13_t2.dat",
                            "Coral_2s2_3_5_DLAB13_t3.dat",
                            "Coral_2s2_3_5_DLAB13_t4.dat",
                            "Coral_2s2_3_5_DLAB13_t5.dat",
                            "Coral_2s2_3_5_DLAB13_t6.dat",
                            "NoData.dat", #"Coral_2s2_3_5_DLAB13_t7.dat" excluded based on cross-validation of extension and density data with alternative cross-section of CT-scan
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat",
                            "NoData.dat")

DLAB8 = coretransect_scan3("Coral_2s2_3_6_DLAB8-1_t1.dat",
                           "Coral_2s2_3_6_DLAB8-1_t2.dat",
                           "Coral_2s2_3_6_DLAB8-1_t3.dat",
                           "Coral_2s2_3_6_DLAB8-1_t4.dat",
                           "Coral_2s2_3_6_DLAB8-1_t5.dat",
                           "Coral_2s2_3_7_DLAB8-2_t1.dat",
                           "Coral_2s2_3_7_DLAB8-2_t2.dat",
                           "Coral_2s2_3_7_DLAB8-2_t3.dat",
                           "Coral_2s2_3_7_DLAB8-2_t4.dat",
                           "Coral_2s2_3_7_DLAB8-2_t5.dat",
                           "Coral_2s2_3_7_DLAB8-2_t6.dat",
                           "Coral_2s2_3_7_DLAB8-2_t7.dat")
#return back to project working directory
setwd('..')
setwd('..')

#Compile all coral core growth chronologies####################
#construct annual_growth function to streamline data into annual growth summaries
annual_growth=function(coreID){annual_calc = coreID$calc_annual
annual_ext = coreID$ext_annual
annual_dens = coreID$dens_annual
coral = as.data.frame(rep(deparse(substitute(coreID),100)))
year = as.data.frame(2016-seq(0,99,1))
colnames(year)=c("year")
length(annual_calc)=100
length(annual_ext)=100
length(annual_dens)=100
annual_growth = cbind(year,annual_ext,annual_dens,annual_calc,coral)
colnames(annual_growth) = c("year","ext","dens","calc","coral")
return(annual_growth)}

#compile site names to add to the respective core
hogreef = as.data.frame(rep("hogreef",100))
colnames(hogreef)="site"
threehill = as.data.frame(rep("threehill",100))
colnames(threehill)="site"
whalebone = as.data.frame(rep("whalebone",100))
colnames(whalebone)="site"
gurnet = as.data.frame(rep("gurnet",100))
colnames(gurnet)="site"
halfway = as.data.frame(rep("halfway",100))
colnames(halfway)="site"

#append site names to respective annual coral growth data
DLAB1a=cbind(annual_growth(DLAB1),threehill)
DLAB2a=cbind(annual_growth(DLAB2),threehill)
DLAB3a=cbind(annual_growth(DLAB3),threehill)
DLAB4a=cbind(annual_growth(DLAB4),whalebone)
DLAB5a=cbind(annual_growth(DLAB5),whalebone)
DLAB6a=cbind(annual_growth(DLAB6),whalebone)
DLAB7a=cbind(annual_growth(DLAB7),hogreef)
DLAB8a=cbind(annual_growth(DLAB8),hogreef)
DLAB9a=cbind(annual_growth(DLAB9),hogreef)
DLAB10a=cbind(annual_growth(DLAB10),halfway)
DLAB11a=cbind(annual_growth(DLAB11),halfway)
DLAB12a=cbind(annual_growth(DLAB12),halfway)
DLAB13a=cbind(annual_growth(DLAB13),gurnet)
DLAB14a=cbind(annual_growth(DLAB14),gurnet)
DLAB15a=cbind(annual_growth(DLAB15),gurnet)

OFRA1a=cbind(annual_growth(OFRA1),whalebone)
OFRA2a=cbind(annual_growth(OFRA2),whalebone)
OFRA3a=cbind(annual_growth(OFRA3),whalebone)
OFRA4a=cbind(annual_growth(OFRA4),hogreef)
OFRA5a=cbind(annual_growth(OFRA5),hogreef)
OFRA6a=cbind(annual_growth(OFRA6),hogreef)
OFRA7a=cbind(annual_growth(OFRA7),threehill)
OFRA8a=cbind(annual_growth(OFRA8),threehill)
OFRA9a=cbind(annual_growth(OFRA9),threehill)
OFRA10a=cbind(annual_growth(OFRA10),halfway)
OFRA11a=cbind(annual_growth(OFRA11),halfway)
OFRA12a=cbind(annual_growth(OFRA12),halfway)
OFRA13a=cbind(annual_growth(OFRA13),gurnet)
OFRA14a=cbind(annual_growth(OFRA14),gurnet)
OFRA15a=cbind(annual_growth(OFRA15),gurnet)

PSTR1a=cbind(annual_growth(PSTR1),whalebone)
PSTR2a=cbind(annual_growth(PSTR2),whalebone)
PSTR3a=cbind(annual_growth(PSTR3),hogreef)
PSTR4a=cbind(annual_growth(PSTR4),hogreef)
PSTR5a=cbind(annual_growth(PSTR5),hogreef)
PSTR6a=cbind(annual_growth(PSTR6),threehill)
PSTR7a=cbind(annual_growth(PSTR7),threehill)
PSTR8a=cbind(annual_growth(PSTR8),threehill)
PSTR9a=cbind(annual_growth(PSTR9),whalebone)
PSTR10a=cbind(annual_growth(PSTR10),gurnet)
PSTR11a=cbind(annual_growth(PSTR11),gurnet)
PSTR12a=cbind(annual_growth(PSTR12),gurnet)

#compile annual coral growth records into single dataframe
annual_growth=rbind(
  DLAB1a,DLAB2a,DLAB3a,DLAB4a,DLAB5a,DLAB6a,DLAB7a,DLAB8a,DLAB9a,DLAB10a,DLAB11a,DLAB12a,DLAB13a,DLAB14a,DLAB15a,
  OFRA1a,OFRA2a,OFRA3a,OFRA4a,OFRA5a,OFRA6a,OFRA7a,OFRA8a,OFRA9a,OFRA10a,OFRA11a,OFRA12a,OFRA13a,OFRA14a,OFRA15a,
  PSTR1a,PSTR2a,PSTR3a,PSTR4a,PSTR5a,PSTR6a,PSTR7a,PSTR8a,PSTR9a,PSTR10a,PSTR11a,PSTR12a)

#compile species column from coreID
annual_growth$species = substr(annual_growth$coral, start = 1, stop = 4)

#remove incomplete 2016 year to avoid downward bias associated with an incomplete year of growth (cores were collected summer 2016)
annual_growth=subset(annual_growth,year<2016)

#Calculate mean core chronologies####################
#generate mean ± 95% extension chronology for DLAB cores
DLAB_ext = dcast(subset(annual_growth,species=="DLAB"), year ~ coral, value.var="ext")
meanDLAB_ext = as.data.frame(cbind(DLAB_ext[,1],t(apply(DLAB_ext[,-1], 1, MeanCI,na.rm=TRUE,conf.level=0.95))))
colnames(meanDLAB_ext)=c("year","mean","lwr","upr")

#generate mean ± 95% extension chronology for PSTR cores
PSTR_ext = dcast(subset(annual_growth,species=="PSTR"), year ~ coral, value.var="ext")
meanPSTR_ext = as.data.frame(cbind(PSTR_ext[,1],t(apply(PSTR_ext[,-1], 1, MeanCI,na.rm=TRUE,conf.level=0.95))))
colnames(meanPSTR_ext)=c("year","mean","lwr","upr")

#generate mean ± 95% extension chronology for OFRA cores
OFRA_ext = dcast(subset(annual_growth,species=="OFRA"), year ~ coral, value.var="ext")
meanOFRA_ext = as.data.frame(cbind(OFRA_ext[,1],t(apply(OFRA_ext[,-1], 1, MeanCI,na.rm=TRUE,conf.level=0.95))))
colnames(meanOFRA_ext)=c("year","mean","lwr","upr")

#generate mean ± 95% density chronology for DLAB cores
DLAB_dens = dcast(subset(annual_growth,species=="DLAB"), year ~ coral, value.var="dens")
meanDLAB_dens = as.data.frame(cbind(DLAB_dens[,1],t(apply(DLAB_dens[,-1], 1, MeanCI,na.rm=TRUE,conf.level=0.95))))
colnames(meanDLAB_dens)=c("year","mean","lwr","upr")

#generate mean ± 95% density chronology for PSTR cores
PSTR_dens = dcast(subset(annual_growth,species=="PSTR"), year ~ coral, value.var="dens")
meanPSTR_dens = as.data.frame(cbind(PSTR_dens[,1],t(apply(PSTR_dens[,-1], 1, MeanCI,na.rm=TRUE,conf.level=0.95))))
colnames(meanPSTR_dens)=c("year","mean","lwr","upr")

#generate mean ± 95% density chronology for OFRA cores
OFRA_dens = dcast(subset(annual_growth,species=="OFRA"), year ~ coral, value.var="dens")
meanOFRA_dens = as.data.frame(cbind(OFRA_dens[,1],t(apply(OFRA_dens[,-1], 1, MeanCI,na.rm=TRUE,conf.level=0.95))))
colnames(meanOFRA_dens)=c("year","mean","lwr","upr")

#generate mean ± 95% calcification chronology for DLAB cores
DLAB_calc = dcast(subset(annual_growth,species=="DLAB"), year ~ coral, value.var="calc")
meanDLAB_calc = as.data.frame(cbind(DLAB_calc[,1],t(apply(DLAB_calc[,-1], 1, MeanCI,na.rm=TRUE,conf.level=0.95))))
colnames(meanDLAB_calc)=c("year","mean","lwr","upr")

#generate mean ± 95% calcification chronology for PSTR cores
PSTR_calc = dcast(subset(annual_growth,species=="PSTR"), year ~ coral, value.var="calc")
meanPSTR_calc = as.data.frame(cbind(PSTR_calc[,1],t(apply(PSTR_calc[,-1], 1, MeanCI,na.rm=TRUE,conf.level=0.95))))
colnames(meanPSTR_calc)=c("year","mean","lwr","upr")

#generate mean ± 95% calcification chronology for OFRA cores
OFRA_calc = dcast(subset(annual_growth,species=="OFRA"), year ~ coral, value.var="calc")
meanOFRA_calc = as.data.frame(cbind(OFRA_calc[,1],t(apply(OFRA_calc[,-1], 1, MeanCI,na.rm=TRUE,conf.level=0.95))))
colnames(meanOFRA_calc)=c("year","mean","lwr","upr")

#Construct NAO timeseries and merge with growth chronologies####################
#NAO data downloaded from https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/nao.shtml
NAO_data <- read_table2("NAO_data.txt", col_names = FALSE)
colnames(NAO_data)=c("year","month","NAO")
NAO_data$decy = NAO_data$year + (NAO_data$month-0.5)/12 #compute decimal year centered at the middle of each month
colnames(NAO_data)=c("year","month","NAO", "decy")

NAO_sub = subset(NAO_data, month < 4 | month == 12) #subset December, January, February, and March to construct NAO (DJFM) winter index
date_winter = NAO_sub$year + NAO_sub$month/12 #construct annual index for NAO (DJFM)
NAO_winter=as.data.frame(cbind(round(date_winter,0),NAO_sub$NAO)) #round to the nearest year to calculate NAO_winter from mean NAO for DJFM
colnames(NAO_winter)=c("year","NAO")

#compute means±CI of NAO winter index for each winter
NAO_winter_means=
  NAO_winter %>%
  group_by(year) %>%
  dplyr::summarize(NAO_winter = MeanCI(NAO,na.rm=TRUE,conf.level=0.95)[1], NAO_winter_lwr = MeanCI(NAO,na.rm=TRUE,conf.level=0.95)[2],NAO_winter_upr = MeanCI(NAO,na.rm=TRUE,conf.level=0.95)[3])

#merge NAO data with coral annual growth data for subsequent linear mixed effects modeling
growth_NAO=merge(annual_growth,NAO_winter_means,by='year')

#Evaluate mean calcification as a function of winter NAO####################
#construct simple null models and add NAO and random effects
c1 = lm(calc~1,data=growth_NAO,na.action=na.omit)
c2 = lme(calc~1,random=~1|coral,data=growth_NAO,method="ML",na.action=na.omit)
c3 = lm(calc~NAO_winter,data=growth_NAO,na.action=na.omit)
c4 = lme(calc~NAO_winter,random=~1|coral,data=growth_NAO,method="ML",na.action=na.omit)
c5 = lme(calc~NAO_winter,random=~NAO_winter|coral,data=growth_NAO,method="ML",na.action=na.omit)

AIC(c1,c2,c3,c4,c5) #AIC selects c5 as the optimal model moving forward

c6 = lme(calc~NAO_winter+species,random=~NAO_winter|coral,data=growth_NAO,method="ML",na.action=na.omit)
c7 = lme(calc~NAO_winter+species+NAO_winter:species,random=~NAO_winter|coral,data=growth_NAO,method="ML",na.action=na.omit)
c8 = lme(calc~NAO_winter+site,random=~NAO_winter|coral,data=growth_NAO,method="ML",na.action=na.omit)
c9 = lme(calc~NAO_winter+site+NAO_winter:site,random=~NAO_winter|coral,data=growth_NAO,method="ML",na.action=na.omit)
c10 = lme(calc~NAO_winter+species+site,random=~NAO_winter|coral,data=growth_NAO,method="ML",na.action=na.omit)
c11 = lme(calc~NAO_winter+species+NAO_winter:species+site,random=~NAO_winter|coral,data=growth_NAO,method="ML",na.action=na.omit)
c12 = lme(calc~NAO_winter+site+NAO_winter:site+species,random=~NAO_winter|coral,data=growth_NAO,method="ML",na.action=na.omit)
c13 = lme(calc~NAO_winter+species+site+NAO_winter:species+NAO_winter:site,random=~NAO_winter|coral,data=growth_NAO,method="ML",na.action=na.omit)

AIC(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13) #AIC selects c13 as the model with lowest AIC

summary(c13) #summary of model 13
Anova(c13,type=3) #all predictors are significant at p<0.05 except for the main effect of site

qqnorm(resid(c13)) #visual assessment of normality
plot(resid(c13)) #visual assessment of homoscedasticity

c13_species_int=emtrends(c13, pairwise ~ species, var="NAO_winter",adjust="tukey") #only DLAB is negative and is sig diff from OFRA and PSTR
c13_site_int=emtrends(c13, pairwise ~ site, var="NAO_winter",adjust="tukey") #all sites overlap with zero, Whalebone sig more positive than Hog Reef
c13_species=emmeans(c13, pairwise~species, adjust="tukey")
c13_site=emmeans(c13, pairwise~site, adjust="tukey")

#Compile mean core metadata####################
#load coring metadata file
metadata=read_csv("coringmetadata.csv")
colnames(metadata)=c("species","coral","site","depth","tissue")

#calculate mean extension rate data for 2010-2015 for each core
mean_core_ext = subset(annual_growth,year>2009) %>%
  group_by(coral) %>%
  summarize(ext = mean(ext))
#append mean core extension rates to core metadata
metadata_ext<-merge(x=metadata,y=mean_core_ext,by="coral")

#generate mean densities for each core using previous import functions for the entire cored material data
setwd("Scan_1/Coral_Densities")
OFRA4_dens=corestandards1(as.data.frame(mean(importCoralXDS('OFRA4_Coral_1_1_6.dat')$X2)))
OFRA5_dens=corestandards1(as.data.frame(mean(importCoralXDS('OFRA5_Coral_1_4_6.dat')$X2)))
OFRA6_dens=corestandards1(as.data.frame(mean(importCoralXDS('OFRA6_Coral_1_4_7.dat')$X2)))
OFRA9_dens=corestandards1(as.data.frame(mean(importCoralXDS('OFRA9_Coral_1_2_6_BONE.dat')$X2)))
OFRA10_dens=corestandards1(as.data.frame(mean(importCoralXDS('OFRA10_Coral_1_3_6.dat')$X2)))
OFRA11_dens=corestandards1(as.data.frame(mean(importCoralXDS('OFRA11_Coral_1_3_5.dat')$X2)))
OFRA14_dens=corestandards1(as.data.frame(mean(importCoralXDS('OFRA14_Coral_1_1_5.dat')$X2)))
PSTR3_dens=corestandards1(as.data.frame(mean(importCoralXDS('PSTR3_Coral_1_4_8.dat')$X2)))
PSTR4_dens=corestandards1(as.data.frame(mean(importCoralXDS('PSTR4_Coral_1_1_7.dat')$X2)))
PSTR5_dens=corestandards1(as.data.frame(mean(importCoralXDS('PSTR5_Coral_1_2_7_BONE.dat')$X2)))
PSTR6_dens=corestandards1(as.data.frame(mean(importCoralXDS('PSTR6_Coral_1_3_7.dat')$X2)))
PSTR7_dens=corestandards1(as.data.frame(mean(importCoralXDS('PSTR7_Coral_1_2_5_BONE.dat')$X2)))
PSTR9_dens=corestandards1(as.data.frame(mean(importCoralXDS('PSTR9_Coral_1_4_1.dat')$X2)))
PSTR10_dens=corestandards1(as.data.frame(mean(importCoralXDS('PSTR10_Coral_1_1_1.dat')$X2)))
PSTR11_dens=corestandards1(as.data.frame(mean(importCoralXDS('PSTR11_Coral_1_3_1.dat')$X2)))
PSTR12_dens=corestandards1(as.data.frame(mean(importCoralXDS('PSTR12_Coral_1_2_1_BONE.dat')$X2)))
setwd('..')
setwd('..')
setwd("Scan_2/Coral_Densities")
PSTR1_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_1_1_PSTR1.dat')$X2)))
DLAB10_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_1_6_DLAB10.dat')$X2)))
OFRA3_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_1_7_OFRA3.dat')$X2)))
PSTR8_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_1_8_PSTR8_1.dat')$X2)))
OFRA1_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_2_2_OFRA1.dat')$X2)))
DLAB1_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_2_3_DLAB1.dat')$X2)))
OFRA13_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_2_4_OFRA13.dat')$X2)))
PSTR2_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_2_5_PSTR2_1.dat')$X2)))
OFRA2_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_2_6_OFRA2_1.dat')$X2)))
OFRA12_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_2_10_OFRA12.dat')$X2)))
OFRA8_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_2_11_OFRA8.dat')$X2)))
OFRA7_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_3_1_OFRA7.dat')$X2)))
OFRA15_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_3_5_OFRA15_img2.dat')$X2)))
DLAB6_dens=corestandards2(as.data.frame(mean(importCoralXDS('Coral_2_3_6_DLAB6_1.dat')$X2)))
setwd('..')
setwd('..')
setwd("Scan_3/Coral_Densities")
DLAB9_dens=corestandards3(as.data.frame(mean(importCoralXDS('Coral_2s2_1_1_DLAB9.dat')$X2)))
DLAB5_dens=corestandards3(as.data.frame(mean(importCoralXDS('Coral_2s2_1_6_DLAB5-1.dat')$X2)))
DLAB12_dens=corestandards3(as.data.frame(mean(importCoralXDS('Coral_2s2_1_8_DLAB12.dat')$X2)))
DLAB4_dens=corestandards3(as.data.frame(mean(importCoralXDS('Coral_2s2_1_9_DLAB4-1.dat')$X2)))
DLAB3_dens=corestandards3(as.data.frame(mean(rbind(importCoralXDS('Coral_2s2_1_14_DLAB3-1.dat'),importCoralXDS('Coral_2s2_2_1_DLAB3-2.dat'))$X2)))
DLAB2_dens=corestandards3(as.data.frame(mean(importCoralXDS('Coral_2s2_2_2_DLAB2-1.dat')$X2)))
DLAB15_dens=corestandards3(as.data.frame(mean(importCoralXDS('Coral_2s2_2_3_DLAB15.dat')$X2)))
DLAB14_dens=corestandards3(as.data.frame(mean(importCoralXDS('Coral_2s2_2_7_DLAB14.dat')$X2)))
DLAB11_dens=corestandards3(as.data.frame(mean(importCoralXDS('Coral_2s2_2_8_DLAB11.dat')$X2)))
DLAB7_dens=corestandards3(as.data.frame(mean(importCoralXDS('Coral_2s2_3_1_DLAB7.dat')$X2)))
DLAB13_dens=corestandards3(as.data.frame(mean(importCoralXDS('Coral_2s2_3_5_DLAB13.dat')$X2)))
DLAB8_dens=corestandards3(as.data.frame(mean(rbind(importCoralXDS('Coral_2s2_3_6_DLAB8-1.dat'),importCoralXDS('Coral_2s2_3_7_DLAB8-2.dat'))$X2)))
setwd('..')
setwd('..')

#compile all mean core density data into dataframe
cores=rbind("DLAB1","DLAB2","DLAB3","DLAB4","DLAB5","DLAB6","DLAB7","DLAB8","DLAB9","DLAB10","DLAB11","DLAB12","DLAB13","DLAB14","DLAB15",
            "OFRA1","OFRA2","OFRA3","OFRA4","OFRA5","OFRA6","OFRA7","OFRA8","OFRA9","OFRA10","OFRA11","OFRA12","OFRA13","OFRA14","OFRA15",
            "PSTR1","PSTR2","PSTR3","PSTR4","PSTR5","PSTR6","PSTR7","PSTR8","PSTR9","PSTR10","PSTR11","PSTR12")

core_dens=rbind(DLAB1_dens,DLAB2_dens,DLAB3_dens,DLAB4_dens,DLAB5_dens,DLAB6_dens,DLAB7_dens,DLAB8_dens,DLAB9_dens,DLAB10_dens,DLAB11_dens,DLAB12_dens,DLAB13_dens,DLAB14_dens,DLAB15_dens,
                OFRA1_dens,OFRA2_dens,OFRA3_dens,OFRA4_dens,OFRA5_dens,OFRA6_dens,OFRA7_dens,OFRA8_dens,OFRA9_dens,OFRA10_dens,OFRA11_dens,OFRA12_dens,OFRA13_dens,OFRA14_dens,OFRA15_dens,
                PSTR1_dens,PSTR2_dens,PSTR3_dens,PSTR4_dens,PSTR5_dens,PSTR6_dens,PSTR7_dens,PSTR8_dens,PSTR9_dens,PSTR10_dens,PSTR11_dens,PSTR12_dens)

mean_core_densities = as.data.frame(cbind(cores,as.numeric(as.character(core_dens$fit))))
colnames(mean_core_densities)=c("coral","dens")

#Calculate uncertainty in densities from standard density vs intensity curves
mean_dens_uncertainty = mean(core_dens$fit-core_dens$lwr)
percent_dens_uncertainty = (mean_dens_uncertainty/mean(core_dens$fit))*100

#append mean core densities to coring metadata
metadata_growth<-merge(x=metadata_ext,y=mean_core_densities,by="coral")

#generate mean core calcification rates
metadata_growth$dens=as.numeric(as.character(metadata_growth$dens))
metadata_growth$calc = metadata_growth$ext*metadata_growth$dens

#summarize mean core physiological measurements by species
mean_core_metadata = metadata_growth %>%
  group_by(species) %>%
  summarize(tissue_mean = mean(tissue), tissue_sd = sd(tissue), ext_mean= mean(ext), ext_sd= sd(ext), dens_mean= mean(dens), dens_sd = sd(dens), calc_mean = mean(calc), calc_sd = sd(calc), n = length(calc))

#determine mean, max, and min presence of years in time series for each species
time_series_lengths = annual_growth[complete.cases(annual_growth),] %>%
  group_by(species) %>%
  summarize(mean_year = mean(year), max_year = min(year), min_year = max(year))

#compile all core meta data
mean_core_metadata = merge(mean_core_metadata,time_series_lengths,by="species")

#summarize mean core physiological measurements by species and sites
mean_core_metadata_sites=
  metadata_growth %>%
  group_by(species,site) %>%
  summarize(tissue_mean = mean(tissue), tissue_sd = sd(tissue), ext_mean= mean(ext), ext_sd= sd(ext), dens_mean= mean(dens), dens_sd = sd(dens), calc_mean = mean(calc), calc_sd = sd(calc), n = length(calc))

#Evaluate mean core metadata####################
#construct mean core physiology vs. species + site models
tissue_lm = lm(tissue~species+site,data=metadata_growth)
ext_lm = lm(ext~species+site,data=metadata_growth)
dens_lm = lm(dens~species+site,data=metadata_growth)
calc_lm = lm(calc~species+site,data=metadata_growth)

anova(tissue_lm) #species is significant, site is not at p=0.05 threshold
emmeans(tissue_lm, pairwise~species, adjust="tukey") #DLAB = PSTR tissue thickness > OFRA

anova(ext_lm) #species is significant, site is not at p=0.05 threshold
emmeans(ext_lm, pairwise~species, adjust="tukey") #DLAB = PSTR extension rates > OFRA

anova(dens_lm) #species and site are significant predictors of core densities
emmeans(dens_lm, pairwise~species, adjust="tukey") #DLAB density < PSTR < OFRA 
emmeans(dens_lm, pairwise~site, adjust="tukey") #halfway density > gurnet, threehill, whalebone

anova(calc_lm) #species is significant, but site is not at p=0.05 threshold
emmeans(calc_lm, pairwise~species, adjust="tukey") #DLAB = OFRA calc > PSTR

#Determine extension, density, and calcification rate pointer years####################
#Extension: construct data frames of the format required for pointer() function for each species
DLAB_ep = DLAB_ext[,-1]
rownames(DLAB_ep)=DLAB_ext[,1]

PSTR_ep = PSTR_ext[,-1]
rownames(PSTR_ep)=PSTR_ext[,1]

OFRA_ep = OFRA_ext[,-1]
rownames(OFRA_ep)=OFRA_ext[,1]

#Determine pointer years for all cores from each species
DLAB_epointer=pointer(DLAB_ep, rgv.thresh=10, nseries.thresh=75)
PSTR_epointer=pointer(PSTR_ep, rgv.thresh=10, nseries.thresh=75)
OFRA_epointer=pointer(OFRA_ep, rgv.thresh=10, nseries.thresh=75)

#Subset pointer years based on additional criteria that at least 1/3 of chronologies are present for each species
DLAB_ext_point=subset(DLAB_epointer, Nature != 0 & Nb.series > max(Nb.series)/3) #no detectable pointer years
PSTR_ext_point=subset(PSTR_epointer, Nature != 0 & Nb.series > max(Nb.series)/3) #no detectable pointer years
OFRA_ext_point=subset(OFRA_epointer, Nature != 0 & Nb.series > max(Nb.series)/3) #no detectable pointer years

#Density: construct data frames of the format required for pointer() function for each species
DLAB_dp = DLAB_dens[,-1]
rownames(DLAB_dp)=DLAB_dens[,1]

PSTR_dp = PSTR_dens[,-1]
rownames(PSTR_dp)=PSTR_dens[,1]

OFRA_dp = OFRA_dens[,-1]
rownames(OFRA_dp)=OFRA_dens[,1]

#Determine pointer years for all cores from each species
DLAB_dpointer=pointer(DLAB_dp, rgv.thresh=10, nseries.thresh=75)
PSTR_dpointer=pointer(PSTR_dp, rgv.thresh=10, nseries.thresh=75)
OFRA_dpointer=pointer(OFRA_dp, rgv.thresh=10, nseries.thresh=75)

#Subset pointer years based on additional criteria that at least 1/3 of chronologies are present for each species
DLAB_dens_point=subset(DLAB_dpointer, Nature != 0 & Nb.series > max(Nb.series)/3) #pointer years in 1988 and 1989
PSTR_dens_point=subset(PSTR_dpointer, Nature != 0 & Nb.series > max(Nb.series)/3) #no detectable pointer years
OFRA_dens_point=subset(OFRA_dpointer, Nature != 0 & Nb.series > max(Nb.series)/3) #no detectable pointer years

#Calcification: construct data frames of the format required for pointer() function for each species
DLAB_cp = DLAB_calc[,-1]
rownames(DLAB_cp)=DLAB_calc[,1]

PSTR_cp = PSTR_calc[,-1]
rownames(PSTR_cp)=PSTR_calc[,1]

OFRA_cp = OFRA_calc[,-1]
rownames(OFRA_cp)=OFRA_calc[,1]

#Determine pointer years for all cores from each species
DLAB_pointer=pointer(DLAB_cp, rgv.thresh=10, nseries.thresh=75)
PSTR_pointer=pointer(PSTR_cp, rgv.thresh=10, nseries.thresh=75)
OFRA_pointer=pointer(OFRA_cp, rgv.thresh=10, nseries.thresh=75)

#Subset pointer years based on additional criteria that at least 1/3 of chronologies are present for each species
DLAB_calc_point=subset(DLAB_pointer, Nature != 0 & Nb.series > max(Nb.series)/3) #pointer years in 1988 and 1989
PSTR_calc_point=subset(PSTR_pointer, Nature != 0 & Nb.series > max(Nb.series)/3) #no detectable pointer years
OFRA_calc_point=subset(OFRA_pointer, Nature != 0 & Nb.series > max(Nb.series)/3) #no detectable pointer years

#lower threshold to 0 to determine how many DLAB cores exhibited positive or negative calcification anomalies in 1988 and 1989
DLAB_calc_point_zero_1988 = pointer(DLAB_cp, rgv.thresh=0, nseries.thresh=75)[71,] #100% of cores had negative calcification anomalies
DLAB_calc_point_zero_1989 = pointer(DLAB_cp, rgv.thresh=0, nseries.thresh=75)[72,] #83% of cores had positive calcification anomalies

#Figure 2####################
#Calcification time series plots w.r.t. NAO index

#establish common colors for each species
cols <- c("D. labyrinthiformis" = "#2c7bb6", "O. franksi" = "#fdae61", "P. strigosa" = "#abd9e9")

NAO_plot=
  ggplot() + 
  geom_hline(yintercept=0,colour="red",linetype="dashed", size=2)+
  geom_ribbon(data=NAO_winter_means,aes(x=year,ymin=NAO_winter_lwr, ymax=NAO_winter_upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=NAO_winter_means, aes(x=year, y=NAO_winter), color="black") +
  ylab("NAO index")+
  xlab("")+
  annotate("text", x=1930, y=2.5, label= expression("(a) DJFM mean"),size=6,hjust = 0)+
  coord_cartesian(ylim = c(-3,3),xlim=c(1930,2015))+
  scale_x_continuous(position="top",expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(-3,3,1))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_line(colour = "grey"),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank())

DLAB_calc_plot=
  ggplot() + 
  geom_segment(data=DLAB_calc_point, mapping=aes(x=Year, y=0, xend=Year, yend=0.8, color = Nature),size=2)+
  geom_line(data=subset(annual_growth, species=="DLAB"), aes(x=year, y=calc, group=coral), color = "grey") +
  geom_ribbon(data=meanDLAB_calc,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanDLAB_calc, aes(x=year, y=mean)) +
  ylab(expression(Calcification~rate))+
  annotate("text", x=1930, y=0.73, label= expression(textstyle(group("(",b,")"))~italic("D. labyrinthiformis")),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,0.8),xlim=c(1930,2015))+
  scale_x_continuous(position="bottom",expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(position="right",expand = c(0,0), breaks = seq(0,0.8,0.2))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_line(colour = "grey"),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.line.x = element_blank())

OFRA_calc_plot=
  ggplot() + 
  geom_line(data=subset(annual_growth, species=="OFRA"), aes(x=year, y=calc, group=coral), color = "grey") +
  geom_ribbon(data=meanOFRA_calc,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanOFRA_calc, aes(x=year, y=mean)) +
  ylab(expression(Calcification~rate))+
  xlab("")+
  annotate("text", x=1930, y=0.73, label= expression(textstyle(group("(",c,")"))~italic("O. franksi")),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,0.8),xlim=c(1930,2015))+
  scale_x_continuous(position="bottom",expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,0.8,0.2))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_line(colour = "grey"),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.line.x = element_blank())

PSTR_calc_plot=
  ggplot() + 
  geom_line(data=subset(annual_growth, species=="PSTR"), aes(x=year, y=calc, group=coral), color = "grey") +
  geom_ribbon(data=meanPSTR_calc,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanPSTR_calc, aes(x=year, y=mean)) +
  ylab(expression(Calcification~rate))+
  annotate("text", x=1930, y=0.73, label= expression(textstyle(group("(",d,")"))~italic("P. strigosa")),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,0.8),xlim=c(1930,2015))+
  scale_x_continuous(expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(position="right",expand = c(0,0), breaks = seq(0,0.8,0.2))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_line(colour = "grey"),
        panel.grid.major.y = element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_blank())

nCoresplot=
  ggplot() + 
  geom_line(data=DLAB_pointer, aes(x=Year, y=Nb.series,colour="D. labyrinthiformis"), size = 1) +
  geom_line(data=OFRA_pointer, aes(x=Year, y=Nb.series,colour="O. franksi"), size = 1) +
  geom_line(data=PSTR_pointer, aes(x=Year, y=Nb.series,colour="P. strigosa"), size = 1) +
  ylab("Number of cores")+
  annotate("text", x=1930, y=15.7, label= paste("(e)"),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,16),xlim=c(1930,2015))+
  scale_x_continuous(expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,15,1))+
  scale_colour_manual(name="",values=cols, guide = guide_legend(override.aes=aes(fill=NA))) +
  theme_classic()+
  theme(text = element_text(size=16),
        #legend.position = "top",
        legend.position = c(0.24, 0.86),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_line(colour = "grey"),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank())

Figure2 = (NAO_plot/DLAB_calc_plot/OFRA_calc_plot/PSTR_calc_plot)|nCoresplot
Figure2
ggsave("Figure2.pdf", plot=Figure2, width = 10, height = 8, units = "in", dpi = 300)

#Figure 3####################
Calc_species_plot=
  ggplot() +
  geom_hline(yintercept=0,colour="red",linetype="dashed",size=2)+
  geom_linerange(data=as.data.frame(c13_species_int$emtrends), aes(x = species, ymin = lower.CL, ymax = upper.CL), size = 3, color="gray") + 
  geom_point(data=as.data.frame(c13_species_int$emtrends), aes(x = species, y = NAO_winter.trend), size=6, color="black", fill="white", shape=22) + 
  annotate("text", x=0.45, y=0.037, label= "(a) Species",size=6,hjust=0)+
  annotate("text", x=0.965, y=-0.001, label= "x",size=6,hjust=0)+
  annotate("text", x=1.965, y=0.018, label= "y",size=6,hjust=0)+
  annotate("text", x=2.965, y=0.0185, label= "y",size=6,hjust=0)+
  scale_x_discrete(labels= c("D. labyrinthiformis","O. franksi","P. strigosa"))+
  ylab("Slope of Calcification vs. NAO (DJFM)")+
  coord_cartesian(ylim = c(-0.03, 0.04))+
  scale_y_continuous(breaks = seq(-0.03,0.04,0.01), expand = c(0, 0))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black",face="italic"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x=element_blank(),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

Calc_site_plot=
  ggplot() +
  geom_hline(yintercept=0,colour="red",linetype="dashed",size=2)+
  geom_linerange(data=as.data.frame(c13_site_int$emtrends), aes(x = site, ymin = lower.CL, ymax = upper.CL), size = 3, color="gray") + 
  geom_point(data=as.data.frame(c13_site_int$emtrends), aes(x = site, y = NAO_winter.trend), size=6, color="black", fill="white", shape=22) + 
  annotate("text", x=0.45, y=0.037, label= "(b) Site",size=6,hjust=0)+
  annotate("text", x=3.86, y=0.0055, label= "xyz",size=6,hjust=0)+
  annotate("text", x=4.91, y=0.0345, label= "xy",size=6,hjust=0)+
  annotate("text", x=2.92, y=0.0035, label= "xz",size=6,hjust=0)+
  annotate("text", x=1.86, y=0.0205, label= "xyz",size=6,hjust=0)+
  annotate("text", x=0.84, y=0.0175, label= "xyz",size=6,hjust=0)+
  scale_x_discrete(labels= c("Gurnet\nRock","Halfway\nFlats","Hog\nReef","Three Hill\nShoals","Whalebone\nBay"))+
  ylab("Slope of Calcification vs. NAO (DJFM)")+
  coord_cartesian(ylim = c(-0.03, 0.04))+
  scale_y_continuous(position="right",breaks = seq(-0.03,0.04,0.01), expand = c(0, 0))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x=element_blank(),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

Figure3 = (Calc_species_plot | Calc_site_plot)
ggsave("Figure3.pdf", plot=Figure3 , width = 12, height = 5, units = "in", dpi = 300)

#Figure S1####################
Scan_1_plot=
  ggplot() + 
  geom_point(data=standards_scan_1, aes(intensity,density),size=3, alpha = 0.5) +
  annotate("text", x=141, y=2.35, label= "(a) Scan 1",size=6)+
  geom_line(data=modelfit1,aes(x=intensity,y=fit))+
  geom_ribbon(data=modelfit1,aes(x=intensity,ymin = lwr, ymax = upr), alpha = .15)+
  ylab(expression(Density~(g~cm^-3)))+xlab("Luminance")+
  coord_cartesian(ylim = c(0, 2.5))+
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5),expand = c(0, 0))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

Scan_2_plot=
  ggplot() + 
  geom_point(data=standards_scan_2, aes(intensity,density),size=3, alpha = 0.5) +
  annotate("text", x=141, y=2.35, label= "(b) Scan 2",size=6)+
  geom_line(data=modelfit2,aes(x=intensity,y=fit))+
  geom_ribbon(data=modelfit2,aes(x=intensity,ymin = lwr, ymax = upr), alpha = .15)+
  ylab("")+
  xlab("Luminance")+
  coord_cartesian(ylim = c(0, 2.5))+
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5),expand = c(0, 0))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

Scan_3_plot=
  ggplot() + 
  geom_point(data=standards_scan_3, aes(intensity,density),size=3, alpha = 0.5) +
  annotate("text", x=141, y=2.35, label= "(c) Scan 3",size=6)+
  geom_line(data=modelfit3,aes(x=intensity,y=fit))+
  geom_ribbon(data=modelfit3,aes(x=intensity,ymin = lwr, ymax = upr), alpha = .15)+
  ylab("")+
  xlab("Luminance")+
  coord_cartesian(ylim = c(0, 2.5))+
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5),expand = c(0, 0))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

FigureS1 = (Scan_1_plot | Scan_2_plot | Scan_3_plot)

ggsave("FigureS1.pdf", plot=FigureS1 , width = 13.5, height = 4, units = "in", dpi = 300)

#Figure S2####################
tissue_plot=
  ggplot() +
  geom_linerange(data=mean_core_metadata, aes(x = species, ymin = tissue_mean-tissue_sd, ymax = tissue_mean+tissue_sd), size = 3, alpha = 0.5) + 
  geom_point(data=mean_core_metadata, aes(x = species, y = tissue_mean), size=6, color="black", fill="white", shape=22) + 
  geom_jitter(data=metadata_growth, aes(x = species, y = tissue, colour = site), width=0.15,alpha=0.4, size = 4) +
  xlab("")+annotate("text", x=0.55, y=10.5, label= "(a)",size=6)+
  scale_x_discrete(labels= c("D. labyrinthiformis","O. franksi","P. strigosa"))+
  ylab("Tissue Thickness (mm)")+
  coord_cartesian(ylim = c(0, 11))+
  scale_y_continuous(breaks = c(2,4,6,8,10),expand = c(0, 0))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black",face="italic"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

ext_plot=
  ggplot() +
  geom_linerange(data=mean_core_metadata, aes(x = species, ymin = ext_mean-ext_sd, ymax = ext_mean+ext_sd), size = 3, alpha = 0.5) + 
  geom_point(data=mean_core_metadata, aes(x = species, y = ext_mean), size=6, color="black", fill="white", shape=22) + 
  geom_jitter(data=metadata_growth, aes(x = species, y = ext, colour = site), width=0.15,alpha=0.4, size = 4) +
  xlab("")+annotate("text", x=0.55, y=0.57, label= "(b)",size=6)+
  scale_x_discrete(labels= c("D. labyrinthiformis","O. franksi","P. strigosa"))+
  ylab(expression(Extension~(cm~yr^-1)))+
  coord_cartesian(ylim = c(0, 0.6))+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black",face="italic"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

dens_plot=
  ggplot() +
  geom_linerange(data=mean_core_metadata, aes(x = species, ymin = dens_mean-dens_sd, ymax = dens_mean+dens_sd), size = 3, alpha = 0.5) + 
  geom_point(data=mean_core_metadata, aes(x = species, y = dens_mean), size=6, color="black", fill="white", shape=22) + 
  geom_jitter(data=metadata_growth, aes(x = species, y = dens, colour = site), width=0.15,alpha=0.4, size = 4) +
  xlab("")+annotate("text", x=0.55, y=2.38, label= "(c)",size=6)+
  scale_x_discrete(labels= c("D. labyrinthiformis","O. franksi","P. strigosa"))+
  ylab(expression(Density~(g~cm^-3)))+
  coord_cartesian(ylim = c(0, 2.5))+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black",face="italic"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

calc_plot=
  ggplot() +
  geom_linerange(data=mean_core_metadata, aes(x = species, ymin = calc_mean-calc_sd, ymax = calc_mean+calc_sd), size = 3, alpha = 0.5) + 
  geom_point(data=mean_core_metadata, aes(x = species, y = calc_mean), size=6, color="black", fill="white", shape=22) + 
  geom_jitter(data=metadata_growth, aes(x = species, y = calc, colour = site), width=0.15,alpha=0.4, size = 4) +
  xlab("")+annotate("text", x=0.55, y=0.76, label= "(d)",size=6)+
  scale_x_discrete(labels= c("D. labyrinthiformis","O. franksi","P. strigosa"))+
  ylab(expression(Calcification~(g~cm^-2~yr^-1)))+
  coord_cartesian(ylim = c(0, 0.8))+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black",face="italic"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

FigureS2 = (tissue_plot | ext_plot) / (dens_plot | calc_plot)

ggsave("FigureS2.pdf", plot=FigureS2 , width = 10, height = 7.5, units = "in", dpi = 300)

#Figure S3####################

#Extension time series plots
DLAB_ext_plot=
  ggplot() + 
  geom_line(data=subset(annual_growth, species=="DLAB"), aes(x=year, y=ext, group=coral), color = "grey") +
  geom_ribbon(data=meanDLAB_ext,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanDLAB_ext, aes(x=year, y=mean)) +
  ylab(expression(Extension~(cm~yr^-1)))+
  xlab("")+
  annotate("text", x=1930, y=0.94, label= expression("(a)"),size=6,hjust = 0)+
  annotate("text", x=1937, y=0.94, label= expression(paste(italic("D. labyrinthiformis"))),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,1),xlim=c(1930,2015))+
  scale_x_continuous(expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.2))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

OFRA_ext_plot=
  ggplot() + 
  geom_line(data=subset(annual_growth, species=="OFRA"), aes(x=year, y=ext, group=coral), color = "grey") +
  geom_ribbon(data=meanOFRA_ext,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanOFRA_ext, aes(x=year, y=mean)) +
  ylab("")+
  xlab("")+
  annotate("text", x=1930, y=0.94, label= expression("(b)"),size=6,hjust = 0)+
  annotate("text", x=1937, y=0.94, label= expression(paste(italic("O. franksi"))),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,1),xlim=c(1930,2015))+
  scale_x_continuous(expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.2))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

PSTR_ext_plot=
  ggplot() + 
  geom_line(data=subset(annual_growth, species=="PSTR"), aes(x=year, y=ext, group=coral), color = "grey") +
  geom_ribbon(data=meanPSTR_ext,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanPSTR_ext, aes(x=year, y=mean)) +
  ylab("")+
  xlab("")+
  annotate("text", x=1930, y=0.94, label= expression("(c)"),size=6,hjust = 0)+
  annotate("text", x=1937, y=0.94, label= expression(paste(italic("P. strigosa"))),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,1),xlim=c(1930,2015))+
  scale_x_continuous(expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.2))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

#Density time series plots
DLAB_dens_plot=
  ggplot() + 
  geom_line(data=subset(annual_growth, species=="DLAB"), aes(x=year, y=dens, group=coral), color = "grey") +
  geom_ribbon(data=meanDLAB_dens,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanDLAB_dens, aes(x=year, y=mean)) +
  ylab(expression(Density~(g~cm^-3)))+
  xlab("")+
  annotate("text", x=1930, y=2.35, label= expression("(d)"),size=6,hjust = 0)+
  annotate("text", x=1937, y=2.35, label= expression(paste(italic("D. labyrinthiformis"))),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,2.5),xlim=c(1930,2015))+
  scale_x_continuous(expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,2.5,0.5))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

OFRA_dens_plot=
  ggplot() + 
  geom_line(data=subset(annual_growth, species=="OFRA"), aes(x=year, y=dens, group=coral), color = "grey") +
  geom_ribbon(data=meanOFRA_dens,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanOFRA_dens, aes(x=year, y=mean)) +
  ylab("")+
  xlab("")+
  annotate("text", x=1930, y=2.35, label= expression("(e)"),size=6,hjust = 0)+
  annotate("text", x=1937, y=2.35, label= expression(paste(italic("O. franksi"))),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,2.5),xlim=c(1930,2015))+
  scale_x_continuous(expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,2.5,0.5))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

PSTR_dens_plot=
  ggplot() + 
  geom_line(data=subset(annual_growth, species=="PSTR"), aes(x=year, y=dens, group=coral), color = "grey") +
  geom_ribbon(data=meanPSTR_dens,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanPSTR_dens, aes(x=year, y=mean)) +
  ylab("")+
  xlab("")+
  annotate("text", x=1930, y=2.35, label= expression("(f)"),size=6,hjust = 0)+
  annotate("text", x=1937, y=2.35, label= expression(paste(italic("P. strigosa"))),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,2.5),xlim=c(1930,2015))+
  scale_x_continuous(expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,2.5,0.5))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

#Calcification time series plots
DLAB_calc_plot=
  ggplot() + 
  geom_segment(data=DLAB_calc_point, mapping=aes(x=Year, y=0, xend=Year, yend=0.8, color = Nature),size=2)+
  geom_line(data=subset(annual_growth, species=="DLAB"), aes(x=year, y=calc, group=coral), color = "grey") +
  geom_ribbon(data=meanDLAB_calc,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanDLAB_calc, aes(x=year, y=mean)) +
  ylab(expression(Calcification~rate~(g~cm^-2~yr^-1)))+
  xlab("")+
  annotate("text", x=1930, y=0.752, label= expression("(g)"),size=6,hjust = 0)+
  annotate("text", x=1937, y=0.752, label= expression(paste(italic("D. labyrinthiformis"))),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,0.8),xlim=c(1930,2015))+
  scale_x_continuous(expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,0.8,0.1))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

OFRA_calc_plot=
  ggplot() + 
  geom_line(data=subset(annual_growth, species=="OFRA"), aes(x=year, y=calc, group=coral), color = "grey") +
  geom_ribbon(data=meanOFRA_calc,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanOFRA_calc, aes(x=year, y=mean)) +
  ylab("")+
  xlab("")+
  annotate("text", x=1930, y=0.752, label= expression("(h)"),size=6,hjust = 0)+
  annotate("text", x=1937, y=0.752, label= expression(paste(italic("O. franksi"))),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,0.8),xlim=c(1930,2015))+
  scale_x_continuous(expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,0.8,0.1))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

PSTR_calc_plot=
  ggplot() + 
  geom_line(data=subset(annual_growth, species=="PSTR"), aes(x=year, y=calc, group=coral), color = "grey") +
  geom_ribbon(data=meanPSTR_calc,aes(x=year,ymin=lwr, ymax=upr), alpha=0.2, fill = "#56B4E9")+
  geom_line(data=meanPSTR_calc, aes(x=year, y=mean)) +
  ylab("")+
  xlab("")+
  annotate("text", x=1930, y=0.752, label= expression("(i)"),size=6,hjust = 0)+
  annotate("text", x=1937, y=0.752, label= expression(paste(italic("P. strigosa"))),size=6,hjust = 0)+
  coord_cartesian(ylim = c(0,0.8),xlim=c(1930,2015))+
  scale_x_continuous(expand = c(0,1), breaks = seq(1930,2020,10))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,0.8,0.1))+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.background = element_blank(),
        axis.ticks = element_blank())

FigureS3=(DLAB_ext_plot|OFRA_ext_plot|PSTR_ext_plot) /
  (DLAB_dens_plot|OFRA_dens_plot|PSTR_dens_plot) /
  (DLAB_calc_plot|OFRA_calc_plot|PSTR_calc_plot)

ggsave("FigureS3.pdf", plot=FigureS3, width = 14, height = 10, units = "in", dpi = 300)

#Data export for BCO-DMO####################
coring_locations=read_csv("site_coordinates.csv") #read in lat lon for each site
core_data=annual_growth[complete.cases(annual_growth),] #remove years with missing data
core_data$ext=format(round(core_data$ext,digits=2),digits=2) #round and format extension rates to 2 decimal places
core_data$dens=format(round(core_data$dens,digits=2),digits=2) #round and format skeletal density values to 2 decimal places
core_data$calc=format(round(core_data$calc,digits=2),digits=2) #round and format calcification rates to 2 decimal places
core_data_coords=merge(core_data,coring_locations,by="site",all.x=TRUE) #merge lat long with coral growth data for each site
core_data_coords$species=str_replace(core_data_coords$species,"DLAB","Diploria_labyrinthiformis") #rename species with full name
core_data_coords$species=str_replace(core_data_coords$species,"PSTR","Psuedodiploria_strigosa") #rename species with full name
core_data_coords$species=str_replace(core_data_coords$species,"OFRA","Orbicella_franksi") #rename species with full name
col_order <- c("year","ext","dens","calc","coral","site","species","site_lat","site_lon") #arrange order of columns
coresBCODMO=core_data_coords[order(core_data_coords$coral,core_data_coords$year),col_order] #sort data by core ID and year
write.csv(coresBCODMO,"coresBCODMO.csv",row.names=FALSE) #export .csv of data for archival in BCO-DMO
