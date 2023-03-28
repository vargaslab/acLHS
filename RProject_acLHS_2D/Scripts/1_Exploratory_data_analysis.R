#### ------- RGeoestad script in 2D ------- ####
#--------------------------------------------------------------#
#                   Data Manipulation                          #
#--------------------------------------------------------------#

root_dir<-getwd()

# Creates a folder to store results for a case
dir.create(paste(getwd(),"/Results/spatial_2D", sep=""))

result_dir<-paste(root_dir,"/Results/spatial_2D",sep="")

# Creates a folder to store results for the Exploratory Data Analysis (EDA)
dir.create(paste(getwd(),"/Results/spatial_2D/EDA", sep=""))

aed_dir<-paste(result_dir,"/EDA",sep="")


# Reading Data File in ASCII format space separated (.txt)
# -999.25 for non available values
#Data_File <- read.table(file=file.choose(),header=TRUE,na.strings="-999.25")

# Reading Data File in ASCII format comma separated (.csv)
# -999.25 for non available values
Data_Frame <- read.csv(file=file.choose(),header=T,na.strings="-999.25")

X<-Data_Frame[,1]
Y<-Data_Frame[,2]
Temp<-Data_Frame[,3]    #Temp_2012_conus_laea_1km
Rs<-Data_Frame[,4]    #stell_CONUS_Rs_1km #layer  

Xlabel <- expression(bold("X [km]")) 
Ylabel <- expression(bold("Y [km]")) 
# Rslabel <- expression(bold(paste(C ~ flux ~ "["~ g ~ C ~ m^-2 ~"]"))) 
Rslabel <- expression(bold(paste(Rs ~ "["~ g ~ C ~ m^-2 ~"]"))) 
Templabel <- expression(bold("Temp [°C]"))
  
#--------------------------------------------------------------#
#               Exploratory Data Analysis                      #
#--------------------------------------------------------------#
###------------ Univariate Data Analysis---------------------###

# Basic Statistics Estimation 

# Basic Statistics
X_Stat<-Estadisticas(X)
Y_Stat<-Estadisticas(Y)
Rs_Stat<-Estadisticas(Rs)
Temp_Stat<-Estadisticas(Temp)
#n_bins=trunc(sqrt(Rs_Stat[1,2]))
n_bins = nclass.Sturges(Rs)
# nclass.Sturges(Temp)
# nclass.Sturges(X)
#n_bins = 10

#--------------Rs data-----------------#
# Histogram and Boxplot
# Save plot as image in png format
png(paste(aed_dir,"/Rs_HistBoxPlotCounts.png",sep=""), bg = "white", width = 1500, height = 1500, res = 250)
HistBoxplot(x=Rs, mean = Rs_Stat[5,2], median = Rs_Stat[4,2], main = "",  
            xlab = Rslabel, ylab = expression(bold("Absolute frequency [count]")), AbsFreq = TRUE, PercentFreq = FALSE,
            nbin = n_bins)
dev.off()
png(paste(aed_dir,"/Rs_HistBoxPlotFreq.png",sep=""), bg = "white",  width = 1500, height = 1500, res = 250)
histf_Rs <- HistBoxplot(x=Rs, mean = Rs_Stat[5,2], median = Rs_Stat[4,2], main ="", 
            xlab = Rslabel, ylab = expression(bold("Relative frequency [%]")), AbsFreq = FALSE, PercentFreq = TRUE,
            nbin =n_bins)
dev.off()

#-----------Temp data-----------------#
# Basic Statistics
png(paste(aed_dir,"/Temp_HistBoxPlotFreq.png",sep=""), bg = "white",  width = 1500, height = 1500, res = 250)
HistBoxplot(x=Temp, mean = Temp_Stat[5,2], median = Temp_Stat[4,2], main ="", 
            xlab = Templabel, ylab = expression(bold("Relative frequency [%]")), AbsFreq = FALSE, PercentFreq = TRUE,
            nbin = n_bins)
dev.off()
png(paste(aed_dir,"/Temp_HistBoxPlotCounts.png",sep=""), bg = "white",  width = 1500, height = 1500, res = 250)
histf_Temp <- HistBoxplot(x=Temp, mean = Temp_Stat[5,2], median = Temp_Stat[4,2], main ="", 
            xlab = Templabel, ylab = expression(bold("Absolute frequency [count]")), AbsFreq = TRUE, PercentFreq = FALSE,
            nbin = n_bins)
dev.off()

###------------ Bivariate Data Analysis---------------------###
# Scatterplot with linear correlation coefficient
# Temp is the independent variable (x-axis)
# Rs is the dependent variable (y-axis)
png(paste(aed_dir,"/Rs-Temp_ScatterPlot.png",sep=""), bg = "white",  width = 1500, height = 1500, res = 220)
ScatterPlot(Temp, Rs,  n_bins, 
            Ymin = Rs_Stat[2,2], Ymax = Rs_Stat[7,2], 
            Xmin = Temp_Stat[2,2], Xmax = Temp_Stat[7,2], 
            YLAB = Rslabel, XLAB = Templabel)
dev.off()
