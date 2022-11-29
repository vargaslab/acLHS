#### ------- RGeoestad script in 2D ------- ####
#--------------------------------------------------------------#
#                   Data Manipulation                          #
#--------------------------------------------------------------#

root_dir<-getwd()

# Reading Data File in ASCII format space separated (.txt)
# -999.25 for non available values
#Data_File <- read.table(file=file.choose(),header=TRUE,na.strings="-999.25")

# Reading Data File in ASCII format comma separated (.csv)
# -999.25 for non available values
Data_Frame <- read.csv(file=file.choose(),header=T,na.strings="-999.25")

# Create a folder to store results for a case
dir.create(paste(getwd(),"/Results/Time_Series", sep=""))

result_dir<-paste(root_dir,"/Results/Time_Series",sep="")

# Create a folder to store results for the Exploratory Data Analysis (EDA)
dir.create(paste(getwd(),"/Results/Time_Series/EDA", sep=""))

aed_dir<-paste(result_dir,"/EDA",sep="")

Time<-Data_Frame$Time
Temp<- Data_Frame$Temp
CO2<- Data_Frame$CO2

Timelabel <- expression(bold("Time [days]")) 
CO2label <- expression(bold(paste(CO[2] ~ efflux ~ "["~μmol ~ m^-2 ~ s^-1~"]"))) 
Templabel <- expression(bold("Temperature [°C]"))
  
#--------------------------------------------------------------#
#               Exploratory Data Analysis                      #
#--------------------------------------------------------------#
###------------ Univariate Data Analysis---------------------###

# Basic Statistics Estimation 

# Basic Statistics
Time_Stat<-Estadisticas(Time)
CO2_Stat<-Estadisticas(CO2)
Temp_Stat<-Estadisticas(Temp)
#n_bins=trunc(sqrt(CO2_Stat[1,2]))
n_bins = nclass.Sturges(CO2)
# nclass.Sturges(Temp)
# nclass.Sturges(Time)
#n_bins = 10

#--------------CO2 data-----------------#
# Histogram and Boxplot
# Save plot as image in png format
png(paste(aed_dir,"/CO2_HistBoxPlotCounts.png",sep=""), bg = "white", width = 1500, height = 1500, res = 250)
HistBoxplot(x=CO2, mean = CO2_Stat[5,2], median = CO2_Stat[4,2], main = "",  
            xlab = CO2label, ylab = expression(bold("Absolute frequency [count]")), AbsFreq = TRUE, PercentFreq = FALSE,
            nbin = n_bins)
dev.off()
png(paste(aed_dir,"/CO2_HistBoxPlotFreq.png",sep=""), bg = "white",  width = 1500, height = 1500, res = 250)
histf_CO2 <- HistBoxplot(x=CO2, mean = CO2_Stat[5,2], median = CO2_Stat[4,2], main ="", 
            xlab = CO2label, ylab = expression(bold("Relative frequency [%]")), AbsFreq = FALSE, PercentFreq = TRUE,
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
# CO2 is the dependent variable (y-axis)
png(paste(aed_dir,"/CO2-Temp_ScatterPlot.png",sep=""), bg = "white",  width = 1500, height = 1500, res = 220)
ScatterPlot(Temp, CO2,  n_bins, 
            Ymin = CO2_Stat[2,2], Ymax = CO2_Stat[7,2], 
            Xmin = Temp_Stat[2,2], Xmax = Temp_Stat[7,2], 
            YLAB = CO2label, XLAB = Templabel)
dev.off()
