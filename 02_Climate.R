# Authors: Jacob Moutouama
# Date last modified (Y-M-D): 2024-08-03
rm(list = ls())
# load packages
# install.packages("xts")
# install.packages("tsbox")
# install.packages("forecast")
library(ggplot2)
library(ggsci)
# load the data
Ker<-read.csv("https://www.dropbox.com/scl/fi/y6jw1ildsp95y25lethwt/Ker.csv?rlkey=ukcw1ngt6i23y98ubtfm9to8e&dl=1",as.is=T) 
View(Ker)

figtempsite<-ggplot(Ker, aes(x=as.Date(DATE, format= "%Y - %m - %d"), y=TMAX))+
  scale_fill_jco()+
  theme_bw()

+ 
  theme(legend.position = "none",
        axis.text.x = element_text(size=4.5,color="black", angle=0))+
  labs( y="Daily temperature  (°C)", x="")


# convert to a Time Series
summary(Ker)


Ker_ts<-xts(Ker[,c("TMAX","TMIN","PRCP")], order.by=as.Date(Ker$DATE))
Ker_ts <-ts_regular(Ker_ts)

historical = na.fill(historical, "extend")

historical = window(historical, start=as.Date("2000-01-01"), end=as.Date("2020-12-31"))