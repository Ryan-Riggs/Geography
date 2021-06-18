library(sf)
library(geosphere)
library(tidyverse)
library(sf)
library(sp)
library(ggplot2)
library(rgeos)
library(dplyr)
library(gmt)
library(vroom)
library(BAMMtools)
library(hydroGOF)
library(data.table)
library(parallel)
library(foreach)
library(reshape)
library(ggplot2)
require(gridExtra)
library(cowplot)
library(randomForest)

if (!"foreign" %in% rownames(installed.packages())){
  install.packages("foreign")}; require(foreign)
if (!"rgdal" %in% rownames(installed.packages())){
  install.packages("rgdal")}; require(rgdal)
if (!"shapefiles" %in% rownames(installed.packages())){
  install.packages("shapefiles")}; require(shapefiles)
if (!"RColorBrewer" %in% rownames(installed.packages())){
  install.packages("RColorBrewer")}; require(RColorBrewer)
if (!"zyp" %in% rownames(installed.packages())){
  install.packages("zyp")}; require(zyp)

##Error functions. 
source("E:/research/2019_08_30_rivObs/git/src/Error_stats_functions.R")


###Read in effective widths csv files. 
Eff_widths = map_df(list.files("E:\\research\\RatingCurveAnalysis\\Obs\\Gauges_30m_widths", full.names = TRUE, pattern = "In"), ~vroom(.x))
Eff_widths$Date = as.Date(as.POSIXct(Eff_widths$`system:time_start`/1000, origin = "1970-01-01"))

######################################################################################################################
###Set up validation data. 
Site_number_xsections = unique(Eff_widths$ID)

##################################################################################################
###Functions. 
auto.legend.pos <- function(x,y,xlim=NULL,ylim=NULL) {
  if (dev.cur() > 1) {
    p <- par('usr')
    if (is.null(xlim)) xlim <- p[1:2]
    if (is.null(ylim)) ylim <- p[3:4]
  } else {
    if (is.null(xlim)) xlim <- range(x, na.rm = TRUE)
    if (is.null(ylim)) ylim <- range(y, na.rm = TRUE)
  }
  countIt <- function(a) {
    tl <- sum(x <= xlim[1]*(1-a)+xlim[2]*a & y >= ylim[1]*a+ylim[2]*(1-a))
    tr <- sum(x >= xlim[1]*a+xlim[2]*(1-a) & y >= ylim[1]*a+ylim[2]*(1-a))
    bl <- sum(x <= xlim[1]*(1-a)+xlim[2]*a & y <= ylim[1]*(1-a)+ylim[2]*a)
    br <- sum(x >= xlim[1]*a+xlim[2]*(1-a) & y <= ylim[1]*(1-a)+ylim[2]*a)
    c(topleft=tl,topright=tr,bottomleft=bl,bottomright=br)
  }
  for (k in seq(0.5,0.1,by=-0.05)) {
    a <- countIt(k)
    if (sum(a==0)>0) break
    #if (!is.na(sum(a))) break
    
  }
  names(a)[which(a==0)][1]   # may delete "[1]"
}

is.error <- function(
  expr,
  tell=FALSE,
  force=FALSE
)
{
  expr_name <- deparse(substitute(expr))
  test <- try(expr, silent=TRUE)
  iserror <- inherits(test, "try-error")
  if(tell) if(iserror) message("Note in is.error: ", test)
  if(force) if(!iserror) stop(expr_name, " is not returning an error.", call.=FALSE)
  # output:
  iserror
}


validation = function(sim, obs){
  rrmse = RRMSE(sim, obs)
  nse = NSE(sim, obs)
  kge = KGE(sim, obs)
  nrmse = NRMSE(sim, obs)
  rbias = rBias(sim, obs)
  return(c(rrmse, nse, kge, nrmse, rbias))
}
usgs_q_processing = function(usgs_q){
  q_v = as.vector(usgs_q[,4])
  q_c = as.character(usgs_q[4])
  q_n = as.numeric(q_v)
  q= q_n *0.02832
  usgs_q = cbind(usgs_q, q)
  as.character(usgs_q$datetime)
  return(usgs_q)
}
#################################################################################################
start = as.Date("1979-01-01")
data_val = Eff_widths
RC_End = as.Date("2014-12-31") ###WAs 2014
RC_year = format(RC_End, "%Y")
RC_year_1 = format(RC_End+1, "%Y")
data_val$ID = data_val$ID
data_val$ID_2 = data_val$ID
data_val$calc_mean = data_val$Effective_width

#data_val$Date = as.character(data_val$Date)
tab$ID = tab$id
tab$ID_2 = tab$id
tab$width_m = data_val$width_m[match(tab$ID, data_val$ID)]
#tab$change= tab$median
#tab$change[mapply(is.na, tab$change)] <- 0
tab$width_m = data$width_m[match(tab$ID_2, data$ID_2)]
xSecq=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecw=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecIDcol=grep("V", names(Site_number_xsections))
mInd = array(5, dimnames = NULL)
rangedf_1 = as.data.frame(matrix(numeric(), nrow = 1, ncol = 4))
gage_stats = as.data.frame(matrix(numeric(), nrow =length(Site_number_xsections), ncol = 23))
gage_stats_GRADES = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 14))
colnames(gage_stats_GRADES)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "p_val","Bias", "RRMSE", "avg_std", "change", 'RRMSE_median', "std_Q", "STDE", "mode")
as.data.frame(gage_stats_GRADES)
l_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
Mean_grades = as.vector(nrow(Site_number_xsections))
u_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals_1 = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
width_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
gage_quants_q = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
gage_quants_w = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
colnames(gage_stats)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "mode","Bias", "RRMSE","avg_std", "change", "RRMSE_median", "std_Q","STDE", "KGE", "NSE", "rBias",
                        "SDRR", "MRR", "NRMSE", "Q_50", "W_50")
as.data.frame(gage_stats)
gage_stats_col1 = as.vector(1)
gage_stats_col2 = as.vector(1)
gage_stats_GRADES_col1 = as.vector(1)
gage_stats_GRADES_col2 = as.vector(1)
paired_df_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
rmse = 1
width_grouping = 30
percentiles = c(0.05, 0.95)
training = 0.7
data_val = as.data.table(data_val)
setkey(data_val, ID)

library(profvis)
library(caret)
library(dataRetrieval)
#profvis({
for(i in 1:length(Site_number_xsections)){
  id = Site_number_xsections[i]
  paired_df = data_val[.(Site_number_xsections[i])]
  paired_df_modis = paired_df[!is.na(paired_df$c),]
  paired_df_ls = paired_df[is.na(paired_df$c),]
  paired_df = paired_df_ls
  if(nchar(Site_number_xsections)<8){
    id = paste0("0", id)
  }
  usgs_q = try(usgs_q_processing(rawDailyData <- readNWISdv(as.character(id),"00060","1984-01-01"))) ##Usgs Make a switch for USGS vs Canadian. 
  if(!is.error(usgs_q)){
  paired_df = inner_join(paired_df, usgs_q)
  #paired_df = distinct(paired_df, Date, .keep_all = TRUE)
    #########################################################################################################    
  paired_df$Q = paired_df$q  
  paired_df = paired_df[!is.na(paired_df$Q),]
    
    if(nrow(paired_df)<2){next}
    set.seed(1)
    trainIndex = createDataPartition(paired_df$Q, p = training,
                                      list = FALSE)
    Train = paired_df[ trainIndex,]
    Valid = paired_df[-trainIndex,]
    
    if(nrow(Valid)<3){next}
    
    y = quantile(Train$Q, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    x = quantile(Train$calc_mean, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    if(length(unique(x))>1){
      ###Create function to estimate Q from a Landsat width based on quantile pair up. 
      spl = approxfun(x, y) #, method = "hyman")) #####either usse 'approxfun' or use splinefun with method = "hyman"
    } else{next}


      par(pty = "s")
      plot(c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)),c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)), type = "n", xlab = "In situ Discharge (cms)", ylab = "Landsat Discharge (cms)", main = "In situ vs Landsat")
      points(Valid$Q, spl(Valid$calc_mean), col = "blue")
      abline(0,1)
      
      names = c("rrmse", "nse", "kge", "nrmse", "rbias")
      
      ##Quantile validation 
      Valid$model = spl(Valid$calc_mean)
      qVal = validation(Valid$model, Valid$Q)
      qVal = as.data.frame(t(qVal))
      colnames(qVal) = names
      
      
      ##RF validation
      paired_df = merge(paired_df_ls, paired_df_modis,by = "Date", all.x = TRUE)#inner_join(paired_df, paired_df_modis, by = "Date")
      paired_df = paired_df[!is.na(paired_df$calc_mean.x)&!is.na(paired_df$c.y),]
      paired_df = merge(paired_df, usgs_q)
      paired_df = paired_df[!is.na(paired_df$q)]
      paired_df$Q = paired_df$q
      paired_df = paired_df%>%select("Q", "calc_mean.x", "c.y")
      c_value = mean(paired_df$c.y, na.rm = TRUE)
      w_value = mean(paired_df$calc_mean.x, na.rm = TRUE)
      
      paired_df = paired_df[!is.na(paired_df$c.y)&!is.na(paired_df$calc_mean.x),]
      
      paired_df$calc_mean.x[is.na(paired_df$calc_mean.x)]= w_value
      paired_df$c.y[is.na(paired_df$c.y)] = c_value
      
      # set.seed(5)
      # pd_preprocess = preProcess(paired_df, "knnImpute")
      # paired_df = predict(pd_preprocess, paired_df)
      # set.seed(2)
      
      if(nrow(paired_df)>1){
      trainIndex = createDataPartition(paired_df$Q, p = training,
                                       list = FALSE)
      Train = paired_df[ trainIndex,]
      Valid = paired_df[-trainIndex,]
      }else{next}
      if(nrow(Valid)<3){next}
      rf = randomForest(Q ~calc_mean.x+c.y,data = Train, ntree= 5000, mtry = 1)
      out = predict(rf, Valid)
      Valid$rf = out
      rfVal = validation(Valid$rf, Valid$Q)
      
      rfVal = as.data.frame(t(rfVal))
      colnames(rfVal) = names
      
      points(Valid$Q, Valid$rf, pch = 19)
      
      plsFit = try(train(Q ~calc_mean.x+c.y, 
                         data = Train,
                         method = "knn",
                         preProc=c("center", "scale"),
                         trControl = trainControl(method = "cv", number = 50, search = "random")))
        
        preds <- predict(plsFit, newdata = Valid)
        dlVal = validation(preds, Valid$Q)
        dlVal = as.data.frame(t(dlVal))
        colnames(dlVal) = names
        
      points(Valid$Q, preds, col = "red")
      best = c(qVal$nse, rfVal$nse, dlVal$nse)
      max(best)
      
      input = which(best==max(best))
      comb = list(qVal, rfVal, dlVal)
      out = as.data.frame(comb[input])
      gage_stats$mode[i] = input
      
      rrmse = out$rrmse 
      r = out$r
      nse = out$nse
      kge = out$kge
      nrmse = out$nrmse
      rbias = out$rbias
        
      bestFit = c("QNT", "RF", "KNN")
    
      legend1 = c(bestFit[input], paste0("RRMSE=",signif(rrmse, 3)),paste0("NSE=",signif(nse, 3)),
                  paste0("rBias=",signif(rbias, 3)), paste0("KGE=",signif(kge, 3)), paste0("NRMSE=",signif(nrmse, 3)))
      legend("topleft", legend1, xpd = TRUE, bty = "n", adj = 1, inset = c(-.05, -.1))
      

      gage_stats$RRMSE[i] = rrmse
      gage_stats$NSE[i] = nse
      gage_stats$KGE[i] = kge
      gage_stats$NRMSE[i] = nrmse
      gage_stats$rBias[i] = rbias


  } else{next}
}
#})
gage_stats_can = gage_stats
apply(gage_stats, 2, median, na.rm = TRUE)
length(na.omit(gage_stats$KGE))


ggplot(gage_stats, aes(x = 1:nrow(gage_stats), y = NSE, color = as.factor(mode)))+geom_point()






























grdc_path = "E:\\research\\RatingCurveAnalysis\\GaugeLocations\\DischargeDatasets\\GRDC\\"
grdc_files = gsub("_grdc", "_Q_Day.Cmd.txt", Site_number_xsections)
grdc_files = paste0(grdc_path, grdc_files)

for(i in 1:length(Site_number_xsections)){
  id = Site_number_xsections[i]
  paired_df = data_val[.(Site_number_xsections[i])]
  usgs_q = try(read.table(grdc_files[i], stringsAsFactors = FALSE))
  if(!is.error(usgs_q)){
    usgs_q$V1 = substr(usgs_q$V1, 0, 10)
    usgs_q$datetime = usgs_q$V1
    usgs_q$q = as.numeric(usgs_q$V2)
    usgs_q$Date = as.Date(usgs_q$datetime, format = "%Y-%m-%d")
    usgs_q = usgs_q[usgs_q$q>=0,]
    paired_df = inner_join(paired_df, usgs_q)
    paired_df$Q = paired_df$q
    
    if(nrow(paired_df)<3){next}
    set.seed(1)
    trainIndex = createDataPartition(paired_df$Q, p = training,
                                     list = FALSE)
    Train = paired_df[ trainIndex,]
    Valid = paired_df[-trainIndex,]
    
    if(nrow(Valid)<3){next}
    
    y = quantile(Train$Q, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    x = quantile(Train$calc_mean, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    if(length(unique(x))>1){
      ###Create function to estimate Q from a Landsat width based on quantile pair up. 
      spl = approxfun(x, y) #, method = "hyman")) #####either usse 'approxfun' or use splinefun with method = "hyman"
    } else{next}
    
    
    par(pty = "s")
    plot(c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)),c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)), type = "n", xlab = "In situ Discharge (cms)", ylab = "Landsat Discharge (cms)", main = "In situ vs Landsat")
    points(Valid$Q, spl(Valid$calc_mean), col = "blue")
    abline(0,1)
    
    names = c("rrmse", "nse", "kge", "nrmse", "rbias")
    
    ##Quantile validation 
    Valid$model = spl(Valid$calc_mean)
    qVal = validation(Valid$model, Valid$Q)
    qVal = as.data.frame(t(qVal))
    colnames(qVal) = names
    
    gage_stats$RRMSE[i] = qVal$rrmse
    gage_stats$NSE[i] = qVal$nse
    gage_stats$KGE[i] = qVal$kge
    gage_stats$NRMSE[i] = qVal$nrmse
    gage_stats$rBias[i] = qVal$rbias
    gage_stats$Site_number[i] = paired_df$ID[1]
    print(i)
  } else{next}
}

gage_file = st_read("E:\\research\\RatingCurveAnalysis\\GaugeLocations\\30_m_1984\\combined.shp")
gage_stats$GRWL_width_m = gage_file$GRWL_wd[match(gage_stats$Site_number, gage_file$Sttn_Nm)]

gage_stats = as.data.frame(gage_stats)
apply(gage_stats, 2, median, na.rm= TRUE)
gage_stats = as.data.frame(gage_stats)
boxplot(gage_stats$RRMSE[gage_stats$GRWL_width_m>90], ylim = c(0,200))

gage_file = gage_file[gage_file$Sttn_Nm%in%gage_stats$Site_number,]
gage_stats$Sttn_Nm = gage_stats$Site_number
gage_stats_vals = merge(gage_file, gage_stats,by = 'Sttn_Nm')

library(ggplot2)
library(tmap)

tmap_mode("view")

tm_shape(gage_stats_vals)+
  tm_bubbles(col = "NSE",size = 0.2, breaks = c(-1,-.5, 0,.5, 1))
tm_shape(gage_stats_vals)+
  tm_bubbles(col = "KGE",size = 0.02, breaks = c(-1, -.4, 1))
tm_shape(gage_stats_vals)+
  tm_bubbles(col = "mode",size = 0.02)


tm_shape(gage_stats_vals)+
  tm_bubbles(col = "rBias",size = 0.02, breaks = c(-100, 0, 100))
tm_shape(gage_stats_vals[gage_stats_vals$GRWL_wd>90,])+
  tm_bubbles(col = "RRMSE",size = 0.2, breaks = c(10, 35, 50, 75, 100, 150))


tm_shape(gage_stats_vals)+
  tm_bubbles(col = "n_Landsat_obs",size = 0.02, breaks = c(0, 100, 500, 1000, 2000))



gage_stats_vals1 = gage_stats_vals%>%select(mode, geometry, rBias, NRMSE, RRMSE, KGE, NSE, n_Landsat_obs, Site_number, lakeFlg, ds2GRWL, GRWL_wd)
st_write(gage_stats_vals1, "E:\\research\\RatingCurveAnalysis\\stats\\grdc_prelim1.shp", append = FALSE)


all=bind_rows(canada, brazil, thailand, russia, grdc)
tm_shape(all)+
  tm_bubbles(col = "NSE",size = 0.2, breaks = c(-1,-.5, 0,.5, 1))








##Knn and/or RF works the best it seems. 





########################################################################################################################
##Use best output for each gauge location. GRDc.  
########################################################################################################################
grdc_path = "E:\\research\\RatingCurveAnalysis\\GaugeLocations\\DischargeDatasets\\GRDC\\"
grdc_files = gsub("_grdc", "_Q_Day.Cmd.txt", Site_number_xsections)
grdc_files = paste0(grdc_path, grdc_files)



start = as.Date("1979-01-01")
data_val = Eff_widths
RC_End = as.Date("2014-12-31") ###WAs 2014
RC_year = format(RC_End, "%Y")
RC_year_1 = format(RC_End+1, "%Y")
data_val$ID = data_val$ID
data_val$ID_2 = data_val$ID
data_val$calc_mean = data_val$Effective_width

#data_val$Date = as.character(data_val$Date)
tab$ID = tab$id
tab$ID_2 = tab$id
tab$width_m = data_val$width_m[match(tab$ID, data_val$ID)]
#tab$change= tab$median
#tab$change[mapply(is.na, tab$change)] <- 0
tab$width_m = data$width_m[match(tab$ID_2, data$ID_2)]
xSecq=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecw=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecIDcol=grep("V", names(Site_number_xsections))
mInd = array(5, dimnames = NULL)
rangedf_1 = as.data.frame(matrix(numeric(), nrow = 1, ncol = 4))
gage_stats = as.data.frame(matrix(numeric(), nrow =length(Site_number_xsections), ncol = 23))
gage_stats_GRADES = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 14))
colnames(gage_stats_GRADES)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "p_val","Bias", "RRMSE", "avg_std", "change", 'RRMSE_median', "std_Q", "STDE", "mode")
as.data.frame(gage_stats_GRADES)
l_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
Mean_grades = as.vector(nrow(Site_number_xsections))
u_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals_1 = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
width_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
gage_quants_q = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
gage_quants_w = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
colnames(gage_stats)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "mode","Bias", "RRMSE","avg_std", "change", "RRMSE_median", "std_Q","STDE", "KGE", "NSE", "rBias",
                        "SDRR", "MRR", "NRMSE", "Q_50", "W_50")
as.data.frame(gage_stats)
gage_stats_col1 = as.vector(1)
gage_stats_col2 = as.vector(1)
gage_stats_GRADES_col1 = as.vector(1)
gage_stats_GRADES_col2 = as.vector(1)
paired_df_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
rmse = 1
width_grouping = 30
percentiles = c(0.05, 0.95)
training = 0.7
data_val = as.data.table(data_val)
setkey(data_val, ID)

library(profvis)
library(caret)
library(dataRetrieval)
#profvis({
for(i in 1:length(Site_number_xsections)){
  id = Site_number_xsections[i]
  paired_df = data_val[.(Site_number_xsections[i])]
  usgs_q = try(read.table(grdc_files[i], stringsAsFactors = FALSE))
  if(!is.error(usgs_q)){
    usgs_q$V1 = substr(usgs_q$V1, 0, 10)
    usgs_q$datetime = usgs_q$V1
    usgs_q$q = as.numeric(usgs_q$V2)
    usgs_q$Date = as.Date(usgs_q$datetime, format = "%Y-%m-%d")
    usgs_q = usgs_q[usgs_q$q>=0,]
    paired_df = inner_join(paired_df, usgs_q)
    paired_df$Q = paired_df$q
    
    paired_df = inner_join(paired_df, usgs_q)
    #paired_df = distinct(paired_df, Date, .keep_all = TRUE)
    #########################################################################################################    
    paired_df$Q = paired_df$q  
    paired_df = paired_df[!is.na(paired_df$Q),]
    
    if(nrow(paired_df)<3){next}
    set.seed(1)
    trainIndex = createDataPartition(paired_df$Q, p = training,
                                     list = FALSE)
    
    Train = paired_df[ trainIndex,]
    Valid = paired_df[-trainIndex,]
    # years = unique(year(paired_df$Date))
    # train = round(length(years)*.6)
    # paired_df$Year = year(paired_df$Date)
    # set.seed(5)
    # t = sample(years, train)
    # TrainingSet = paired_df$Year%in%t
    # Train = paired_df[TrainingSet,]
    # Valid = paired_df[!TrainingSet,]
    
    
    
    
    
    # Train = paired_df[paired_df$Date>=as.Date(as.character('2015-01-01', format = "%Y-%m-%d")),]
    # Valid = paired_df[paired_df$Date<as.Date(as.character('2015-01-01', format = "%Y-%m-%d")),]
    
    if(nrow(Valid)<3|nrow(Train)<3){next}
    
    y = quantile(Train$Q, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    x = quantile(Train$calc_mean, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    if(length(unique(x))>1){
      ###Create function to estimate Q from a Landsat width based on quantile pair up. 
      spl = approxfun(x, y) #, method = "hyman")) #####either usse 'approxfun' or use splinefun with method = "hyman"
    } else{next}
    
    
    # par(pty = "s")
    # plot(c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)),c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)), type = "n", xlab = "In situ Discharge (cms)", ylab = "Landsat Discharge (cms)", main = "In situ vs Landsat")
    # points(Valid$Q, spl(Valid$calc_mean), col = "blue")
    # abline(0,1)
    
    plot(paired_df$Date, paired_df$q, type = "l")
    points(Valid$Date, spl(Valid$calc_mean), col = "green")
    
    names = c("rrmse", "nse", "kge", "nrmse", "rbias")
    
    ##Quantile validation 
    Valid$model = spl(Valid$calc_mean)
    qVal = validation(Valid$model, Valid$Q)
    qVal = as.data.frame(t(qVal))
    colnames(qVal) = names
    
    
    rf = randomForest(Q ~calc_mean+Date,data = Train, ntree= nrow(Train), mtry = 1)
    out = predict(rf, Valid)
    Valid$rf = out
    rfVal = validation(Valid$rf, Valid$Q)
    
    rfVal = as.data.frame(t(rfVal))
    colnames(rfVal) = names
    
    #points(Valid$Q, Valid$rf, pch = 19)
    points(Valid$Date, Valid$rf, col = "blue")
    
    plsFit = try(train(Q ~calc_mean+Date, 
                       data = Train,
                       method = "knn",
                       preProc=c("center", "scale"),
                       trControl = trainControl(method = "cv", number = c(5, 10, 20, 50, 100),search = "random")))
    
    preds <- predict(plsFit, newdata = Valid)
    dlVal = validation(preds, Valid$Q)
    dlVal = as.data.frame(t(dlVal))
    colnames(dlVal) = names
    
    #points(Valid$Q, preds, col = "red")
    points(Valid$Date, preds, col = "red")
    best = c(qVal$nse, rfVal$nse, dlVal$nse)
    max(best, na.rm = TRUE)
    
    input = which(best==max(best, na.rm = TRUE))
    if(length(input) ==0){next}
    comb = list(qVal, rfVal, dlVal)
    out = as.data.frame(comb[input])
    gage_stats$mode[i] = input
    
    rrmse = out$rrmse 
    r = out$r
    nse = out$nse
    kge = out$kge
    nrmse = out$nrmse
    rbias = out$rbias
    
    bestFit = c("QNT", "RF", "KNN")
    
    legend1 = c(bestFit[input], paste0("RRMSE=",signif(rrmse, 3)),paste0("NSE=",signif(nse, 3)),
                paste0("rBias=",signif(rbias, 3)), paste0("KGE=",signif(kge, 3)), paste0("NRMSE=",signif(nrmse, 3)))
    legend("topleft", legend1, xpd = TRUE, bty = "n", adj = 1, inset = c(-.05, -.1))
    
    
    gage_stats$RRMSE[i] = rrmse
    gage_stats$NSE[i] = nse
    gage_stats$KGE[i] = kge
    gage_stats$NRMSE[i] = nrmse
    gage_stats$rBias[i] = rbias
    gage_stats$Site_number[i] = paired_df$ID[1]
    gage_stats$n_Landsat_obs[i] = nrow(paired_df)
    print(i)
    
  } else{next}
}








################################################################################################################################
##India gauges. 
################################################################################################################################
########################################################################################################################
India_discharge = list.files("E:\\research\\RatingCurveAnalysis\\GaugeLocations\\DischargeDatasets\\India\\", full.names = TRUE)
t = readxl::read_xlsx(India_discharge[1])
sites_file = fread("E:\\research\\RatingCurveAnalysis\\GaugeLocations\\India_sameOrderasQfiles.csv")
sites = as.data.frame(Site_number_xsections)
t[1,]
year = seq.Date(as.Date('1984-01-01'), as.Date('1984-12-31'), 1)
levelCols= grep("Level", t[2,])
t = t%>%select(-levelCols)
flowCols = grep("Flow", t[2,])
df = t[,flowCols]
df = df[3:nrow(df),]
colnames(df) = year
df = as.data.table(df)
df = sapply(df, as.numeric)
df$ID = sites_file$Sttn_Nm
df = melt(df, id.vars = "ID", measure.vars = year)
df$Date = rownames(df)
df$q = as.numeric(df$V1)






# t$station_id = gsim_all$gsim.no[match(t$`River Point Level & Flow Report (ALL AGENCIES) From 19840101 to 19841231`, gsim_all$station)]
# 

# colnames(India) = c("River", "Lat", "Long")
# India = India[3:nrow(India),]
# India$Sttn_Nm = paste0(1:nrow(India), "_In")
# 





start = as.Date("1979-01-01")
data_val = Eff_widths
RC_End = as.Date("2014-12-31") ###WAs 2014
RC_year = format(RC_End, "%Y")
RC_year_1 = format(RC_End+1, "%Y")
data_val$ID = data_val$ID
data_val$ID_2 = data_val$ID
data_val$calc_mean = data_val$Effective_width

#data_val$Date = as.character(data_val$Date)
tab$ID = tab$id
tab$ID_2 = tab$id
tab$width_m = data_val$width_m[match(tab$ID, data_val$ID)]
#tab$change= tab$median
#tab$change[mapply(is.na, tab$change)] <- 0
tab$width_m = data$width_m[match(tab$ID_2, data$ID_2)]
xSecq=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecw=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecIDcol=grep("V", names(Site_number_xsections))
mInd = array(5, dimnames = NULL)
rangedf_1 = as.data.frame(matrix(numeric(), nrow = 1, ncol = 4))
gage_stats = as.data.frame(matrix(numeric(), nrow =length(Site_number_xsections), ncol = 23))
gage_stats_GRADES = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 14))
colnames(gage_stats_GRADES)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "p_val","Bias", "RRMSE", "avg_std", "change", 'RRMSE_median', "std_Q", "STDE", "mode")
as.data.frame(gage_stats_GRADES)
l_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
Mean_grades = as.vector(nrow(Site_number_xsections))
u_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals_1 = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
width_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
gage_quants_q = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
gage_quants_w = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
colnames(gage_stats)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "mode","Bias", "RRMSE","avg_std", "change", "RRMSE_median", "std_Q","STDE", "KGE", "NSE", "rBias",
                        "SDRR", "MRR", "NRMSE", "Q_50", "W_50")
as.data.frame(gage_stats)
gage_stats_col1 = as.vector(1)
gage_stats_col2 = as.vector(1)
gage_stats_GRADES_col1 = as.vector(1)
gage_stats_GRADES_col2 = as.vector(1)
paired_df_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
rmse = 1
width_grouping = 30
percentiles = c(0.05, 0.95)
training = 0.7
data_val = as.data.table(data_val)
setkey(data_val, ID)

library(profvis)
library(caret)
library(dataRetrieval)
#profvis({
for(i in 1:length(Site_number_xsections)){
  id = Site_number_xsections[i]
  paired_df = data_val[.(Site_number_xsections[i])]
  usgs_q = try(read.table(grdc_files[i], stringsAsFactors = FALSE))
  if(!is.error(usgs_q)){
    usgs_q$V1 = substr(usgs_q$V1, 0, 10)
    usgs_q$datetime = usgs_q$V1
    usgs_q$q = as.numeric(usgs_q$V2)
    usgs_q$Date = as.Date(usgs_q$datetime, format = "%Y-%m-%d")
    usgs_q = usgs_q[usgs_q$q>=0,]
    paired_df = inner_join(paired_df, usgs_q)
    paired_df$Q = paired_df$q
    
    paired_df = inner_join(paired_df, usgs_q)
    #paired_df = distinct(paired_df, Date, .keep_all = TRUE)
    #########################################################################################################    
    paired_df$Q = paired_df$q  
    paired_df = paired_df[!is.na(paired_df$Q),]
    
    if(nrow(paired_df)<3){next}
    set.seed(1)
    trainIndex = createDataPartition(paired_df$Q, p = training,
                                     list = FALSE)
    Train = paired_df[ trainIndex,]
    Valid = paired_df[-trainIndex,]
    
    if(nrow(Valid)<3){next}
    
    y = quantile(Train$Q, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    x = quantile(Train$calc_mean, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    if(length(unique(x))>1){
      ###Create function to estimate Q from a Landsat width based on quantile pair up. 
      spl = approxfun(x, y) #, method = "hyman")) #####either usse 'approxfun' or use splinefun with method = "hyman"
    } else{next}
    
    
    par(pty = "s")
    plot(c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)),c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)), type = "n", xlab = "In situ Discharge (cms)", ylab = "Landsat Discharge (cms)", main = "In situ vs Landsat")
    points(Valid$Q, spl(Valid$calc_mean), col = "blue")
    abline(0,1)
    
    names = c("rrmse", "nse", "kge", "nrmse", "rbias")
    
    ##Quantile validation 
    Valid$model = spl(Valid$calc_mean)
    qVal = validation(Valid$model, Valid$Q)
    qVal = as.data.frame(t(qVal))
    colnames(qVal) = names
    
    
    rf = randomForest(Q ~calc_mean,data = Train, ntree= 5000, mtry = 1)
    out = predict(rf, Valid)
    Valid$rf = out
    rfVal = validation(Valid$rf, Valid$Q)
    
    rfVal = as.data.frame(t(rfVal))
    colnames(rfVal) = names
    
    points(Valid$Q, Valid$rf, pch = 19)
    
    plsFit = try(train(Q ~calc_mean, 
                       data = Train,
                       method = "knn",
                       preProc=c("center", "scale"),
                       trControl = trainControl(method = "cv", number = 10, search = "random")))
    
    preds <- predict(plsFit, newdata = Valid)
    dlVal = validation(preds, Valid$Q)
    dlVal = as.data.frame(t(dlVal))
    colnames(dlVal) = names
    
    points(Valid$Q, preds, col = "red")
    best = c(qVal$nse, rfVal$nse, dlVal$nse)
    max(best)
    
    input = which(best==max(best, na.rm = TRUE))
    if(length(input) ==0){next}
    comb = list(qVal, rfVal, dlVal)
    out = as.data.frame(comb[input])
    gage_stats$mode[i] = input
    
    rrmse = out$rrmse 
    r = out$r
    nse = out$nse
    kge = out$kge
    nrmse = out$nrmse
    rbias = out$rbias
    
    bestFit = c("QNT", "RF", "KNN")
    
    legend1 = c(bestFit[input], paste0("RRMSE=",signif(rrmse, 3)),paste0("NSE=",signif(nse, 3)),
                paste0("rBias=",signif(rbias, 3)), paste0("KGE=",signif(kge, 3)), paste0("NRMSE=",signif(nrmse, 3)))
    legend("topleft", legend1, xpd = TRUE, bty = "n", adj = 1, inset = c(-.05, -.1))
    
    
    gage_stats$RRMSE[i] = rrmse
    gage_stats$NSE[i] = nse
    gage_stats$KGE[i] = kge
    gage_stats$NRMSE[i] = nrmse
    gage_stats$rBias[i] = rbias
    gage_stats$Site_number[i] = paired_df$ID[1]
    gage_stats$n_Landsat_obs[i] = nrow(paired_df)
    print(i)
    
  } else{next}
}
################################################################################################################################
##Canada gauges. 
################################################################################################################################
########################################################################################################################
library(hydat)
library(tidyhydat)
gsim_all = read.csv("E:\\research\\GSIM\\GSIM_metadata\\GSIM_catalog\\GSIM_metadata.csv")
sites = as.data.frame(Site_number_xsections)
sites$gsim = gsub("_gsim", "", sites$Site_number_xsections)
sites$files= gsim_all$reference.no[match(sites$gsim, gsim_all$gsim.no)]




start = as.Date("1979-01-01")
data_val = Eff_widths
RC_End = as.Date("2014-12-31") ###WAs 2014
RC_year = format(RC_End, "%Y")
RC_year_1 = format(RC_End+1, "%Y")
data_val$ID = data_val$ID
data_val$ID_2 = data_val$ID
data_val$calc_mean = data_val$Effective_width

#data_val$Date = as.character(data_val$Date)
tab$ID = tab$id
tab$ID_2 = tab$id
tab$width_m = data_val$width_m[match(tab$ID, data_val$ID)]
#tab$change= tab$median
#tab$change[mapply(is.na, tab$change)] <- 0
tab$width_m = data$width_m[match(tab$ID_2, data$ID_2)]
xSecq=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecw=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecIDcol=grep("V", names(Site_number_xsections))
mInd = array(5, dimnames = NULL)
rangedf_1 = as.data.frame(matrix(numeric(), nrow = 1, ncol = 4))
gage_stats = as.data.frame(matrix(numeric(), nrow =length(Site_number_xsections), ncol = 23))
gage_stats_GRADES = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 14))
colnames(gage_stats_GRADES)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "p_val","Bias", "RRMSE", "avg_std", "change", 'RRMSE_median', "std_Q", "STDE", "mode")
as.data.frame(gage_stats_GRADES)
l_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
Mean_grades = as.vector(nrow(Site_number_xsections))
u_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals_1 = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
width_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
gage_quants_q = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
gage_quants_w = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
colnames(gage_stats)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "mode","Bias", "RRMSE","avg_std", "change", "RRMSE_median", "std_Q","STDE", "KGE", "NSE", "rBias",
                        "training", "validation", "NRMSE", "Q_50", "W_50")
as.data.frame(gage_stats)
gage_stats_col1 = as.vector(1)
gage_stats_col2 = as.vector(1)
gage_stats_GRADES_col1 = as.vector(1)
gage_stats_GRADES_col2 = as.vector(1)
paired_df_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
rmse = 1
width_grouping = 30
percentiles = c(0.05, 0.95)
training = 0.7
data_val = as.data.table(data_val)
setkey(data_val, ID)

library(profvis)
library(caret)
library(dataRetrieval)
library(lubridate)
#profvis({
for(i in 1:length(Site_number_xsections)){
  id = Site_number_xsections[i]
  paired_df = data_val[.(Site_number_xsections[i])]
  usgs_q = try(hy_daily(sites$files[i]))
  usgs_q = as.data.frame(usgs_q)
  if(nrow(usgs_q)>0){
    usgs_q$q = as.numeric(usgs_q$Value)
    usgs_q = usgs_q[usgs_q$q>=0,]
    #usgs_q_level  = usgs_q[usgs_q$Parameter =="Level",]
    usgs_q= usgs_q[usgs_q$Parameter=="Flow",]
    # colnames(usgs_q_level) = c(colnames(usgs_q[1:5]), "level")
    # usgs_q = left_join(usgs_q, usgs_q_level, by = "Date")
    # usgs_q = usgs_q%>%select(STATION_NUMBER.x, Date, q, level)
    paired_df = inner_join(paired_df, usgs_q, by="Date")
    paired_df$Q = paired_df$q

    paired_df = paired_df[paired_df$q>0,]
    #########################################################################################################    
    paired_df$Q = paired_df$q  
    paired_df = paired_df[!is.na(paired_df$Q),]
    
    if(nrow(paired_df)<3){next}
    set.seed(1)
    trainIndex = createDataPartition(paired_df$Q, p = training,
                                     list = FALSE)

    Train = paired_df[ trainIndex,]
    Valid = paired_df[-trainIndex,]
    # years = unique(year(paired_df$Date))
    # train = round(length(years)*.6)
    # paired_df$Year = year(paired_df$Date)
    # set.seed(5)
    # t = sample(years, train)
    # TrainingSet = paired_df$Year%in%t
    # Train = paired_df[TrainingSet,]
    # Valid = paired_df[!TrainingSet,]

    
    
    
    
    # Train = paired_df[paired_df$Date>=as.Date(as.character('2015-01-01', format = "%Y-%m-%d")),]
    # Valid = paired_df[paired_df$Date<as.Date(as.character('2015-01-01', format = "%Y-%m-%d")),]
    
    if(nrow(Valid)<3|nrow(Train)<3){next}
    
    y = quantile(Train$Q, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    x = quantile(Train$calc_mean, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    if(length(unique(x))>1){
      ###Create function to estimate Q from a Landsat width based on quantile pair up. 
      spl = approxfun(x, y) #, method = "hyman")) #####either usse 'approxfun' or use splinefun with method = "hyman"
    } else{next}
    
    
    # par(pty = "s")
    # plot(c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)),c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)), type = "n", xlab = "In situ Discharge (cms)", ylab = "Landsat Discharge (cms)", main = "In situ vs Landsat")
    # points(Valid$Q, spl(Valid$calc_mean), col = "blue")
    # abline(0,1)
    
    plot(paired_df$Date, paired_df$q, type = "l")
    points(Valid$Date, spl(Valid$calc_mean), col = "green")
    
    names = c("rrmse", "nse", "kge", "nrmse", "rbias")
    
    ##Quantile validation 
    Valid$model = spl(Valid$calc_mean)
    qVal = validation(Valid$model, Valid$Q)
    qVal = as.data.frame(t(qVal))
    colnames(qVal) = names
    
    
    rf = randomForest(Q ~calc_mean+Date,data = Train, ntree= nrow(Train), mtry = 1)
    out = predict(rf, Valid)
    Valid$rf = out
    rfVal = validation(Valid$rf, Valid$Q)
    
    rfVal = as.data.frame(t(rfVal))
    colnames(rfVal) = names
    
    #points(Valid$Q, Valid$rf, pch = 19)
    points(Valid$Date, Valid$rf, col = "blue")
    
    plsFit = try(train(Q ~calc_mean+Date, 
                       data = Train,
                       method = "knn",
                       preProc=c("center", "scale"),
                       trControl = trainControl(method = "cv", number = c(5, 10, 20, 50, 100),search = "random")))
    
    preds <- predict(plsFit, newdata = Valid)
    dlVal = validation(preds, Valid$Q)
    dlVal = as.data.frame(t(dlVal))
    colnames(dlVal) = names
    
    #points(Valid$Q, preds, col = "red")
    points(Valid$Date, preds, col = "red")
    best = c(qVal$nse, rfVal$nse, dlVal$nse)
    max(best, na.rm = TRUE)
    
    input = which(best==max(best, na.rm = TRUE))
    if(length(input) ==0){next}
    comb = list(qVal, rfVal, dlVal)
    out = as.data.frame(comb[input])
    gage_stats$mode[i] = input
    
    rrmse = out$rrmse 
    r = out$r
    nse = out$nse
    kge = out$kge
    nrmse = out$nrmse
    rbias = out$rbias
    
    bestFit = c("QNT", "RF", "KNN")
    
    legend1 = c(bestFit[input], paste0("RRMSE=",signif(rrmse, 3)),paste0("NSE=",signif(nse, 3)),
                paste0("rBias=",signif(rbias, 3)), paste0("KGE=",signif(kge, 3)), paste0("NRMSE=",signif(nrmse, 3)))
    legend("topleft", legend1, xpd = TRUE, bty = "n", adj = 1, inset = c(-.05, -.1))
    
    
    gage_stats$RRMSE[i] = rrmse
    gage_stats$NSE[i] = nse
    gage_stats$KGE[i] = kge
    gage_stats$NRMSE[i] = nrmse
    gage_stats$rBias[i] = rbias
    gage_stats$Site_number[i] = paired_df$ID[1]
    gage_stats$n_Landsat_obs[i] = nrow(paired_df)
    gage_stats$validation[i] = nrow(Valid)
    gage_stats$training[i] = nrow(Train)
    print(i)
    
  } else{next}
}
################################################################################################################################
##Russian gauges. 
################################################################################################################################
########################################################################################################################
RussiaQ = fread("E:\\research\\RatingCurveAnalysis\\GaugeLocations\\DischargeDatasets\\Russia\\Russian_daily_Q.txt", fill = TRUE)
RussiaQ = RussiaQ[RussiaQ$Year>=1984,]
tab = melt(RussiaQ, id.vars = c("Code", "Day", "Year"), measure.vars = c(month.abb))
tab$Month = match(tab$variable, month.abb)
tab$Date = paste(tab$Month, tab$Day, tab$Year, sep = "/")
tab$Date = as.Date(tab$Date, format = "%m/%d/%Y")

gsim_all = read.csv("E:\\research\\GSIM\\GSIM_metadata\\GSIM_catalog\\GSIM_metadata.csv")
sites = as.data.frame(Site_number_xsections)
sites$gsim = gsub("_gsim", "", sites$Site_number_xsections)
sites$files= gsim_all$reference.no[match(sites$gsim, gsim_all$gsim.no)]




start = as.Date("1979-01-01")
data_val = Eff_widths
RC_End = as.Date("2014-12-31") ###WAs 2014
RC_year = format(RC_End, "%Y")
RC_year_1 = format(RC_End+1, "%Y")
data_val$ID = data_val$ID
data_val$ID_2 = data_val$ID
data_val$calc_mean = data_val$Effective_width

#data_val$Date = as.character(data_val$Date)
tab$ID = tab$id
tab$ID_2 = tab$id
tab$width_m = data_val$width_m[match(tab$ID, data_val$ID)]
#tab$change= tab$median
#tab$change[mapply(is.na, tab$change)] <- 0
tab$width_m = data$width_m[match(tab$ID_2, data$ID_2)]
xSecq=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecw=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecIDcol=grep("V", names(Site_number_xsections))
mInd = array(5, dimnames = NULL)
rangedf_1 = as.data.frame(matrix(numeric(), nrow = 1, ncol = 4))
gage_stats = as.data.frame(matrix(numeric(), nrow =length(Site_number_xsections), ncol = 23))
gage_stats_GRADES = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 14))
colnames(gage_stats_GRADES)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "p_val","Bias", "RRMSE", "avg_std", "change", 'RRMSE_median', "std_Q", "STDE", "mode")
as.data.frame(gage_stats_GRADES)
l_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
Mean_grades = as.vector(nrow(Site_number_xsections))
u_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals_1 = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
width_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
gage_quants_q = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
gage_quants_w = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
colnames(gage_stats)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "mode","Bias", "RRMSE","avg_std", "change", "RRMSE_median", "std_Q","STDE", "KGE", "NSE", "rBias",
                        "SDRR", "MRR", "NRMSE", "Q_50", "W_50")
as.data.frame(gage_stats)
gage_stats_col1 = as.vector(1)
gage_stats_col2 = as.vector(1)
gage_stats_GRADES_col1 = as.vector(1)
gage_stats_GRADES_col2 = as.vector(1)
paired_df_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
rmse = 1
width_grouping = 30
percentiles = c(0.05, 0.95)
training = 0.7
data_val = as.data.table(data_val)
setkey(data_val, ID)

library(profvis)
library(caret)
library(dataRetrieval)
#profvis({
for(i in 1:length(Site_number_xsections)){
  id = Site_number_xsections[i]
  paired_df = data_val[.(Site_number_xsections[i])]
  usgs_q = try(tab[tab$Code==sites$files[i],])
  if(nrow(usgs_q)>0){
    usgs_q$q = as.numeric(usgs_q$value)
    usgs_q = usgs_q[usgs_q$q>=0,]
    paired_df = inner_join(paired_df, usgs_q)
    paired_df$Q = paired_df$q
    
    paired_df = inner_join(paired_df, usgs_q)
    #paired_df = distinct(paired_df, Date, .keep_all = TRUE)
    #########################################################################################################    
    paired_df$Q = paired_df$q  
    paired_df = paired_df[!is.na(paired_df$Q),]
    
    if(nrow(paired_df)<3){next}
    set.seed(1)
    trainIndex = createDataPartition(paired_df$Q, p = training,
                                     list = FALSE)
    
    Train = paired_df[ trainIndex,]
    Valid = paired_df[-trainIndex,]
    # years = unique(year(paired_df$Date))
    # train = round(length(years)*.6)
    # paired_df$Year = year(paired_df$Date)
    # set.seed(5)
    # t = sample(years, train)
    # TrainingSet = paired_df$Year%in%t
    # Train = paired_df[TrainingSet,]
    # Valid = paired_df[!TrainingSet,]
    
    
    
    
    
    # Train = paired_df[paired_df$Date>=as.Date(as.character('2015-01-01', format = "%Y-%m-%d")),]
    # Valid = paired_df[paired_df$Date<as.Date(as.character('2015-01-01', format = "%Y-%m-%d")),]
    
    if(nrow(Valid)<3|nrow(Train)<3){next}
    
    y = quantile(Train$Q, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    x = quantile(Train$calc_mean, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    if(length(unique(x))>1){
      ###Create function to estimate Q from a Landsat width based on quantile pair up. 
      spl = approxfun(x, y) #, method = "hyman")) #####either usse 'approxfun' or use splinefun with method = "hyman"
    } else{next}
    
    
    # par(pty = "s")
    # plot(c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)),c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)), type = "n", xlab = "In situ Discharge (cms)", ylab = "Landsat Discharge (cms)", main = "In situ vs Landsat")
    # points(Valid$Q, spl(Valid$calc_mean), col = "blue")
    # abline(0,1)
    
    plot(paired_df$Date, paired_df$q, type = "l")
    points(Valid$Date, spl(Valid$calc_mean), col = "green")
    
    names = c("rrmse", "nse", "kge", "nrmse", "rbias")
    
    ##Quantile validation 
    Valid$model = spl(Valid$calc_mean)
    qVal = validation(Valid$model, Valid$Q)
    qVal = as.data.frame(t(qVal))
    colnames(qVal) = names
    
    
    rf = randomForest(Q ~calc_mean+Date,data = Train, ntree= nrow(Train), mtry = 1)
    out = predict(rf, Valid)
    Valid$rf = out
    rfVal = validation(Valid$rf, Valid$Q)
    
    rfVal = as.data.frame(t(rfVal))
    colnames(rfVal) = names
    
    #points(Valid$Q, Valid$rf, pch = 19)
    points(Valid$Date, Valid$rf, col = "blue")
    
    plsFit = try(train(Q ~calc_mean+Date, 
                       data = Train,
                       method = "knn",
                       preProc=c("center", "scale"),
                       trControl = trainControl(method = "cv", number = c(5, 10, 20, 50, 100),search = "random")))
    
    preds <- predict(plsFit, newdata = Valid)
    dlVal = validation(preds, Valid$Q)
    dlVal = as.data.frame(t(dlVal))
    colnames(dlVal) = names
    
    #points(Valid$Q, preds, col = "red")
    points(Valid$Date, preds, col = "red")
    best = c(qVal$nse, rfVal$nse, dlVal$nse)
    max(best, na.rm = TRUE)
    
    input = which(best==max(best, na.rm = TRUE))
    if(length(input) ==0){next}
    comb = list(qVal, rfVal, dlVal)
    out = as.data.frame(comb[input])
    gage_stats$mode[i] = input
    
    rrmse = out$rrmse 
    r = out$r
    nse = out$nse
    kge = out$kge
    nrmse = out$nrmse
    rbias = out$rbias
    
    bestFit = c("QNT", "RF", "KNN")
    
    legend1 = c(bestFit[input], paste0("RRMSE=",signif(rrmse, 3)),paste0("NSE=",signif(nse, 3)),
                paste0("rBias=",signif(rbias, 3)), paste0("KGE=",signif(kge, 3)), paste0("NRMSE=",signif(nrmse, 3)))
    legend("topleft", legend1, xpd = TRUE, bty = "n", adj = 1, inset = c(-.05, -.1))
    
    
    gage_stats$RRMSE[i] = rrmse
    gage_stats$NSE[i] = nse
    gage_stats$KGE[i] = kge
    gage_stats$NRMSE[i] = nrmse
    gage_stats$rBias[i] = rbias
    gage_stats$Site_number[i] = paired_df$ID[1]
    gage_stats$n_Landsat_obs[i] = nrow(paired_df)
    print(i)
    
  } else{next}
}
################################################################################################################################
##Thailand gauges. 
################################################################################################################################
########################################################################################################################
thailandQ = "E:\\research\\RatingCurveAnalysis\\GaugeLocations\\DischargeDatasets\\Thailand\\"
gsim_all = read.csv("E:\\research\\GSIM\\GSIM_metadata\\GSIM_catalog\\GSIM_metadata.csv")
sites = as.data.frame(Site_number_xsections)
sites$gsim = gsub("_gsim", "", sites$Site_number_xsections)
sites$files= gsim_all$reference.no[match(sites$gsim, gsim_all$gsim.no)]



start = as.Date("1979-01-01")
data_val = Eff_widths
RC_End = as.Date("2014-12-31") ###WAs 2014
RC_year = format(RC_End, "%Y")
RC_year_1 = format(RC_End+1, "%Y")
data_val$ID = data_val$ID
data_val$ID_2 = data_val$ID
data_val$calc_mean = data_val$Effective_width

#data_val$Date = as.character(data_val$Date)
tab$ID = tab$id
tab$ID_2 = tab$id
tab$width_m = data_val$width_m[match(tab$ID, data_val$ID)]
#tab$change= tab$median
#tab$change[mapply(is.na, tab$change)] <- 0
tab$width_m = data$width_m[match(tab$ID_2, data$ID_2)]
xSecq=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecw=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecIDcol=grep("V", names(Site_number_xsections))
mInd = array(5, dimnames = NULL)
rangedf_1 = as.data.frame(matrix(numeric(), nrow = 1, ncol = 4))
gage_stats = as.data.frame(matrix(numeric(), nrow =length(Site_number_xsections), ncol = 23))
gage_stats_GRADES = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 14))
colnames(gage_stats_GRADES)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "p_val","Bias", "RRMSE", "avg_std", "change", 'RRMSE_median', "std_Q", "STDE", "mode")
as.data.frame(gage_stats_GRADES)
l_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
Mean_grades = as.vector(nrow(Site_number_xsections))
u_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals_1 = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
width_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
gage_quants_q = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
gage_quants_w = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
colnames(gage_stats)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "mode","Bias", "RRMSE","avg_std", "change", "RRMSE_median", "std_Q","STDE", "KGE", "NSE", "rBias",
                        "SDRR", "MRR", "NRMSE", "Q_50", "W_50")
as.data.frame(gage_stats)
gage_stats_col1 = as.vector(1)
gage_stats_col2 = as.vector(1)
gage_stats_GRADES_col1 = as.vector(1)
gage_stats_GRADES_col2 = as.vector(1)
paired_df_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
rmse = 1
width_grouping = 30
percentiles = c(0.05, 0.95)
training = 0.7
data_val = as.data.table(data_val)
setkey(data_val, ID)

library(profvis)
library(caret)
library(dataRetrieval)
#profvis({
for(i in 1:length(Site_number_xsections)){
  id = Site_number_xsections[i]
  paired_df = data_val[.(Site_number_xsections[i])]
  t = list.files(thailandQ, sites$files[1], full.names = TRUE)
  usgs_q = try(rbindlist(lapply(t, fread)))
  if(!is.error(usgs_q)){
    usgs_q$q = as.numeric(usgs_q$V2)
    usgs_q$Date = as.Date(usgs_q$V1, format = "%Y/%m/%d")
    usgs_q = usgs_q[usgs_q$q>=0,]
    paired_df = inner_join(paired_df, usgs_q)
    paired_df$Q = paired_df$q
    
    paired_df = inner_join(paired_df, usgs_q)
    #paired_df = distinct(paired_df, Date, .keep_all = TRUE)
    #########################################################################################################    
    paired_df$Q = paired_df$q  
    paired_df = paired_df[!is.na(paired_df$Q),]
    
    if(nrow(paired_df)<3){next}
    set.seed(1)
    trainIndex = createDataPartition(paired_df$Q, p = training,
                                     list = FALSE)
    
    Train = paired_df[ trainIndex,]
    Valid = paired_df[-trainIndex,]
    # years = unique(year(paired_df$Date))
    # train = round(length(years)*.6)
    # paired_df$Year = year(paired_df$Date)
    # set.seed(5)
    # t = sample(years, train)
    # TrainingSet = paired_df$Year%in%t
    # Train = paired_df[TrainingSet,]
    # Valid = paired_df[!TrainingSet,]
    
    
    
    
    
    # Train = paired_df[paired_df$Date>=as.Date(as.character('2015-01-01', format = "%Y-%m-%d")),]
    # Valid = paired_df[paired_df$Date<as.Date(as.character('2015-01-01', format = "%Y-%m-%d")),]
    
    if(nrow(Valid)<3|nrow(Train)<3){next}
    
    y = quantile(Train$Q, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    x = quantile(Train$calc_mean, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    if(length(unique(x))>1){
      ###Create function to estimate Q from a Landsat width based on quantile pair up. 
      spl = approxfun(x, y) #, method = "hyman")) #####either usse 'approxfun' or use splinefun with method = "hyman"
    } else{next}
    
    
    # par(pty = "s")
    # plot(c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)),c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)), type = "n", xlab = "In situ Discharge (cms)", ylab = "Landsat Discharge (cms)", main = "In situ vs Landsat")
    # points(Valid$Q, spl(Valid$calc_mean), col = "blue")
    # abline(0,1)
    
    plot(paired_df$Date, paired_df$q, type = "l")
    points(Valid$Date, spl(Valid$calc_mean), col = "green")
    
    names = c("rrmse", "nse", "kge", "nrmse", "rbias")
    
    ##Quantile validation 
    Valid$model = spl(Valid$calc_mean)
    qVal = validation(Valid$model, Valid$Q)
    qVal = as.data.frame(t(qVal))
    colnames(qVal) = names
    
    
    rf = randomForest(Q ~calc_mean+Date,data = Train, ntree= nrow(Train), mtry = 1)
    out = predict(rf, Valid)
    Valid$rf = out
    rfVal = validation(Valid$rf, Valid$Q)
    
    rfVal = as.data.frame(t(rfVal))
    colnames(rfVal) = names
    
    #points(Valid$Q, Valid$rf, pch = 19)
    points(Valid$Date, Valid$rf, col = "blue")
    
    plsFit = try(train(Q ~calc_mean+Date, 
                       data = Train,
                       method = "knn",
                       preProc=c("center", "scale"),
                       trControl = trainControl(method = "cv", number = c(5, 10, 20, 50, 100),search = "random")))
    
    preds <- predict(plsFit, newdata = Valid)
    dlVal = validation(preds, Valid$Q)
    dlVal = as.data.frame(t(dlVal))
    colnames(dlVal) = names
    
    #points(Valid$Q, preds, col = "red")
    points(Valid$Date, preds, col = "red")
    best = c(qVal$nse, rfVal$nse, dlVal$nse)
    max(best, na.rm = TRUE)
    
    input = which(best==max(best, na.rm = TRUE))
    if(length(input) ==0){next}
    comb = list(qVal, rfVal, dlVal)
    out = as.data.frame(comb[input])
    gage_stats$mode[i] = input
    
    rrmse = out$rrmse 
    r = out$r
    nse = out$nse
    kge = out$kge
    nrmse = out$nrmse
    rbias = out$rbias
    
    bestFit = c("QNT", "RF", "KNN")
    
    legend1 = c(bestFit[input], paste0("RRMSE=",signif(rrmse, 3)),paste0("NSE=",signif(nse, 3)),
                paste0("rBias=",signif(rbias, 3)), paste0("KGE=",signif(kge, 3)), paste0("NRMSE=",signif(nrmse, 3)))
    legend("topleft", legend1, xpd = TRUE, bty = "n", adj = 1, inset = c(-.05, -.1))
    
    
    gage_stats$RRMSE[i] = rrmse
    gage_stats$NSE[i] = nse
    gage_stats$KGE[i] = kge
    gage_stats$NRMSE[i] = nrmse
    gage_stats$rBias[i] = rbias
    gage_stats$Site_number[i] = paired_df$ID[1]
    gage_stats$n_Landsat_obs[i] = nrow(paired_df)
    gage_stats$validation[i] = nrow(Valid)
    gage_stats$training[i] = nrow(Train)
    print(i)
    
  } else{next}
}
################################################################################################################################
##Brazil gauges. 
################################################################################################################################
########################################################################################################################
br_files = "E:\\research\\RatingCurveAnalysis\\GaugeLocations\\DischargeDatasets\\Brazil\\"
gsim = read.csv("E:\\research\\RatingCurveAnalysis\\GaugeLocations\\30_m_1984\\GSIM_GRWL_30m.csv")
gsim_all = read.csv("E:\\research\\GSIM\\GSIM_metadata\\GSIM_catalog\\GSIM_metadata.csv")
gsim_br = gsim[grep("BR_", gsim$Station_Num), ]
gsim_br$gage = gsim_all$reference.no[match(gsim_br$Station_Num, gsim_all$gsim.no)]
sites = as.data.frame(Site_number_xsections)
sites$files= gsim_br$gage[match(sites[,1], paste0(gsim_br$Station_Num, "_gsim"))]
  
stations = paste0(sites$files, ".csv")
br_files = paste0(br_files, stations)







start = as.Date("1979-01-01")
data_val = Eff_widths
RC_End = as.Date("2014-12-31") ###WAs 2014
RC_year = format(RC_End, "%Y")
RC_year_1 = format(RC_End+1, "%Y")
data_val$ID = data_val$ID
data_val$ID_2 = data_val$ID
data_val$calc_mean = data_val$Effective_width

#data_val$Date = as.character(data_val$Date)
tab$ID = tab$id
tab$ID_2 = tab$id
tab$width_m = data_val$width_m[match(tab$ID, data_val$ID)]
#tab$change= tab$median
#tab$change[mapply(is.na, tab$change)] <- 0
tab$width_m = data$width_m[match(tab$ID_2, data$ID_2)]
xSecq=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecw=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecIDcol=grep("V", names(Site_number_xsections))
mInd = array(5, dimnames = NULL)
rangedf_1 = as.data.frame(matrix(numeric(), nrow = 1, ncol = 4))
gage_stats = as.data.frame(matrix(numeric(), nrow =length(Site_number_xsections), ncol = 23))
gage_stats_GRADES = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 14))
colnames(gage_stats_GRADES)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "p_val","Bias", "RRMSE", "avg_std", "change", 'RRMSE_median', "std_Q", "STDE", "mode")
as.data.frame(gage_stats_GRADES)
l_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
Mean_grades = as.vector(nrow(Site_number_xsections))
u_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals_1 = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
width_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
gage_quants_q = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
gage_quants_w = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
colnames(gage_stats)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "mode","Bias", "RRMSE","avg_std", "change", "RRMSE_median", "std_Q","STDE", "KGE", "NSE", "rBias",
                        "SDRR", "MRR", "NRMSE", "Q_50", "W_50")
as.data.frame(gage_stats)
gage_stats_col1 = as.vector(1)
gage_stats_col2 = as.vector(1)
gage_stats_GRADES_col1 = as.vector(1)
gage_stats_GRADES_col2 = as.vector(1)
paired_df_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
rmse = 1
width_grouping = 30
percentiles = c(0.05, 0.95)
training = 0.7
data_val = as.data.table(data_val)
setkey(data_val, ID)

library(profvis)
library(caret)
library(dataRetrieval)
#profvis({
for(i in 1:length(Site_number_xsections)){
  id = Site_number_xsections[i]
  paired_df = data_val[.(Site_number_xsections[i])]
  usgs_q = try(read.csv(br_files[i]))
  if(!is.error(usgs_q)){
    tab2 = melt(usgs_q)
    tab2 = tab2[grep("Vazao", tab2$variable),]
    tab2 = tab2[-grep("Status", tab2$variable),]
    substrRight <- function(x, n){
      substr(x, nchar(x)-n+1, nchar(x))
    }
    tab2$variable = as.character(tab2$variable)
    tab2$Day = substrRight(tab2$variable, 2)
    tab2$Day = as.numeric(tab2$Day)
    tab2$Date = as.Date(as.character(tab2$Data), format = "%d/%m/%Y")
    tab2$updatedDate = tab2$Date+(tab2$Day-1)
    tab2$Date = tab2$updatedDate
    usgs_q = tab2
    usgs_q$q = as.numeric(usgs_q$value)
    usgs_q$Date = usgs_q$Date
    usgs_q = usgs_q[usgs_q$q>=0,]
    paired_df = inner_join(paired_df, usgs_q)
    paired_df$Q = paired_df$q
    
    paired_df = inner_join(paired_df, usgs_q)
    #paired_df = distinct(paired_df, Date, .keep_all = TRUE)
    #########################################################################################################    
    paired_df$Q = paired_df$q  
    paired_df = paired_df[!is.na(paired_df$Q),]
    
    if(nrow(paired_df)<3){next}
    set.seed(1)
    trainIndex = createDataPartition(paired_df$Q, p = training,
                                     list = FALSE)
    
    Train = paired_df[ trainIndex,]
    Valid = paired_df[-trainIndex,]
    # years = unique(year(paired_df$Date))
    # train = round(length(years)*.6)
    # paired_df$Year = year(paired_df$Date)
    # set.seed(5)
    # t = sample(years, train)
    # TrainingSet = paired_df$Year%in%t
    # Train = paired_df[TrainingSet,]
    # Valid = paired_df[!TrainingSet,]
    
    
    
    
    
    # Train = paired_df[paired_df$Date>=as.Date(as.character('2015-01-01', format = "%Y-%m-%d")),]
    # Valid = paired_df[paired_df$Date<as.Date(as.character('2015-01-01', format = "%Y-%m-%d")),]
    
    if(nrow(Valid)<3|nrow(Train)<3){next}
    
    y = quantile(Train$Q, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    x = quantile(Train$calc_mean, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    if(length(unique(x))>1){
      ###Create function to estimate Q from a Landsat width based on quantile pair up. 
      spl = approxfun(x, y) #, method = "hyman")) #####either usse 'approxfun' or use splinefun with method = "hyman"
    } else{next}
    
    
    # par(pty = "s")
    # plot(c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)),c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)), type = "n", xlab = "In situ Discharge (cms)", ylab = "Landsat Discharge (cms)", main = "In situ vs Landsat")
    # points(Valid$Q, spl(Valid$calc_mean), col = "blue")
    # abline(0,1)
    
    plot(paired_df$Date, paired_df$q, type = "l")
    points(Valid$Date, spl(Valid$calc_mean), col = "green")
    
    names = c("rrmse", "nse", "kge", "nrmse", "rbias")
    
    ##Quantile validation 
    Valid$model = spl(Valid$calc_mean)
    qVal = validation(Valid$model, Valid$Q)
    qVal = as.data.frame(t(qVal))
    colnames(qVal) = names
    
    
    rf = randomForest(Q ~calc_mean+Date,data = Train, ntree= nrow(Train), mtry = 1)
    out = predict(rf, Valid)
    Valid$rf = out
    rfVal = validation(Valid$rf, Valid$Q)
    
    rfVal = as.data.frame(t(rfVal))
    colnames(rfVal) = names
    
    #points(Valid$Q, Valid$rf, pch = 19)
    points(Valid$Date, Valid$rf, col = "blue")
    
    plsFit = try(train(Q ~calc_mean+Date, 
                       data = Train,
                       method = "knn",
                       preProc=c("center", "scale"),
                       trControl = trainControl(method = "cv", number = c(5, 10, 20, 50, 100),search = "random")))
    
    preds <- predict(plsFit, newdata = Valid)
    dlVal = validation(preds, Valid$Q)
    dlVal = as.data.frame(t(dlVal))
    colnames(dlVal) = names
    
    #points(Valid$Q, preds, col = "red")
    points(Valid$Date, preds, col = "red")
    best = c(qVal$nse, rfVal$nse, dlVal$nse)
    max(best, na.rm = TRUE)
    
    input = which(best==max(best, na.rm = TRUE))
    if(length(input) ==0){next}
    comb = list(qVal, rfVal, dlVal)
    out = as.data.frame(comb[input])
    gage_stats$mode[i] = input
    
    rrmse = out$rrmse 
    r = out$r
    nse = out$nse
    kge = out$kge
    nrmse = out$nrmse
    rbias = out$rbias
    
    bestFit = c("QNT", "RF", "KNN")
    
    legend1 = c(bestFit[input], paste0("RRMSE=",signif(rrmse, 3)),paste0("NSE=",signif(nse, 3)),
                paste0("rBias=",signif(rbias, 3)), paste0("KGE=",signif(kge, 3)), paste0("NRMSE=",signif(nrmse, 3)))
    legend("topleft", legend1, xpd = TRUE, bty = "n", adj = 1, inset = c(-.05, -.1))
    
    
    gage_stats$RRMSE[i] = rrmse
    gage_stats$NSE[i] = nse
    gage_stats$KGE[i] = kge
    gage_stats$NRMSE[i] = nrmse
    gage_stats$rBias[i] = rbias
    gage_stats$Site_number[i] = paired_df$ID[1]
    gage_stats$n_Landsat_obs[i] = nrow(paired_df)
    gage_stats$validation[i] = nrow(Valid)
    gage_stats$training[i] = nrow(Train)
    print(i)
    
  } else{next}
}



gage_file = st_read("E:\\research\\RatingCurveAnalysis\\GaugeLocations\\30_m_1984\\combined.shp")
gage_stats$GRWL_width_m = gage_file$GRWL_wd[match(gage_stats$Site_number, gage_file$Sttn_Nm)]

gage_stats = as.data.frame(gage_stats)
apply(gage_stats, 2, median, na.rm= TRUE)
gage_stats = as.data.frame(gage_stats)
boxplot(gage_stats$RRMSE[gage_stats$GRWL_width_m>90], ylim = c(0,200))

gage_file = gage_file[gage_file$Sttn_Nm%in%gage_stats$Site_number,]
gage_stats$Sttn_Nm = gage_stats$Site_number
gage_stats_vals = merge(gage_file, gage_stats,by = 'Sttn_Nm')

library(ggplot2)
library(tmap)

tmap_mode("view")

tm_shape(gage_stats_vals)+
  tm_bubbles(col = "NSE",size = 0.1, breaks = c(-1,-.5, 0,.5, 1))
tm_shape(gage_stats_vals)+
  tm_bubbles(col = "KGE",size = 0.02, breaks = c(-1, -.4, 1))
tm_shape(gage_stats_vals)+
  tm_bubbles(col = "mode",size = 0.02)


tm_shape(gage_stats_vals)+
  tm_bubbles(col = "rBias",size = 0.02, breaks = c(-100, 0, 100))



tm_shape(gage_stats_vals)+
  tm_bubbles(col = "n_Landsat_obs",size = 0.02, breaks = c(0, 100, 500, 1000, 2000))



gage_stats_vals1 = gage_stats_vals%>%select(mode, geometry, rBias, NRMSE, RRMSE, KGE, NSE, n_Landsat_obs, Site_number, lakeFlg, ds2GRWL, GRWL_wd)
st_write(gage_stats_vals1, "E:\\research\\RatingCurveAnalysis\\stats\\grdc_prelim1.shp")










######See if GRDC record is same or different for same gauges. 
gage_stats$gsim = gsub("_gsim", "", gage_stats$Site_number)
gage_stats$GRDC = gsim_all$paired.db.no[match(gage_stats$gsim, gsim_all$gsim.no)]
gage_stats$GRDC = gsim_all$grdb.no[match(gage_stats$GRDC, gsim_all$reference.no)]

gage_stats$ANA = sites$files


##ANA
usgs_q = try(read.csv(br_files[6]))
  tab2 = melt(usgs_q)
  tab2 = tab2[grep("Vazao", tab2$variable),]
  tab2 = tab2[-grep("Status", tab2$variable),]
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  tab2$variable = as.character(tab2$variable)
  tab2$Day = substrRight(tab2$variable, 2)
  tab2$Day = as.numeric(tab2$Day)
  tab2$Date = as.Date(as.character(tab2$Data), format = "%d/%m/%Y")
  tab2$updatedDate = tab2$Date+(tab2$Day-1)
  tab2$Date = tab2$updatedDate
  usgs_q = tab2
  usgs_q$q = as.numeric(usgs_q$value)
  usgs_q$Date = usgs_q$Date
  usgs_q = usgs_q[usgs_q$q>=0,]
  ana = usgs_q
ana = ana[order(ana$Date),]
  
##GRDC
  grdc_path = "E:\\research\\RatingCurveAnalysis\\GaugeLocations\\DischargeDatasets\\GRDC\\"
  grdc_files = paste0("3618951", "_Q_Day.Cmd.txt")
  grdc_files = paste0(grdc_path, grdc_files)
  usgs_q = try(read.table(grdc_files, stringsAsFactors = FALSE))
    usgs_q$V1 = substr(usgs_q$V1, 0, 10)
    usgs_q$datetime = usgs_q$V1
    usgs_q$q = as.numeric(usgs_q$V2)
    usgs_q$Date = as.Date(usgs_q$datetime, format = "%Y-%m-%d")
    usgs_q = usgs_q[usgs_q$q>=0,]



plot(usgs_q$Date, usgs_q$q, type = "l")
lines(ana$Date, ana$q, col = "red")




join = inner_join(ana, usgs_q, by = "Date")
plot(join$q.x, join$q.y)


