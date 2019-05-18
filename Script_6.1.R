# set this directory as the root
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# load("E:/OneDrive - Cal Poly Pomona/Research/Air Data Cities - From 20180919/.RData")

# libraries
library(lubridate)
library(pracma)
library(zoo)
library(imputeTS)
library(tidyverse)
library(foreach)
library(GGally)
library(imputeTS)
library(REdaS)
library(keras)
library(recipes)
library(reshape2)
library(Metrics)
library(hydroGOF)
library(mefa)
library(exactRankTests)


# misc and basic functions
{
  # remove all variables
  rm(list = ls())
  
  # pasting
  Texcat <- function(..., s = "") {
    return(paste(..., sep = s))
  }
  
  # theming
  BasicThemeNoLegend <- theme_bw() + theme(legend.position = "none", axis.text = element_text(size = rel(0.8), margin = margin(r = 6)))
  BasicThemeWithLegend <- theme_bw() + theme(legend.position = "bottom",axis.text = element_text(size = rel(0.8), margin = margin(r = 5)), legend.text = element_text(margin = margin(r = 0.15, unit = "in")))
  
  # function for foreach
  cfun <- function(...) NULL
  
  # Linear variable name
  paramNameToEvalLinear <- "Texcat(dataExtension, '_', city, '_', pol, '_', lag, '_', current, '_', paste(param, collapse = '_'), '_', nTest)"
  
  # LSTM variable name
  paramNameToEvalLSTM <- "Texcat(dataExtension, '_',city, '_', pol, '_', lag, '_', current, '_', paste(param, collapse = '_'), '_', nTest, '_', batch, '_', step, '_',
  epoch, '_', numUnit, '_', stop, '_', split)"
  
  # get parameter text (for saving files)
  GetParamText <- function(model) {
    if (model == "s" | model == "c") {
      return(Texcat(model, ".", eval(parse(text = paramNameToEvalLinear))))
    } else {
      return(Texcat(model, ".", eval(parse(text = paramNameToEvalLSTM))))
    }
  }
  
  # set regression query
  SetQuery <- function(city, pol, lag, current, param, nTest, batch, step, epoch, numUnit, stop, split) {
    assign("city", city, envir = globalenv())
    assign("pol", pol, envir = globalenv())
    assign("lag", lag, envir = globalenv())
    assign("current", current, envir = globalenv())
    assign("param", param, envir = globalenv())
    assign("nTest", nTest, envir = globalenv())
    assign("batch", batch, envir = globalenv())
    assign("step", step, envir = globalenv())
    assign("epoch", epoch, envir = globalenv())
    assign("numUnit", numUnit, envir = globalenv())
    assign("stop", stop, envir = globalenv())
    assign("split", split, envir = globalenv())
  }
  
  # get regression query
  GetQuery <- function() {
    print(c(city, pol, lag, current, param, nTest, batch, step, epoch, numUnit, stop, split))
  }
  
  # get MAPE
  GetMAPE <- function(actual, pred) {
    return(
      sprintf(
        "%2.2f",
        (100 * (mean(abs(pred - actual) / actual, na.rm = TRUE)))
      )
    )
  }
}

# constants
{
  cityNames <- c(
    "Beijing",
    "Baoding",
    "Cangzhou",
    "Chengde",
    "Handan",
    "Hengshui",
    "Langfang",
    "Qinhuangdao",
    "Shijiazhuang",
    "Tangshan",
    "Tianjin",
    "Xingtai",
    "Zhangjiakou"
  )
  
  pollutantNames <- c(
    "CO",
    "NO2",
    "SO2",
    "O3",
    "PM2.5",
    "PM10"
  )
  
  
  UnitsPol <- c(
    expression(CO ~ (mg/m^3)),
    expression(NO2 ~ (ug/m^3)),
    expression(SO2 ~ (ug/m^3)),
    expression(O3 ~ (ug/m^3)),
    expression(PM2.5 ~ (ug/m^3)),
    expression(PM10 ~ (ug/m^3))
  )
  
  
  UnitsEnv <- c(
    " (hPa)",
    " (%)",
    " (C)"
  )
  
  # minimum value for log transformation
  minValue <- 0.001
  
  printPlot <- 0
  dataExtension <- 1
}

MakeFirstData <- function() {
  # read files
  {
    foreach(str = cityNames, .combine = "cfun") %do% {
      name <- Texcat("Env_", str, ".csv")
      data <- read.csv(name)
      data[, 1] <- data[, 1] %>% as.POSIXct(format = "%Y-%m-%d %H:%M")
      assign(str, data, envir = globalenv())
    }
    
    foreach(str = pollutantNames, .combine = "cfun") %do% {
      name <- Texcat("Pol_", str, ".csv")
      data <- read.csv(name)
      colnames(data)[1] <- "Time"
      assign(str, data, envir = globalenv())
    }
    
    angleDistanceMatrix = read.csv("cities.csv", header = F) %>% .[-1, ] %>% .[, -1]
    assign("angleDistanceMatrix", angleDistanceMatrix, envir = globalenv())
  }
  
  # edit format
  {
    foreach(str = cityNames, .combine = "cfun") %do% {
      data <- get(str)
      d <- apply(data[, 2:3], 1, function(x) Texcat(2016, "/", x[1], "/", x[2])) %>% as.POSIXct(format = "%Y/%m/%d") %>% yday()
      int <- seq(as.POSIXct("2016/01/01 00:00:00"), as.POSIXct("2016/12/31 23:59:00"), by = "1 hour")
      h <- (d - 1) * 24 + data[, 4] + 1
      data[, 10] <- h
      newData <- data.frame(matrix(NA, 8784, 6))
      newData[, 1] <- int
      newData[data[, 10], 2:6] <- data.frame(data[, 5:9])
      newData[, 6] <- newData[, 6] %>% deg2rad()
      newData[, 7] <- newData[, 2] * cos(newData[, 6])
      newData[, 8] <- newData[, 2] * sin(newData[, 6])
      dr <- newData[, 6]
      newData <- newData[, c(1, 3:5, 7:8)]
      for (i in 1:length(pollutantNames)) {
        newData[, 6 + i] <- get(pollutantNames[i])[, str]
      }
      colnames(newData) <- c("Time", "AP", "HM", "TP", "WX", "WY", pollutantNames)
      newData[, 2:6] <- apply(newData[, 2:6], 2, function(x) {
        baseRange <- boxplot(x, plot = F) %>% .$stats
        lower <- baseRange[2] - (baseRange[4] - baseRange[2]) * 4.5
        higher <- baseRange[4] + (baseRange[4] - baseRange[2]) * 4.5
        x[x < lower | x > higher] <- NA
        return(x)
      })
      newData[, -1:-6] <- apply(newData[, -1:-6], 2, function(x) {
        x[x < 0] <- 0
        return(x)
      })
      newData[, -1] <- na.interpolation(newData[, -1], option = "linear")
      newData[newData[, "HM"] > 100, "HM"] <- 100
      newData[, "WI"] <- (newData[, "WX"]^2 + newData[, "WY"]^2) %>% sqrt()
      zeroIndex <- newData[, "WI"] == 0
      newData[, "DR"] <- ifelse(zeroIndex, dr, atan2(newData[, "WY"], newData[, "WX"])) %>% na.interpolation(option = "linear") %>% rad2deg()
      newData[, "DR"] <- newData[, "DR"] %% 360
      newData <- newData[, c(1:4, 13:14, 7:12)]
      assign(str, newData, envir = globalenv())
    }
  }
  
  # time
  t <- Beijing[, c("Time")]
  assign("t", t, envir = globalenv())
  
  if (dataExtension == 1) {
    print("Data eXtension is enabled.")
    extension <- read.csv("extension.csv")
    cityIDs <- c(
      1816670, 1816971, 1816080, 2038087, 1808963,
      1808392, 1804540, 1797595, 1795270, 1793346,
      1792947, 1788927, 2033196
    )
    
    for (i in 1:length(cityIDs)) {
      chunk <- extension[extension[, 3] == cityIDs[i], c(2, 7, 10, 13, 14, 15)]
      colnames(chunk) <- c("Time", "AP", "HM", "TP", "WI", "DR")
      chunk$Time <- gsub(" \\+0000 UTC", "", chunk$Time) %>% as.character %>% as.POSIXlt(format = "%Y-%m-%d %H:%M:%OS")
      chunk$Time <- chunk$Time + hours(8)
      chunk = chunk[!duplicated(chunk),]
      
      newData <- data.frame(matrix(NA, 8784, 7))
      colnames(newData) <- c("AP", "HM", "TP", "WX", "WY", "WI", "DR")
      print(dim(chunk)[1])
      newData[, "TP"] <- chunk$AP
      newData$TP <- newData$TP - 273
      newData[, "AP"] <- chunk$HM
      newData[, "HM"] <- chunk$TP
      newData[, "WI"] <- chunk$WI
      newData[, "DR"] <- (-chunk$DR + 90) %>% deg2rad()
      newData[, "WX"] <- newData[, "WI"] * cos(newData[, "DR"])
      newData[, "WY"] <- newData[, "WI"] * sin(newData[, "DR"])
      newData <- na.interpolation(newData, option = "linear")
      newData[, "WI"] <- (newData[, "WX"]^2 + newData[, "WY"]^2) %>% sqrt()
      newData[, "DR"] <- ifelse(zeroIndex, dr, atan2(newData[, "WY"], newData[, "WX"])) %>% na.interpolation(option = "linear") %>% rad2deg()
      newData[, "DR"] <- (newData[, "DR"] + 180) %% 360
      assign(cityNames[i], data.frame(t, newData[, -4:-5], get(cityNames[i])[, -1:-6]), envir = globalenv())
    }
  }
}

# basic plots
Printplots = function(){
  # air pollution
  foreach(city = cityNames, .combine = "cfun") %do% {
    i = 1
    foreach(str = pollutantNames, .combine = "cfun") %do% {
      df <- data.frame(t,get(city) %>% .[,str])
      colnames(df) = c("Time","val")
      assign(Texcat("graph",i),
             ggplot(df) + geom_line(aes(x = Time, y = val, col = "blue")) +
               geom_line(aes(x = Time, y = rollmean(val, 24, na.pad = T), col = "green")) +
               geom_line(aes(x = Time, y = rollmean(val, 168, na.pad = T), col = "red")) +
               BasicThemeNoLegend + labs(y = str) +
               scale_color_manual(values=c( "#BB933E", "#255B87","#CC5287"),
                                  label = c(str, "24-hour average", "1-week average")) + 
               theme(legend.title= element_blank())
      )
      i = i +1
    }
    graph0 = cowplot::plot_grid(plotlist = mget(Texcat("graph",1:6)), ncol = 2)
    legend <- cowplot::get_legend(graph1 + BasicThemeWithLegend + 
                                    theme(legend.title = element_blank(),legend.text = element_text(size = rel(0.95))) +
                                    guides(shape = guide_legend(nrow = 1), col = guide_legend(nrow = 1)) +
                                    scale_color_manual(values=c( "#BB933E", "#255B87","#CC5287"),
                                                       label = c("Concentration", "24-hour average", "1-week average")))
    graph = cowplot::plot_grid(graph0, legend, ncol = 1, rel_heights = c(36, 1))
    ggsave(Texcat("./basic_plots/",city, ".png"),graph , width = 9, height = 7.2, dpi = 200, units = "in")
  }
  
  # air pollution hist
  foreach(city = cityNames, .combine = "cfun") %do% {
    i = 1
    foreach(str = pollutantNames, .combine = "cfun") %do% {
      df <- get(city) %>% .[,str] %>% data.frame
      assign(Texcat("graph",i),
             ggplot(df) + geom_histogram(aes(x = .),fill = "#BB933E") +
               BasicThemeNoLegend + labs(y = UnitsPol[i]) +
               theme(legend.title= element_blank())
      )
      s = Texcat(s,"graph",i,", ")
      i = i +1
    }
    graph0 = cowplot::plot_grid(plotlist = mget(Texcat("graph",1:6)), ncol = 2)
    ggsave(Texcat("./basic_plots/",city, "_hist.png"),graph0 , width = 9, height = 7.2, dpi = 200, units = "in")
  }
  
  # meteorological
  foreach(city = cityNames, .combine = "cfun") %do% {
    i = 1
    foreach(str = c("AP","HM","TP"), .combine = "cfun") %do% {
        df <- data.frame(t,get(city) %>% .[,str])
        colnames(df) = c("Time","val")
        assign(Texcat("graph",i),
               ggplot(df) + geom_line(aes(x = Time, y = val, col = "blue")) +
                 geom_line(aes(x = Time, y = rollmean(val, 24, na.pad = T), col = "green")) +
                 geom_line(aes(x = Time, y = rollmean(val, 168, na.pad = T), col = "red")) +
                 BasicThemeNoLegend + labs(y = str) +
                 scale_color_manual(values=c( "#BB933E", "#255B87","#CC5287"),
                                    label = c(str, "24-hour average", "1-week average")) + 
                 theme(legend.title= element_blank())
        )
        i = i +1
    }
    graph0 = cowplot::plot_grid(plotlist = mget(Texcat("graph",1:3)), ncol = 2)
    legend <- cowplot::get_legend(graph1 + BasicThemeWithLegend + 
                                    theme(legend.title = element_blank(),legend.text = element_text(size = rel(0.95))) +
                                    guides(shape = guide_legend(nrow = 1), col = guide_legend(nrow = 1)) +
                                    scale_color_manual(values=c( "#BB933E", "#255B87","#CC5287"),
                                                       label = c("Concentration", "24-hour average", "1-week average")))
    graph = cowplot::plot_grid(graph0, legend, ncol = 1, rel_heights = c(24, 1))
    ggsave(Texcat("./basic_plots/",city, "_met.png"),graph , width = 9, height = 4.8, dpi = 200, units = "in")
  }
  
  # wind rose
  source("Windrose.R")
  foreach(city = cityNames, .combine = "cfun") %do% {
      df <- get(city) %>% .[,c("WI","DR")]
      graph = plot.windrose(spd = df$WI, dir =  (-df$DR+90)%%360,
                    spdres = 2,
                    dirres = 22.5,
                    spdmin = 0,
                    spdmax = max(df$WI),
                    palette = "Spectral") + theme_bw() + theme(axis.title.x = element_blank()) +  labs(y = "Count (Hours)")
      ggsave(Texcat("./basic_plots/",city,"_WR", ".png"), width = 8, height = 4.5, dpi = 200, units = "in")
  }
  
  # corr in beijing
  {
    source("PointsWithLine.R")
    df <- data.frame(get(cityNames[1]) %>% .[-1:-6])
    graph <- ggpairs(df,
                     lower = list(continuous = PointsWithLine),
                     diag = list(continuous = PlotDiag),
                     upper = list(continuous = wrap("cor", alignPercent = 0.6))
    ) + BasicThemeNoLegend + theme(strip.text = element_text(size = rel(0.95)))
    ggsave(Texcat("./corr/",cityNames[1],"_self.png"), graph, width = 9, height = 6, dpi = 200, units = "in")
  }
  
  # corr in beijing vs. meteorological
  {
    source("PointsWithCorr.R")
    df <- data.frame(get(cityNames[1]) %>% .[,c(2:4,7:12)])
    graph <- ggduo(df,
                   columnsX = 4:9,
                   columnsY = 1:3,
                   types = list(continuous = PointsWithCorr)
    ) + BasicThemeNoLegend + theme(strip.text = element_text(size = rel(0.95)))
    ggsave(Texcat("./corr/",cityNames[1],"_self_met.png"), graph, width = 9, height = 6, dpi = 200, units = "in")
  }
  
  # corr vs. beijing
  source("PointsWithCorr.R")
  foreach(city = cityNames[-1], .combine = "cfun") %do% {
    
    df <- data.frame(get(cityNames[1]) %>% .[-1:-6],get(city) %>% .[-1:-6])
    colnames(df) = rep(pollutantNames,2)
    colnames(df)[-1:-6] = colnames(df)[-1:-6] %>% Texcat("A")
    graph <- ggduo(df,
                   columnsX = 1:6,
                   columnsY = 7:12,
                   types = list(continuous = PointsWithCorr)
    ) + BasicThemeNoLegend + theme(strip.text = element_text(size = rel(0.95)))
    ggsave(Texcat("./corr/",cityNames[1],"_",city, ".png"), graph, width = 9, height = 6, dpi = 200, units = "in")
  }
  
  # corr vs. beijing all
  {
  df = data.frame(get(cityNames[1]) %>% .[-1:-6],get(city) %>% .[-1:-6])
  df = rep(df,12)
  df[,13] = rep(cityNames[-1],each=length(t)[1])
  foreach(city = cityNames[-1], .combine = "cfun") %do% {
    df[df[,13] == city,7:12] <- data.frame(get(city) %>% .[-1:-6])
  }
  df = df[,-13]
  
  colnames(df) = rep(pollutantNames,2)
  colnames(df)[-1:-6] = colnames(df)[-1:-6] %>% Texcat("A")
  graph <- ggduo(df,
                 columnsX = 1:6,
                 columnsY = 7:12, 
                 types = list(continuous = PointsWithCorr)
  ) +
    BasicThemeNoLegend + theme(strip.text = element_text(size = rel(0.95)))
  ggsave(Texcat("./corr/",cityNames[1],"_duo.png"), graph, width = 9, height = 6, dpi = 200, units = "in")
  }
  
  # corr vs. beijing all in one plot
  {
  df = data.frame(0,0,0,0)
  colnames(df) = c("City","Pollutant","Correlation","Distance")
  foreach(str = pollutantNames, .combine = "cfun") %do% {
    chunk <- get(str)[,-1]
    corr <- chunk[,-1] %>% apply(2,function(x) cor(x,chunk[,1]))
    chunk = data.frame(names(corr),rep(str,length(corr)),corr, angleDistanceMatrix[-1,1])
    colnames(chunk) = c("City","Pollutant","Correlation","Distance")
    df = rbind(df,chunk)
  }
  graph <- ggplot(df[-1,]) + geom_point(aes(x= City %>% fct_rev,y=Correlation,color=Pollutant,shape=Pollutant), size = 2, stroke = 1) + coord_flip() +
    BasicThemeWithLegend + guides(color = guide_legend(nrow = 1),shape = guide_legend(nrow = 1)) +
    theme(legend.text = element_text(margin = margin(l = -0.06, r = 0.15, unit = "in")))+
    scale_y_continuous(limits = c(0, 1)) +
    scale_shape_manual(values = 1:6) + labs(x = "City", y = "Correlation vs. Beijing's value")
  ggsave(Texcat("./corr/",cityNames[1],"_all.png"), graph, width = 9, height = 6, dpi = 200, units = "in")
  
  graph <- ggplot(df[-1,]) + geom_point(aes(x= Distance %>% as.numeric() ,y=Correlation,color=Pollutant,shape=Pollutant), size = 2, stroke = 1) + coord_flip() +
    geom_smooth(aes(x= Distance %>% as.numeric() ,y=Correlation,color=Pollutant),method = lm, se = FALSE, size = 0.5) +
    BasicThemeWithLegend + guides(color = guide_legend(nrow = 1),shape = guide_legend(nrow = 1)) +
    theme(legend.text = element_text(margin = margin(l = -0.06, r = 0.15, unit = "in")))+
    scale_y_continuous(limits = c(0, 1)) + 
    scale_shape_manual(values = 1:6) + labs(x = "Distance from Beijing (km)", y = "Correlation vs. Beijing's value")
  ggsave(Texcat("./corr/",cityNames[1],"_alla.png"), graph, width = 9, height = 6, dpi = 200, units = "in")
  }
  
  # corr vs. beijing wind
  PlotCorrWind = function(aData){
    source("PointsWithCorr.R")
    df <- data.frame(get(cityNames[1]) %>% .[-1:-6],aData[,(dim(aData)[2]-12):dim(aData)[2]])
    colnames(df)[1:6] = pollutantNames
    colnames(df)[-1:-6] = cityNames
    graph <- ggduo(df,
                   columnsX = 1:6,
                   columnsY = 7:19,
                   types = list(continuous = PointsWithCorr)
    ) +
      BasicThemeNoLegend + theme(strip.text = element_text(size = rel(0.95)),axis.title = element_text(size = rel(0.8)))
    ggsave(Texcat("./corr/",cityNames[1],"_wind.png"), graph, width = 9, height = 12.1, dpi = 200, units = "in")
  }
  
  # contribution of wind
  PlotContWind = function(aData){
    for (i in 1:12) {
      df = data.frame(t,aData[,Texcat("CW.",i+1)])
      colnames(df) = c("Time","val")
      assign(Texcat("graph",i),
             ggplot(df) + geom_line(aes(x = Time, y = val, col = "blue")) +
               geom_line(aes(x = Time, y = rollmean(val, 24, na.pad = T), col = "green")) +
               geom_line(aes(x = Time, y = rollmean(val, 168, na.pad = T), col = "red")) +
               BasicThemeNoLegend + labs(y = "CW (m/s)") + ggtitle(label = cityNames[i+1]) +
               scale_color_manual(values=c( "#BB933E", "#255B87","#CC5287"),
                                  label = c(str, "24-hour average", "1-week average")) + ylim(0,10) +
               theme(legend.title= element_blank(), plot.title = element_text(size = rel(0.8)))
      )
      i = i +1
    }
    graph0 = cowplot::plot_grid(plotlist = mget(Texcat("graph",1:12)), ncol = 2)
    legend <- cowplot::get_legend(graph1 + BasicThemeWithLegend + 
                                    theme(legend.title = element_blank(),legend.text = element_text(size = rel(0.95))) +
                                    guides(shape = guide_legend(nrow = 1), col = guide_legend(nrow = 1)) +
                                    scale_color_manual(values=c( "#BB933E", "#255B87","#CC5287"),
                                                       label = c("Concentration", "24-hour average", "1-week average")))
    graph = cowplot::plot_grid(graph0, legend, ncol = 1, rel_heights = c(72, 1))
    ggsave(Texcat("./basic_plots/CW.png"),graph , width = 9, height = 12.5, dpi = 200, units = "in")
  }
  
  # acf in beijing
  {
  foreach(str = pollutantNames, .combine = "cfun") %do% {
    i <- which(pollutantNames == str)
    env <- get(str)
    df <- data.frame(-72:72)
    df[, 2] <- ccf(env[, 2], env[, 2], plot = FALSE, lag.max = 72)$acf
    colnames(df) <- c("Lag", "value")
    df <- as.data.frame.array(reshape2::melt(df, id = "Lag"))
    graph <- ggplot(df, aes(x = Lag, y = value)) +
      geom_hline(yintercept = 0) + geom_line(aes(colour = variable)) +
      geom_point(aes(colour = variable, shape = variable), stroke = 0.7) +
      scale_shape_manual(values = 1:4) + BasicThemeNoLegend +
      scale_color_manual(values=c( "#BB933E", "#255B87","#CC5287", "#93CC52")) + 
      labs(x = "Lag (Hour)", y = Texcat("ACF (vs. ", str, ")"))
    assign(Texcat("graph", i), graph)
  }
  graphGrid <- cowplot::plot_grid(graph1, graph2, graph3, graph4, graph5, graph6, align = "v", ncol = 2)
  ggsave("./corr/acf_Beijing.png", width = 9, height = 5.8, dpi = 200, units = "in")
  }
  
  # ccf vs beijing
  for (j in c(2,6,10)) {
    foreach(str = pollutantNames, .combine = "cfun") %do% {
      i <- which(pollutantNames == str)
      env <- get(str)
      df <- data.frame(-72:72)
      for (k in j:(j+3)) {
        df[, k-j+2] <- ccf(env[, 2], env[, 1 + k], plot = FALSE, lag.max = 72)$acf
      }
      colnames(df) <- c("Lag", cityNames[j:(j+3)])
      df <- as.data.frame.array(reshape2::melt(df, id = "Lag"))
      graph <- ggplot(df, aes(x = Lag, y = value)) +
        geom_hline(yintercept = 0) + geom_line(aes(colour = variable)) +
        geom_point(aes(colour = variable, shape = variable), stroke = 0.7) +
        scale_shape_manual(values = 1:4) + BasicThemeNoLegend +
        scale_color_manual(values=c( "#BB933E", "#255B87","#CC5287", "#93CC52")) + 
        labs(x = "Lag (Hour)", y = Texcat("CCF (vs. ", str, ")"))
      assign(Texcat("graph", i), graph)
    }
    graphGrid <- cowplot::plot_grid(graph1, graph2, graph3, graph4, graph5, graph6, align = "v", ncol = 2)
    legend <- cowplot::get_legend(graph1 + theme(
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.95), margin = margin(r = 6)),
      legend.position = "bottom"
    )
    + guides(shape = guide_legend(nrow = 1), col = guide_legend(nrow = 1)))
    graph <- cowplot::plot_grid(graphGrid, legend, ncol = 1, rel_heights = c(30, 1))
    ggsave(Texcat("./corr/ccf_Beijing_",cityNames[j],".png"), width = 9, height = 5.8, dpi = 200, units = "in")
  }
}

# Create data for analysis
CreateData <- function() {
  data <- t
  
  for (i in 1:(cityNames %>% length())) {
    data <- data.frame(data, get(cityNames[i])[, -1])
    colnames(data)[(i * 11 - 9):(i * 11 + 1)] <- Texcat(c("AP", "HM", "TP", "WI", "DR", pollutantNames), ".", i)
  }
  data <- data[, -1]
  return(data)
}

### combine wind data
CombineWind <- function(data, city, index) {
  aData <- data
  for (i in 1:length(cityNames)) {
    if (city != cityNames[i]) {
      angle <- angleDistanceMatrix[min(index, i), max(index, i)] %>% as.character() %>% as.numeric() + ifelse(index < i, 180, 0)
      aData[, Texcat("CW.", i)] <- pmax(
        aData[, Texcat("WI.", i)] *
          cos((aData[, Texcat("DR.", i)] - angle) %>%
                deg2rad()),
        0
      )
    } else {
      aData[, Texcat("CW.", i)] <- aData[, Texcat("WI.", i)]
    }
    aData[, Texcat("CW.", i)] <- aData[, Texcat("CW.", i)] + minValue
  }
  truncateIndex <- colnames(aData) %in% DR | colnames(aData) %in% WI
  aData <- data.frame(aData[, !truncateIndex])
  return(aData)
}

### extract variables
ExtractVariables <- function(aData, param, index) {
  variables <- c("")
  foreach(p = param, .combine = "cfun") %do% {
    if (gsub(".*\\.", "", p) == "S") {
      variables <- c(variables, gsub("\\.S", "", p) %>% Texcat(".", index))
    } else if (gsub(".*\\.", "", p) == "A") {
      variables <- c(variables, gsub("\\.A", "", p) %>% get())
    }
  }
  variables <- variables[-1]
  eData <- aData[, variables]
  return(eData)
}

### get standardization parameter
GetStandardizeParameters <- function(data) {
  logIndex <- colnames(data)[gsub("\\.[0-9]{1,2}$", "", colnames(data)) %in% c("CW", pollutantNames)]
  sData <- data
  sData[, logIndex] <- sData[, logIndex] %>% log()
  scale <- sData %>% apply(2, sd)
  mean <- sData %>% apply(2, mean)
  return(list(scale, mean, logIndex))
}

### standardization
Standardize <- function(data, scale, mean, logIndex) {
  name <- intersect(colnames(data), names(scale))
  extractNames <- intersect(colnames(data), logIndex)
  sData <- data
  sData[, extractNames] <- sData[, extractNames] %>% log()
  for (i in 1:length(name)) {
    sData[, name[i]] <- (sData[, name[i]] - mean[name[i]]) / scale[name[i]]
  }
  return(sData)
}

### create lag data according to the wind
CreateLagData <- function(eData, index, lag, target) {
  ln <- dim(eData)[1]
  polColumns <- colnames(eData)[!(gsub("\\.[0-9]{1,2}$", "", colnames(eData)) %in% c("AP", "HM", "TP", "CW") |
                                    gsub(".*\\.([0-9]{1,2})$", "\\1", colnames(eData)) == index)]
  
  if (length(polColumns) == 0) {
    return("")
  }
  
  df <- data.frame(t)
  foreach(str = polColumns, .combine = "cfun") %do% {
    polIntensity <- eData[, str]
    targetIntensity <- eData[, target]
    
    interval <- 288
    nr <- ((ln - lag) %% interval + 1):(ln - lag)
    suggestedLag <- rep(NA, length(ln))
    for (i in seq((ln - lag) %% interval, ln - lag - interval, interval)) {
      range <- (i + 1):(i + interval)
      ccf <- ccf(targetIntensity[range], (polIntensity)[range], lag = interval / 2, plot = F) %>%
        .$acf %>%
        .[((interval + 1) / 2):(interval + 1)]
      suggestedLag[range] <- which(ccf == max(ccf))
    }
    suggestedLag <- suggestedLag %>% na.omit()
    suggestedLag <- loess(suggestedLag ~ nr, span = 0.50) %>% .$fitted
    suggestedLag <- c(rep(suggestedLag[1], nr[1] - 1), suggestedLag, rep(suggestedLag %>% tail(1), ln - nr %>% tail(1)))
    suggestedLag <- suggestedLag %>% pmax(0) %>% pmin(144) %>% floor()
    df[, str] <- suggestedLag
    
    newData <- rep(NA, ln + max(suggestedLag))
    newData[(1:ln) + suggestedLag] <- polIntensity
    eData[, str] <- newData[1:ln] %>% na.interpolation(option = "spline") %>% pmax(0) + minValue
    eData[(ln - lag + 1):ln, str] <- NA
  }
  colnames(df) <- c("Time", cityNames[-index])
  df <- df %>% reshape2::melt("Time")
  ggplot(df) + geom_line(aes(x = Time, y = value, color = variable)) +
    BasicThemeWithLegend + theme(strip.text = element_text(size = rel(0.95), margin = margin(r = 6))) +
    labs(y = Texcat("Suggested Lag (vs. Beijing's ", gsub("\\.[0-9]{1,2}$", "", target), ") / Hours")) +
    guides(color = guide_legend(title = "City", nrow = 2))
  ggsave(Texcat("./suggested_lag/", GetParamText("c"), ".png"), last_plot(), width = 10, height = 6, dpi = 200, units = "in")
  
  extractNames <- colnames(eData)[gsub("\\.[0-9]{1,2}$", "", colnames(eData)) %in%
                                    c(
                                      gsub("\\.[0-9]{1,2}$", "", polColumns) %>% unique(),
                                      "AP", "HM", "TP"
                                    ) | colnames(eData) == Texcat("CW.", index)]
  
  return(eData[, extractNames])
}


### Regression
DoRegressions <- function(city, pol, lag, current, param, nTest, batch, step, epoch, numUnit, stop, split, s = F, c = F, l = F) {
  data <- CreateData()
  index <- which(city == cityNames)
  target <- Texcat(pol, ".", index)
  
  aData <- CombineWind(data, city, index)
  PlotCorrWind(aData)
  eData <- ExtractVariables(aData, param, index)
  scaler <- GetStandardizeParameters(eData)
  
  if (s == T) {
    SimpleLinear(eData, scaler, index, target, lag, current, nTest)
  }
  if (length("CW.A" == param) > 0) {
    if (c == T) {
      lData <- CreateLagData(eData, index, lag, target)
      CompoundLinear(lData, scaler, index, target, lag, current, nTest)
    }
  }
  if (l == T) {
    LSTM(eData, scaler, index, target, lag, current, nTest, batch, step, epoch, numUnit, stop, split)
  }
}


### Simple Linear
SimpleLinear <- function(eData, scaler, city, target, lag, current, nTest) {
  all <- Standardize(eData, scaler[[1]], scaler[[2]], scaler[[3]])
  y <- all[, target] %>% tail(-lag)
  x <- all %>% head(-lag)
  if (current == F) {
    x <- x[, !(colnames(x) %in% target)]
  }
  train <- data.frame(x, y) %>% head(-nTest)
  test <- data.frame(x, y) %>% tail(nTest)
  model <- lm(y ~ ., train)
  assign(Texcat("model.", GetParamText("s")), model, envir = globalenv())
  saveRDS(model, Texcat("./model/", GetParamText("s"), ".RDS"))
  
  pred <- c(rep(NA, lag), (predict(model, x) * scaler[[1]][target] + scaler[[2]][target]) %>% exp())
  assign(Texcat("prdic.", GetParamText("s")), pred, envir = globalenv())
  saveRDS(pred, Texcat("./prdic/", GetParamText("s"), ".RDS"))
}

### Linear
CompoundLinear <- function(lData, scaler, index, target, lag, current, nTest) {
  all <- Standardize(lData, scaler[[1]], scaler[[2]], scaler[[3]])
  y <- all[, target] %>% tail(-lag)
  x <- all %>% head(-lag)
  if (current == F) {
    x <- x[, !(colnames(x) %in% target)]
  }
  
  train <- data.frame(x, y) %>% head(-nTest)
  test <- data.frame(x, y) %>% tail(nTest)
  model <- lm(y ~ ., train)
  assign(Texcat("model.", GetParamText("c")), model, envir = globalenv())
  saveRDS(model, Texcat("./model/", GetParamText("c"), ".RDS"))
  
  pred <- c(rep(NA, lag), (predict(model, x) * scaler[[1]][target] + scaler[[2]][target]) %>% exp())
  assign(Texcat("prdic.", GetParamText("ca")), pred, envir = globalenv())
  saveRDS(pred, Texcat("./prdic/", GetParamText("c"), ".RDS"))
}

### LSTM
LSTM <- function(eData, scaler, index, target, lag, current, nTest, batch, step, epoch, numUnit, stop, split) {
  # set random seed
  use_session_with_seed(144, disable_gpu = FALSE)
  
  all <- Standardize(eData, scaler[[1]], scaler[[2]], scaler[[3]])
  sampleDiff <- ifelse(dim(all)[1] %% batch > lag + step - 1,
                       dim(all)[1] %% batch,
                       (dim(all)[1] - (lag + step - 1)) %% batch + (lag + step - 1)
  )
  ln <- dim(all)[1] - sampleDiff
  variables <- colnames(all)
  if (current == F) {
    variables <- variables[-which(variables == target)]
  }
  d <- length(variables)
  
  train.x <- array(NA, c(ln - nTest, step, d))
  test.x <- array(NA, c(nTest, step, d))
  
  y <- tail(all[, target], ln)
  train.y <- head(y, -nTest)
  test.y <- tail(y, nTest)
  
  for (k in 1:length(variables)) {
    for (j in 1:step) {
      chunk <- unlist(all[, variables[k]] %>% tail(ln + lag + j - 1) %>% head(ln))
      train.x[, step - j + 1, k] <- chunk %>% head(-nTest)
      test.x[, step - j + 1, k] <- chunk %>% tail(nTest)
    }
  }
  
  model <- keras_model_sequential()
  
  model %>%
    bidirectional(
      layer_lstm(
      units = length(variables) * numUnit,
      return_sequences = FALSE,
      stateful = T
      ),
    input_shape = c(step, d),
    batch_size = batch) %>%
    layer_dropout(0.05) %>%
    layer_dense(units = 1)
  
  if (stop == 0) {
    history = model %>%
      compile(loss = "mean_squared_error", optimizer = "RMSProp") %>%
      fit(
        x = train.x,
        y = train.y,
        batch_size = batch,
        epochs = epoch,
        verbose = 2,
        shuffle = FALSE
      )
  } else {
    history = model %>%
      compile(loss = "mean_squared_error", optimizer = "RMSProp") %>%
      fit(
        x = train.x,
        y = train.y,
        batch_size = batch,
        epochs = epoch,
        verbose = 2,
        shuffle = FALSE,
        validation_split = ifelse(split == 0, batch / (ln - nTest), batch * split / (ln - nTest)),
        callbacks = callback_early_stopping(patience = stop)
      )
  }
  
  
  pred.train <- model %>%
    predict(train.x, batch_size = batch) %>%
    .[, 1]
  pred.test <- model %>%
    predict(test.x, batch_size = batch) %>%
    .[, 1]
  pred <- c(rep(NA, dim(all)[1] - ln), pred.train, pred.test) * scaler[[1]][target] %>% as.numeric() + scaler[[2]][target] %>% as.numeric()
  pred <- pred %>% exp()
  
  # assign(Texcat("model.", GetParamText("l")), model, envir = globalenv())
  # assign(Texcat("prdic.", GetParamText("l")), pred, envir = globalenv())
  
  
  save_model_hdf5(
    model
    , Texcat("./model/", GetParamText("l"), ".hdf5"),
    overwrite = TRUE,
    include_optimizer = TRUE
  )
  saveRDS(pred, Texcat("./prdic/", GetParamText("l"), ".RDS"))
  saveRDS(history, Texcat("./history/", GetParamText("l"), ".RDS"))
  k_clear_session()
}

MakeFirstData()
# variables for variable selection
{
  CO <- c("CO.1", "CO.2", "CO.3", "CO.4", "CO.5", "CO.6", "CO.7", "CO.8", "CO.9", "CO.10", "CO.11", "CO.12", "CO.13")
  NO2 <- c("NO2.1", "NO2.2", "NO2.3", "NO2.4", "NO2.5", "NO2.6", "NO2.7", "NO2.8", "NO2.9", "NO2.10", "NO2.11", "NO2.12", "NO2.13")
  SO2 <- c("SO2.1", "SO2.2", "SO2.3", "SO2.4", "SO2.5", "SO2.6", "SO2.7", "SO2.8", "SO2.9", "SO2.10", "SO2.11", "SO2.12", "SO2.13")
  O3 <- c("O3.1", "O3.2", "O3.3", "O3.4", "O3.5", "O3.6", "O3.7", "O3.8", "O3.9", "O3.10", "O3.11", "O3.12", "O3.13")
  PM2.5 <- c("PM2.5.1", "PM2.5.2", "PM2.5.3", "PM2.5.4", "PM2.5.5", "PM2.5.6", "PM2.5.7", "PM2.5.8", "PM2.5.9", "PM2.5.10", "PM2.5.11", "PM2.5.12", "PM2.5.13")
  PM10 <- c("PM10.1", "PM10.2", "PM10.3", "PM10.4", "PM10.5", "PM10.6", "PM10.7", "PM10.8", "PM10.9", "PM10.10", "PM10.11", "PM10.12", "PM10.13")
  AP <- c("AP.1", "AP.2", "AP.3", "AP.4", "AP.5", "AP.6", "AP.7", "AP.8", "AP.9", "AP.10", "AP.11", "AP.12", "AP.13")
  HM <- c("HM.1", "HM.2", "HM.3", "HM.4", "HM.5", "HM.6", "HM.7", "HM.8", "HM.9", "HM.10", "HM.11", "HM.12", "HM.13")
  TP <- c("TP.1", "TP.2", "TP.3", "TP.4", "TP.5", "TP.6", "TP.7", "TP.8", "TP.9", "TP.10", "TP.11", "TP.12", "TP.13")
  CW <- c("CW.1", "CW.2", "CW.3", "CW.4", "CW.5", "CW.6", "CW.7", "CW.8", "CW.9", "CW.10", "CW.11", "CW.12", "CW.13")
  WI <- c("WI.1", "WI.2", "WI.3", "WI.4", "WI.5", "WI.6", "WI.7", "WI.8", "WI.9", "WI.10", "WI.11", "WI.12", "WI.13")
  DR <- c("DR.1", "DR.2", "DR.3", "DR.4", "DR.5", "DR.6", "DR.7", "DR.8", "DR.9", "DR.10", "DR.11", "DR.12", "DR.13")
}

foreach(pol = pollutantNames, .combine = "cfun") %do% {
  # foreach(city = cityNames, .combine = "cfun") %do% {
  city <- "Beijing"
  # pol <- "NO2"
  param <- c("AP.S", "HM.S", "TP.S", "CW.A",  Texcat(pol, ".A"))
  SetQuery(city, pol, 1, T, param, 576, 64, 1, 30, 16,0,0)
  DoRegressions(city, pol, lag, current , param, nTest, batch, step, epoch, numUnit, stop, split, s = F, c=T, l=F)
  SetQuery(city, pol, 8, T, param, 576, 64, 8, 30, 16, 0, 0)
  DoRegressions(city, pol, lag, current , param, nTest, batch, step, epoch, numUnit, stop, split, s = F, c=T, l=F)
  SetQuery(city, pol, 24, T, param, 576, 64, 24, 30, 16, 0,0)
  DoRegressions(city, pol, lag, current , param, nTest, batch, step, epoch, numUnit, stop, split, s = F, c=T, l=F)
  # }
}

PrintPred = function(df,str){
  graph = ggplot(df) + geom_line(aes(x = Time, y = value, color = variable)) + BasicThemeNoLegend +
    labs(y = UnitsPol[which(pollutantNames == pol)]) +
    scale_color_manual(values=c( "#BB933E", "#255B87","#CC5287", "#93CC52"),
                       label = c("Actual","Linear", "C. Linear", "LSTM-RNN")) +
    scale_x_datetime(breaks = seq(as.POSIXct("2016-12-01"), as.POSIXct("2017-01-01"), by="1 week"), date_labels = "%b %d")
  return(graph)
}

PrintQQ = function(df,str){
  graph = ggplot(df) + geom_point(aes(x = actual, y = value, shape = variable, color = variable), alpha = 0.5) +
    scale_color_manual(values=c( "#255B87", "#CC5287","#93CC52"),
                       label = c("Linear", "C. Linear", "LSTM-RNN")) +
    scale_shape_manual(values=15:17,
                       label = c("Linear", "C. Linear", "LSTM-RNN")) +
    geom_abline(intercept = 0, slope = 1, color = "#BB933E") + BasicThemeNoLegend + 
    labs(x = Texcat("Actual ", pol), y = Texcat("Predicted ", pol)) +
    guides(shape = guide_legend(title = "Model", nrow = 1), col = guide_legend(title = "Model", nrow = 1)) 
  return(graph)
}

mape = data.frame(c(" ",pollutantNames),
          c("Linear", rep("",length(pollutantNames))),
          c("C. Linear", rep("",length(pollutantNames))),
          c("LSTM-RNN", rep("",length(pollutantNames))),
          stringsAsFactors = FALSE)
rmse = mape



st = c(16,32,64)
ep = c(100,30,10)
sp = c(30,15,5)

foreach(l = c(1, 8, 24), .combine = "cfun") %do% {
  i <- 1
  foreach(pol = pollutantNames, .combine = "cfun") %do% {
    index <- which(l == c(1, 8, 24))

    # param <- c("AP.S", "HM.S", "TP.S", "CW.A", Texcat(pol, ".A"))
    param <- c("AP.S", "HM.S", "TP.S", "CW.A",  Texcat(pollutantNames, ".A"))

    SetQuery(city, pol, l, T, param, 576, 64, l, 100, 16, sp[index], 9)
    # SetQuery(city, pol, 24, T, param, 576, 64, 8, 10, 16, 0, 0)
    # SetQuery(city, pol, l, T, param, 576, st[index], l, ep[index], 16, 0, 0)
    # SetQuery(city, pol, l, T, param, 576, 64, l, 30, 16, 1, 9)
    # SetQuery(city, pol, l, T, param, 576, 64, l, 30, 16, 0, 0)
    # SetQuery(city, pol, l, T, param, 576, 576, l, 100, 16, 0, 0)
    # SetQuery(city, pol, l, T, param, 576, 32, l, 30, 16,0,0)


    index <- which(pol == pollutantNames) + 1
    Time <- t %>% tail(576)
    actual <- get(city)[, pol] %>% tail(576)
    Linear <- readRDS(Texcat("./prdic/", GetParamText("s"), ".RDS")) %>% tail(576)
    Compound_Linear <- readRDS(Texcat("./prdic/", GetParamText("c"), ".RDS")) %>% tail(576)
    LSTM_RNN <- readRDS(Texcat("./prdic/", GetParamText("l"), ".RDS")) %>% tail(576)


    print(
      c(
        d(actual, Linear),
        d(actual, Compound_Linear),
        d(actual, LSTM_RNN)
      )
    )
    print(
      c(
        sprintf("%2.2f", rmse(actual, Linear)),
        sprintf("%2.2f", rmse(actual, Compound_Linear)),
        sprintf("%2.2f", rmse(actual, LSTM_RNN))
      )
    )

    # sm = readRDS(Texcat("./model/", GetParamText("s"), ".RDS"))
    # cm = readRDS(Texcat("./model/", GetParamText("c"), ".RDS"))
    #
    # # print(summary(sm))
    # print(summary(cm))
    # #
    # df <- data.frame(Time, actual, Linear) %>% reshape2::melt(id = "Time")
    # PrintPred(df,"s")
    # df <- data.frame(actual, Linear)%>% reshape2::melt(id = "actual")
    # PrintQQ(df,"s")
    #
    # df <- data.frame(Time, actual, Compound_Linear) %>% reshape2::melt(id = "Time")
    # PrintPred(df,"c")
    # df <- data.frame(actual, Compound_Linear)%>% reshape2::melt(id = "actual")
    # PrintQQ(df,"c")
    #
    # df <- data.frame(Time, actual, LSTM_RNN) %>% reshape2::melt(id = "Time")
    # PrintPred(df,"l")
    # df <- data.frame(actual, LSTM_RNN)%>% reshape2::melt(id = "actual")
    # PrintQQ(df,"l")

    # df <- data.frame(Time, actual, Linear, Compound_Linear, LSTM_RNN) %>% reshape2::melt(id = "Time")
    # assign(Texcat("graph", i),
    #   PrintPred(df, "a"),
    #   envir = globalenv()
    # )
    # df <- data.frame(actual, Linear, Compound_Linear, LSTM_RNN) %>% reshape2::melt(id = "actual")
    # assign(Texcat("graphQ", i),
    #   PrintQQ(df, "a"),
    #   envir = globalenv()
    # )
    # i <- i + 1
  }

  # graph0 <- cowplot::plot_grid(plotlist = mget(Texcat("graph", 1:6)), ncol = 2)
  # legend <- cowplot::get_legend(graph1 + BasicThemeWithLegend +
  #   theme(legend.text = element_text(size = rel(0.95), margin = margin(l = -3, r = 5))) +
  #   guides(shape = guide_legend(title = "Model", nrow = 1), col = guide_legend(title = "Model", nrow = 1)))
  # graph <- cowplot::plot_grid(graph0, legend, ncol = 1, rel_heights = c(36, 1))
  # ggsave(Texcat("./prdic_image/", GetParamText(str), ".png"), graph, width = 9, height = 6, dpi = 200, units = "in")
  # 
  # graph0 <- cowplot::plot_grid(plotlist = mget(Texcat("graphQ", 1:6)), ncol = 2)
  # legend <- cowplot::get_legend(graphQ1 + BasicThemeWithLegend +
  #   theme(legend.text = element_text(size = rel(0.95), margin = margin(l = -3, r = 5))) +
  #   guides(shape = guide_legend(title = "Model", nrow = 1), col = guide_legend(title = "Model", nrow = 1)))
  # graph <- cowplot::plot_grid(graph0, legend, ncol = 1, rel_heights = c(36, 1))
  # ggsave(Texcat("./qq_image/", GetParamText(str), ".png"), graph, width = 9, height = 6, dpi = 200, units = "in")
}


foreach(l = c(1, 8, 24), .combine = "cfun") %do% {
  i <- 1
  cat(l)
  foreach(pol = pollutantNames, .combine = "cfun") %do% {
    cat(pol)
    index <- which(l == c(1, 8, 24))
    
    param <- c("AP.S", "HM.S", "TP.S", "CW.A", Texcat(pol, ".A"))
    SetQuery(city, pol, l, T, param, 576, 64, l, 100, 16, sp[index], 9)
    
    actual <- get(city)[, pol] %>% tail(576)
    
    Linear_S <- (readRDS(Texcat("./prdic/", GetParamText("s"), ".RDS")) %>% tail(576) - actual)^2
    Compound_Linear_S <- (readRDS(Texcat("./prdic/", GetParamText("c"), ".RDS")) %>% tail(576) - actual)^2
    LSTM_RNN_S <- (readRDS(Texcat("./prdic/", GetParamText("l"), ".RDS")) %>% tail(576) - actual)^2
    
    param <- c("AP.S", "HM.S", "TP.S", "CW.A",  Texcat(pollutantNames, ".A"))
    SetQuery(city, pol, l, T, param, 576, 64, l, 100, 16, sp[index], 9)
    
    Linear_A <- (readRDS(Texcat("./prdic/", GetParamText("s"), ".RDS")) %>% tail(576) - actual)^2
    Compound_Linear_A <- (readRDS(Texcat("./prdic/", GetParamText("c"), ".RDS")) %>% tail(576) - actual)^2
    LSTM_RNN_A <- (readRDS(Texcat("./prdic/", GetParamText("l"), ".RDS")) %>% tail(576) - actual)^2
    
    print(c(
    wilcox.exact(Linear_A,Linear_S, alternative = "l", paired = T) %>% .$p.value < 0.05,
    wilcox.exact(Compound_Linear_A,Compound_Linear_S, alternative = "l", paired = T) %>% .$p.value < 0.05,
    wilcox.exact(LSTM_RNN_A,LSTM_RNN_S, alternative = "l", paired = T) %>% .$p.value < 0.05
    ))
    
    param <- c("AP.S", "HM.S", "TP.S", "CW.A", Texcat(pol, ".A"))
    SetQuery(city, pol, l, T, param, 576, 64, l, 100, 16, sp[index], 9)
    
    actual <- get(city)[, pol] %>% tail(576)
    
    Linear_S <- readRDS(Texcat("./prdic/", GetParamText("s"), ".RDS")) %>% tail(576)
    Compound_Linear_S <- readRDS(Texcat("./prdic/", GetParamText("c"), ".RDS")) %>% tail(576)
    LSTM_RNN_S <- readRDS(Texcat("./prdic/", GetParamText("l"), ".RDS")) %>% tail(576)
    
    param <- c("AP.S", "HM.S", "TP.S", "CW.A",  Texcat(pollutantNames, ".A"))
    SetQuery(city, pol, l, T, param, 576, 64, l, 100, 16, sp[index], 9)
    
    Linear_A <- readRDS(Texcat("./prdic/", GetParamText("s"), ".RDS")) %>% tail(576)
    Compound_Linear_A <- readRDS(Texcat("./prdic/", GetParamText("c"), ".RDS")) %>% tail(576) 
    LSTM_RNN_A <- readRDS(Texcat("./prdic/", GetParamText("l"), ".RDS")) %>% tail(576)
    print(
      c(
        sprintf("%2.1f", 100*(d(actual, Linear_A) - d(actual, Linear_S))/d(actual, Linear_S)),
        sprintf("%2.1f", 100*(d(actual, Compound_Linear_A)- d(actual, Compound_Linear_S))/d(actual, Compound_Linear_S)),
        sprintf("%2.1f", 100*(d(actual, LSTM_RNN_A)- d(actual, LSTM_RNN_S))/d(actual, LSTM_RNN_S))
      )
    )
    
    print(
      c(
        sprintf("%2.1f", 100*(rmse(actual, Linear_A) - rmse(actual, Linear_S))/rmse(actual, Linear_S)),
        sprintf("%2.1f", 100*(rmse(actual, Compound_Linear_A)- rmse(actual, Compound_Linear_S))/rmse(actual, Compound_Linear_S)),
        sprintf("%2.1f", 100*(rmse(actual, LSTM_RNN_A)- rmse(actual, LSTM_RNN_S))/rmse(actual, LSTM_RNN_S))
      )
    )
  }
}

foreach(l = c(1, 8, 24), .combine = "cfun") %do% {
  i <- 1
  cat(l)
  foreach(pol = pollutantNames, .combine = "cfun") %do% {
    param <- c("AP.S", "HM.S", "TP.S", "CW.A",  Texcat(pollutantNames, ".A"))
    SetQuery(city, pol, l, T, param, 576, 64, l, 100, 16, sp[index], 9)
    LSTM_RNN_M <- load_model_hdf5(Texcat("./model/", GetParamText("l"), ".hdf5"))
    
  }
}