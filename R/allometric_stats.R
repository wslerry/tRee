# require(rgdal)
# require(raster)
# require(rgeos)
# require(PerformanceAnalytics)
# require(randomForest)
# require(e1071)
# require(caret)
# require(ggplot2)

readvector <- function(vectorpath){
  require(rgdal)
  readOGR(vectorpath)
}

readraster <- function(rasterdir){
  require(raster)
  raster(rasterdir)
}

datacorrelation<-function(data){
  require(PerformanceAnalytics)
  res <- cor(data)
  print(round(res, 2))
  options(repr.plot.width=10, repr.plot.height=10)
  chart.Correlation(data, histogram=TRUE, method = c("pearson","kendall","spearman"),pch=19)
}

V <- function(data,model){
  library(randomForest)
  library(e1071)
  library(caret)
  library(ggplot2)

  #parallelCluster <- parallel::makeCluster(parallel::detectCores())
  #print(parallelCluster)

  ggplotRegression <- function (fit) {
    a <- signif(coef(fit)[1], digits = 5)
    b <- signif(coef(fit)[2], digits = 5)
    if (coef(fit)[2] >= 0)  {
      textlab <- paste("y = ",a," + ",b,"x", sep="")
    } else {
      textlab <- paste("y = ",a," - ",b,"x", sep="")
    }
    options(repr.plot.width=4, repr.plot.height=4)
    ggplot(fit$model,aes_string(x = names(fit$model)[1], y = names(fit$model)[2])) +
      geom_point() +
      geom_smooth(method = "lm", col = "red",size = 0.5, se=TRUE) +
      labs(x="Observations",
           y = "Predictions",
           title = paste("Adj. R2 = ",signif(summary(fit)$adj.r.squared, 5)," | ",
                         textlab
           )
      )+theme(plot.title = element_text(size = 8, face = "bold"))
  }

  traindata <- base::sample(nrow(data), size=0.7*nrow(data), replace = FALSE)

  TrainSet <- data[traindata,]
  ValidSet <- data[-traindata,]
  ctrl <- caret::trainControl(method = "cv", number = 5, savePredictions = TRUE)

  if (model== 'mod1'){
    mod<-train(log(VUB)~I(log(DBH))+I(log(H))+I(log(BA)),
               data=TrainSet, method="lm", trControl = ctrl, metric="Rsquared")
    predictions <- predict(mod, ValidSet)
    predicted_V <- data.frame(log(ValidSet$VUB),predictions)
    fit.mod<-lm(log(ValidSet$VUB)~predictions,data=predicted_V)
    g<-ggplotRegression(fit.mod)
    print(mod$finalModel)
    g

  } else if (model == 'mod2'){
    mod<-train(log(VUB)~I(log(DBH))+I(log(H))+I(log(DBH)^2)+I(log(H)^2),
               data=TrainSet, method="lm", trControl = ctrl, metric="Rsquared")
    predictions <- predict(mod, ValidSet)
    predicted_V <- data.frame(log(ValidSet$VUB),predictions)
    fit.mod<-lm(log(ValidSet$VUB)~predictions,data=predicted_V)
    g<-ggplotRegression(fit.mod)
    print(mod$finalModel)
    g
  } else if (model == 'mod3') {
    mod<-train(log(VUB)~I(log(DBH^2*H)),
               data=TrainSet, method="lm", trControl = ctrl, metric="Rsquared")
    predictions <- predict(mod, ValidSet)
    predicted_V <- data.frame(predictions,log(ValidSet$VUB))
    fit.mod<-lm(log(ValidSet$VUB)~predictions,data=predicted_V)
    g<-ggplotRegression(fit.mod)
    print(mod$finalModel)
    g

  } else if (model == 'mod4') {
    set.seed(10)
    mod<-train(log(VUB)~I(log(DBH^2))+I(log(H^2)),
               data=TrainSet, method="lm", trControl = ctrl, metric="Rsquared")
    predictions <- predict(mod, ValidSet)
    predicted_V <- data.frame(predictions,log(ValidSet$VUB))
    fit.mod<-lm(log(ValidSet$VUB)~predictions,data=predicted_V)
    g<-ggplotRegression(fit.mod)
    print(mod$finalModel)
    g

  } else if (model == 'mod5') {
    set.seed(10)
    mod<-train(log(VUB)~I(log(H)),
               data=TrainSet, method="lm", trControl = ctrl, metric="Rsquared")
    predictions <- predict(mod, ValidSet)
    predicted_V <- data.frame(predictions,log(ValidSet$VUB))
    fit.mod<-lm(log(ValidSet$VUB)~predictions,data=predicted_V)
    g<-ggplotRegression(fit.mod)
    print(mod$finalModel)
    g

  } else if (model == 'modRF') {
    set.seed(10)
    mod<-train(log(VUB)~I(log(H)),
               data=TrainSet, method="rf", ntree=1000,
               trControl = ctrl, metric="Rsquared",
               importance=T,proximity=T)
    predictions <- predict(mod, ValidSet)
    predicted <- data.frame(predictions,log(ValidSet$VUB))
    fit.mod<-lm(log(ValidSet$VUB)~predictions,data=predicted)
    g<-ggplotRegression(fit.mod)
    fit.mod$coefficients
    print(mod$finalModel)
    g

  } else if (model == 'modRF2') {
    rf <- randomForest(VUB~H,
                       TrainSet,
                       ntree=1000,
                       importance=T,proximity=TRUE,
                       keep.forest=TRUE)
    print(rf)
    plot(rf,main='Random Forest')

    predictions <- predict(rf, ValidSet)
    predicted <- data.frame(log(predictions), log(ValidSet$VUB))
    fit.mod<-lm(log(ValidSet$VUB)~log(predictions), data=predicted)
    plot(predicted)+abline(fit.mod)

  } else {stop("No Model selected. Please select model eg., 'mod1','mod2','mod3','mod4','mod5' or 'modRF'")}
}

DBH <- function(data){
  library(randomForest)
  library(caret)
  library(ggplot2)

  ggplotRegression <- function (fit) {
    a <- signif(coef(fit)[1], digits = 5)
    b <- signif(coef(fit)[2], digits = 5)
    if (coef(fit)[2] >= 0)  {
      textlab <- paste("y = ",a," + ",b,"x", sep="")
    } else {
      textlab <- paste("y = ",a," - ",b,"x", sep="")
    }
    options(repr.plot.width=5, repr.plot.height=5)
    ggplot(fit$model,aes_string(x = names(fit$model)[1], y = names(fit$model)[2])) +
      geom_point() +
      geom_smooth(method = "lm", col = "red",size = 0.5, se=TRUE) +
      labs(x="Observations",
           y = "Predictions",
           title = paste("Adj. R2 = ",signif(summary(fit)$adj.r.squared, 5)," | ",
                         textlab
           )
      )+theme(plot.title = element_text(size = 8, face = "bold"))
  }

  res <- cor(data)

  traindata <- base::sample(nrow(data), size=0.6*nrow(data), replace = FALSE)

  TrainSet <- data[traindata,]
  ValidSet <- data[-traindata,]
  ctrl <- caret::trainControl(method = "cv", number = 10, savePredictions = TRUE)

  mod<-train(log(DBH)~I(log(H)),
             data=TrainSet, method="lm",
             trControl = ctrl, metric="Rsquared")
  predictions <- predict(mod, ValidSet)
  predicted <- data.frame(log(ValidSet$DBH),predictions)
  fit.mod<-lm(log(ValidSet$DBH)~predictions,data=predicted)
  g<-ggplotRegression(fit.mod)
  print(round(res, 3))
  print(mod$finalModel)
  g
}

BA <- function(data){
  library(randomForest)
  library(caret)
  library(ggplot2)

  ggplotRegression <- function (fit) {
    a <- signif(coef(fit)[1], digits = 5)
    b <- signif(coef(fit)[2], digits = 5)
    if (coef(fit)[2] >= 0)  {
      textlab <- paste("y = ",a," + ",b,"x", sep="")
    } else {
      textlab <- paste("y = ",a," - ",b,"x", sep="")
    }
    options(repr.plot.width=4, repr.plot.height=4)
    ggplot(fit$model,aes_string(x = names(fit$model)[1], y = names(fit$model)[2])) +
      geom_point() +
      geom_smooth(method = "lm", col = "red",size = 0.5, se=TRUE) +
      labs(x="Observations",
           y = "Predictions",
           title = paste("Adj. R2 = ",signif(summary(fit)$adj.r.squared, 5)," | ",
                         textlab
           )
      )+theme(plot.title = element_text(size = 8, face = "bold"))
  }

  res <- cor(data)

  traindata <- base::sample(nrow(data), size=0.6*nrow(data), replace = FALSE)

  TrainSet <- data[traindata,]
  ValidSet <- data[-traindata,]
  ctrl <- caret::trainControl(method = "cv", number = 10, savePredictions = TRUE)

  mod<-train(log(BA)~I(log(H)),
             data=TrainSet, method="lm",
             trControl = ctrl, metric="Rsquared")
  predictions <- predict(mod, ValidSet)
  predicted_BA <- data.frame(log(ValidSet$BA),predictions)
  fit.mod<-lm(log(ValidSet$BA)~predictions,data=predicted_BA)
  g<-ggplotRegression(fit.mod)
  print(round(res, 3))
  print(mod$finalModel)
  g
}

vph <- function (data,img,poly,model) {
  require(raster)
  require(rgdal)
  require(rgeos)

  traindata <- base::sample(nrow(data), size=0.7*nrow(data), replace = FALSE)

  TrainSet <- data[traindata,]
  ValidSet <- data[-traindata,]

  if (model == 'mod1'){

    mod<-lm(log(VUB)~I(log(DBH))+I(log(H))+I(log(BA)),
            data=data)
    coef(mod)
    vol <- function(x) {
      exp(coef(mod)[1])*(23.1^coef(mod)[2])*(x^coef(mod)[3])*(0.046^coef(mod)[4])
    }

    vol.agg<- aggregate(img, fact = 4/res(img))
    #print(res(vol.agg))

    volimg <- raster::calc(vol.agg, vol)

    x<-xres(volimg)
    y<-yres(volimg)
    ha<-(x*y)/10000
    #print(paste(ha," Ha"))
    volha.func<- function(x){
      x/ha
    }
    #beginCluster()
    volha <- raster::calc(volimg,volha.func)
    #endCluster()

    par(mar = rep(0.5, 4))

    #plot(volha, xlab = "", ylab = "", nc='n', nr = 'n', main="Volume/Ha :\n Model 1")
    #plot(poly, add = TRUE, border =  "darkmagenta", lwd = 2)

    # rasterize compartments
    compsRas <- rasterize(poly, volha)

    # get mean statistic from raster
    zoneStat <- zonal(volha, compsRas, 'mean')

    # Create new 'topHeight' attribute from zonal statistics
    poly[["volumeperha"]] <- zoneStat[,"mean"]

    # Plot result
    colRamp <- colorRampPalette(c('lightgoldenrod1', 'tomato2'))(10)
    polyCols <- colRamp[as.numeric(cut(poly[["volumeperha"]],breaks = 10))]

    options(repr.plot.width=4, repr.plot.height=4)
    plot(poly, col = polyCols, xlab = "", ylab = "", xaxt='n', yaxt = 'n', main="Volume/Ha")
    text(gCentroid(poly, byid = TRUE), round(poly[["volumeperha"]],1), col = "blue", font = 2)

    writeRaster(volha, './output/volperha_mod1.tif', format = "GTiff", overwrite=TRUE)
    writeOGR(poly, "output", "vhamod1", driver ="ESRI Shapefile", overwrite=TRUE)
  }
  else if (model == 'mod2'){
    mod<-lm(log(VUB)~I(log(DBH))+I(log(H))+I(log(DBH)^2)+I(log(H)^2),
            data=data)
    coef(mod)
    vol <- function(x) {
      exp(coef(mod)[1])*(23.1^coef(mod)[2])*(x^coef(mod)[3])*((23.1^2)^coef(mod)[4])*((x^2)^coef(mod)[5])
    }

    vol.agg<- aggregate(img, fact = 3/res(img))
    #print(res(vol.agg))

    volimg <- raster::calc(vol.agg, vol)

    x<-xres(volimg)
    y<-yres(volimg)
    ha<-(x*y)/10000
    #print(paste(ha," Ha"))
    volha.func<- function(x){
      x/ha
    }
    #beginCluster()
    volha <- raster::calc(volimg,volha.func)
    #endCluster()

    par(mar = rep(0.5, 4))

    #plot(volha, xlab = "", ylab = "", nc='n', nr = 'n', main="Volume/Ha :\n Model 1")
    #plot(poly, add = TRUE, border =  "darkmagenta", lwd = 2)

    # rasterize compartments
    compsRas <- rasterize(poly, volha)

    # get mean statistic from raster
    zoneStat <- zonal(volha, compsRas, 'mean')

    # Create new 'topHeight' attribute from zonal statistics
    poly[["volumeperha"]] <- zoneStat[,"mean"]

    # Plot result
    colRamp <- colorRampPalette(c('lightgoldenrod1', 'tomato2'))(10)
    polyCols <- colRamp[as.numeric(cut(poly[["volumeperha"]],breaks = 10))]

    plot(poly, col = polyCols, xlab = "", ylab = "",xaxt='n', yaxt = 'n', main="Volume/Ha")
    text(gCentroid(poly, byid = TRUE), round(poly[["volumeperha"]],1), col = "blue", font = 2)

    writeRaster(volha, './output/volperha_mod2.tif', format = "GTiff", overwrite=TRUE)
    writeOGR(poly, "output", "vhamod2", driver ="ESRI Shapefile", overwrite=TRUE)
  }
  else if (model == 'mod3'){
    mod<-lm(log(VUB)~I(log(DBH^2 * H)),
            data=data)
    coef(mod)
    vol <- function(x) {
      exp(coef(mod)[1])*(((23.1^2) * x)^coef(mod)[2])
    }

    vol.agg<- aggregate(img, fact = 4/res(img))
    #print(res(vol.agg))

    volimg <- raster::calc(vol.agg, vol)

    x<-xres(volimg)
    y<-yres(volimg)
    ha<-(x*y)/10000
    #print(paste(ha," Ha"))
    volha.func<- function(x){
      x/ha
    }
    #beginCluster()
    volha <- raster::calc(volimg,volha.func)
    #endCluster()

    par(mar = rep(0.5, 4))

    #plot(volha, xlab = "", ylab = "", nc='n', nr = 'n', main="Volume/Ha :\n Model 1")
    #plot(poly, add = TRUE, border =  "darkmagenta", lwd = 2)

    # rasterize compartments
    compsRas <- rasterize(poly, volha)

    # get mean statistic from raster
    zoneStat <- zonal(volha, compsRas, 'mean')

    # Create new 'topHeight' attribute from zonal statistics
    poly[["volumeperha"]] <- zoneStat[,"mean"]

    # Plot result
    colRamp <- colorRampPalette(c('lightgoldenrod1', 'tomato2'))(10)
    polyCols <- colRamp[as.numeric(cut(poly[["volumeperha"]],breaks = 10))]

    plot(poly, col = polyCols, xlab = "", ylab = "",xaxt='n', yaxt = 'n', main="Volume/Ha")
    text(gCentroid(poly, byid = TRUE), round(poly[["volumeperha"]],1), col = "blue", font = 2)

    writeRaster(volha, './output/volperha_mod3.tif', format = "GTiff", overwrite=TRUE)
    writeOGR(poly, "output", "vhamod3", driver ="ESRI Shapefile", overwrite=TRUE)
  }
  else if (model == 'mod4'){
    mod<-lm(log(VUB)~I(log(DBH^2))+I(log(H^2)),
            data=data)
    coef(mod)
    vol <- function(x) {
      exp(coef(mod)[1])*((23.1^2)^coef(mod)[2])*((x^2)^coef(mod)[3])
    }

    vol.agg<- aggregate(img, fact = 4/res(img))
    #print(res(vol.agg))

    volimg <- raster::calc(vol.agg, vol)

    x<-xres(volimg)
    y<-yres(volimg)
    ha<-(x*y)/10000
    #print(paste(ha," Ha"))
    volha.func<- function(x){
      x/ha
    }
    #beginCluster()
    volha <- raster::calc(volimg,volha.func)
    #endCluster()

    par(mar = rep(0.5, 4))

    #plot(volha, xlab = "", ylab = "", nc='n', nr = 'n', main="Volume/Ha :\n Model 1")
    #plot(poly, add = TRUE, border =  "darkmagenta", lwd = 2)

    # rasterize compartments
    compsRas <- rasterize(poly, volha)

    # get mean statistic from raster
    zoneStat <- zonal(volha, compsRas, 'mean')

    # Create new 'topHeight' attribute from zonal statistics
    poly[["volumeperha"]] <- zoneStat[,"mean"]

    # Plot result
    colRamp <- colorRampPalette(c('lightgoldenrod1', 'tomato2'))(10)
    polyCols <- colRamp[as.numeric(cut(poly[["volumeperha"]],breaks = 10))]

    plot(poly, col = polyCols, xlab = "", ylab = "",xaxt='n', yaxt = 'n', main="Volume/Ha")
    text(gCentroid(poly, byid = TRUE), round(poly[["volumeperha"]],1), col = "blue", font = 2)

    writeRaster(volha, './output/volperha_mod4.tif', format = "GTiff", overwrite=TRUE)
    writeOGR(poly, "output", "vhamod4", driver ="ESRI Shapefile", overwrite=TRUE)
  }
  else if (model == 'mod5'){
    mod<-lm(log(VUB)~I(log(H)),
            data=data)
    coef(mod)
    vol <- function(x) {
      exp(coef(mod)[1])*(x^coef(mod)[2])
    }

    vol.agg<- aggregate(img, fact = 4/res(img))
    #print(res(vol.agg))

    volimg <- raster::calc(vol.agg, vol)

    x<-xres(volimg)
    y<-yres(volimg)
    ha<-(x*y)/10000
    #print(paste(ha," Ha"))
    volha.func<- function(x){
      x/ha
    }
    #beginCluster()
    volha <- raster::calc(volimg,volha.func)
    #endCluster()

    par(mar = rep(0.5,4))

    #plot(volha, xlab = "", ylab = "", nc='n', nr = 'n', main="Volume/Ha :\n Model 1")
    #plot(poly, add = TRUE, border =  "darkmagenta", lwd = 2)

    # rasterize compartments
    compsRas <- rasterize(poly, volha)

    # get mean statistic from raster
    zoneStat <- zonal(volha, compsRas, 'mean')

    # Create new 'topHeight' attribute from zonal statistics
    poly[["volumeperha"]] <- zoneStat[,"mean"]

    # Plot result
    colRamp <- colorRampPalette(c('lightgoldenrod1', 'tomato2'))(10)
    polyCols <- colRamp[as.numeric(cut(poly[["volumeperha"]],breaks = 10))]

    plot(poly, col = polyCols, xlab = "", ylab = "",xaxt='n', yaxt = 'n', main="Volume/Ha")
    text(gCentroid(poly, byid = TRUE), round(poly[["volumeperha"]],1), col = "blue", font = 2)

    writeRaster(volha, './output/volperha_mod5.tif', format = "GTiff", overwrite=TRUE)
    writeOGR(poly, "output", "vhamod5", driver ="ESRI Shapefile", overwrite=TRUE)
  }
}
