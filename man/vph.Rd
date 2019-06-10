\name{vph}
\alias{vph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Volume Per Hectares Prediction Model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vph(data, img, poly, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{img}{
%%     ~~Describe \code{img} here~~
}
  \item{poly}{
%%     ~~Describe \code{poly} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, img, poly, model)
{
    require(raster)
    require(rgdal)
    require(rgeos)
    traindata <- base::sample(nrow(data), size = 0.7 * nrow(data),
        replace = FALSE)
    TrainSet <- data[traindata, ]
    ValidSet <- data[-traindata, ]
    if (model == "mod1") {
        mod <- lm(log(VUB) ~ I(log(DBH)) + I(log(H)) + I(log(BA)),
            data = data)
        coef(mod)
        vol <- function(x) {
            exp(coef(mod)[1]) * (23.1^coef(mod)[2]) * (x^coef(mod)[3]) *
                (0.046^coef(mod)[4])
        }
        vol.agg <- aggregate(img, fact = 4/res(img))
        volimg <- raster::calc(vol.agg, vol)
        x <- xres(volimg)
        y <- yres(volimg)
        ha <- (x * y)/10000
        volha.func <- function(x) {
            x/ha
        }
        volha <- raster::calc(volimg, volha.func)
        par(mar = rep(0.5, 4))
        compsRas <- rasterize(poly, volha)
        zoneStat <- zonal(volha, compsRas, "mean")
        poly[["volumeperha"]] <- zoneStat[, "mean"]
        colRamp <- colorRampPalette(c("lightgoldenrod1", "tomato2"))(10)
        polyCols <- colRamp[as.numeric(cut(poly[["volumeperha"]],
            breaks = 10))]
        options(repr.plot.width = 4, repr.plot.height = 4)
        plot(poly, col = polyCols, xlab = "", ylab = "", xaxt = "n",
            yaxt = "n", main = "Volume/Ha")
        text(gCentroid(poly, byid = TRUE), round(poly[["volumeperha"]],
            1), col = "blue", font = 2)
        writeRaster(volha, "./output/volperha_mod1.tif", format = "GTiff",
            overwrite = TRUE)
        writeOGR(poly, "output", "vhamod1", driver = "ESRI Shapefile",
            overwrite = TRUE)
    }
    else if (model == "mod2") {
        mod <- lm(log(VUB) ~ I(log(DBH)) + I(log(H)) + I(log(DBH)^2) +
            I(log(H)^2), data = data)
        coef(mod)
        vol <- function(x) {
            exp(coef(mod)[1]) * (23.1^coef(mod)[2]) * (x^coef(mod)[3]) *
                ((23.1^2)^coef(mod)[4]) * ((x^2)^coef(mod)[5])
        }
        vol.agg <- aggregate(img, fact = 3/res(img))
        volimg <- raster::calc(vol.agg, vol)
        x <- xres(volimg)
        y <- yres(volimg)
        ha <- (x * y)/10000
        volha.func <- function(x) {
            x/ha
        }
        volha <- raster::calc(volimg, volha.func)
        par(mar = rep(0.5, 4))
        compsRas <- rasterize(poly, volha)
        zoneStat <- zonal(volha, compsRas, "mean")
        poly[["volumeperha"]] <- zoneStat[, "mean"]
        colRamp <- colorRampPalette(c("lightgoldenrod1", "tomato2"))(10)
        polyCols <- colRamp[as.numeric(cut(poly[["volumeperha"]],
            breaks = 10))]
        plot(poly, col = polyCols, xlab = "", ylab = "", xaxt = "n",
            yaxt = "n", main = "Volume/Ha")
        text(gCentroid(poly, byid = TRUE), round(poly[["volumeperha"]],
            1), col = "blue", font = 2)
        writeRaster(volha, "./output/volperha_mod2.tif", format = "GTiff",
            overwrite = TRUE)
        writeOGR(poly, "output", "vhamod2", driver = "ESRI Shapefile",
            overwrite = TRUE)
    }
    else if (model == "mod3") {
        mod <- lm(log(VUB) ~ I(log(DBH^2 * H)), data = data)
        coef(mod)
        vol <- function(x) {
            exp(coef(mod)[1]) * (((23.1^2) * x)^coef(mod)[2])
        }
        vol.agg <- aggregate(img, fact = 4/res(img))
        volimg <- raster::calc(vol.agg, vol)
        x <- xres(volimg)
        y <- yres(volimg)
        ha <- (x * y)/10000
        volha.func <- function(x) {
            x/ha
        }
        volha <- raster::calc(volimg, volha.func)
        par(mar = rep(0.5, 4))
        compsRas <- rasterize(poly, volha)
        zoneStat <- zonal(volha, compsRas, "mean")
        poly[["volumeperha"]] <- zoneStat[, "mean"]
        colRamp <- colorRampPalette(c("lightgoldenrod1", "tomato2"))(10)
        polyCols <- colRamp[as.numeric(cut(poly[["volumeperha"]],
            breaks = 10))]
        plot(poly, col = polyCols, xlab = "", ylab = "", xaxt = "n",
            yaxt = "n", main = "Volume/Ha")
        text(gCentroid(poly, byid = TRUE), round(poly[["volumeperha"]],
            1), col = "blue", font = 2)
        writeRaster(volha, "./output/volperha_mod3.tif", format = "GTiff",
            overwrite = TRUE)
        writeOGR(poly, "output", "vhamod3", driver = "ESRI Shapefile",
            overwrite = TRUE)
    }
    else if (model == "mod4") {
        mod <- lm(log(VUB) ~ I(log(DBH^2)) + I(log(H^2)), data = data)
        coef(mod)
        vol <- function(x) {
            exp(coef(mod)[1]) * ((23.1^2)^coef(mod)[2]) * ((x^2)^coef(mod)[3])
        }
        vol.agg <- aggregate(img, fact = 4/res(img))
        volimg <- raster::calc(vol.agg, vol)
        x <- xres(volimg)
        y <- yres(volimg)
        ha <- (x * y)/10000
        volha.func <- function(x) {
            x/ha
        }
        volha <- raster::calc(volimg, volha.func)
        par(mar = rep(0.5, 4))
        compsRas <- rasterize(poly, volha)
        zoneStat <- zonal(volha, compsRas, "mean")
        poly[["volumeperha"]] <- zoneStat[, "mean"]
        colRamp <- colorRampPalette(c("lightgoldenrod1", "tomato2"))(10)
        polyCols <- colRamp[as.numeric(cut(poly[["volumeperha"]],
            breaks = 10))]
        plot(poly, col = polyCols, xlab = "", ylab = "", xaxt = "n",
            yaxt = "n", main = "Volume/Ha")
        text(gCentroid(poly, byid = TRUE), round(poly[["volumeperha"]],
            1), col = "blue", font = 2)
        writeRaster(volha, "./output/volperha_mod4.tif", format = "GTiff",
            overwrite = TRUE)
        writeOGR(poly, "output", "vhamod4", driver = "ESRI Shapefile",
            overwrite = TRUE)
    }
    else if (model == "mod5") {
        mod <- lm(log(VUB) ~ I(log(H)), data = data)
        coef(mod)
        vol <- function(x) {
            exp(coef(mod)[1]) * (x^coef(mod)[2])
        }
        vol.agg <- aggregate(img, fact = 4/res(img))
        volimg <- raster::calc(vol.agg, vol)
        x <- xres(volimg)
        y <- yres(volimg)
        ha <- (x * y)/10000
        volha.func <- function(x) {
            x/ha
        }
        volha <- raster::calc(volimg, volha.func)
        par(mar = rep(0.5, 4))
        compsRas <- rasterize(poly, volha)
        zoneStat <- zonal(volha, compsRas, "mean")
        poly[["volumeperha"]] <- zoneStat[, "mean"]
        colRamp <- colorRampPalette(c("lightgoldenrod1", "tomato2"))(10)
        polyCols <- colRamp[as.numeric(cut(poly[["volumeperha"]],
            breaks = 10))]
        plot(poly, col = polyCols, xlab = "", ylab = "", xaxt = "n",
            yaxt = "n", main = "Volume/Ha")
        text(gCentroid(poly, byid = TRUE), round(poly[["volumeperha"]],
            1), col = "blue", font = 2)
        writeRaster(volha, "./output/volperha_mod5.tif", format = "GTiff",
            overwrite = TRUE)
        writeOGR(poly, "output", "vhamod5", driver = "ESRI Shapefile",
            overwrite = TRUE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
