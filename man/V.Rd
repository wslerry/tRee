\name{V}
\alias{V}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Volume prediction model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
V(data, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function (data, model)
{
    library(randomForest)
    library(e1071)
    library(caret)
    library(ggplot2)
    ggplotRegression <- function(fit) {
        a <- signif(coef(fit)[1], digits = 5)
        b <- signif(coef(fit)[2], digits = 5)
        if (coef(fit)[2] >= 0) {
            textlab <- paste("y = ", a, " + ", b, "x", sep = "")
        }
        else {
            textlab <- paste("y = ", a, " - ", b, "x", sep = "")
        }
        options(repr.plot.width = 4, repr.plot.height = 4)
        ggplot(fit$model, aes_string(x = names(fit$model)[1],
            y = names(fit$model)[2])) + geom_point() + geom_smooth(method = "lm",
            col = "red", size = 0.5, se = TRUE) + labs(x = "Observations",
            y = "Predictions", title = paste("Adj. R2 = ", signif(summary(fit)$adj.r.squared,
                5), " | ", textlab)) + theme(plot.title = element_text(size = 8,
            face = "bold"))
    }
    traindata <- base::sample(nrow(data), size = 0.7 * nrow(data),
        replace = FALSE)
    TrainSet <- data[traindata, ]
    ValidSet <- data[-traindata, ]
    ctrl <- caret::trainControl(method = "cv", number = 5, savePredictions = TRUE)
    if (model == "mod1") {
        mod <- train(log(VUB) ~ I(log(DBH)) + I(log(H)) + I(log(BA)),
            data = TrainSet, method = "lm", trControl = ctrl,
            metric = "Rsquared")
        predictions <- predict(mod, ValidSet)
        predicted_V <- data.frame(log(ValidSet$VUB), predictions)
        fit.mod <- lm(log(ValidSet$VUB) ~ predictions, data = predicted_V)
        g <- ggplotRegression(fit.mod)
        print(mod$finalModel)
        g
    }
    else if (model == "mod2") {
        mod <- train(log(VUB) ~ I(log(DBH)) + I(log(H)) + I(log(DBH)^2) +
            I(log(H)^2), data = TrainSet, method = "lm", trControl = ctrl,
            metric = "Rsquared")
        predictions <- predict(mod, ValidSet)
        predicted_V <- data.frame(log(ValidSet$VUB), predictions)
        fit.mod <- lm(log(ValidSet$VUB) ~ predictions, data = predicted_V)
        g <- ggplotRegression(fit.mod)
        print(mod$finalModel)
        g
    }
    else if (model == "mod3") {
        mod <- train(log(VUB) ~ I(log(DBH^2 * H)), data = TrainSet,
            method = "lm", trControl = ctrl, metric = "Rsquared")
        predictions <- predict(mod, ValidSet)
        predicted_V <- data.frame(predictions, log(ValidSet$VUB))
        fit.mod <- lm(log(ValidSet$VUB) ~ predictions, data = predicted_V)
        g <- ggplotRegression(fit.mod)
        print(mod$finalModel)
        g
    }
    else if (model == "mod4") {
        set.seed(10)
        mod <- train(log(VUB) ~ I(log(DBH^2)) + I(log(H^2)),
            data = TrainSet, method = "lm", trControl = ctrl,
            metric = "Rsquared")
        predictions <- predict(mod, ValidSet)
        predicted_V <- data.frame(predictions, log(ValidSet$VUB))
        fit.mod <- lm(log(ValidSet$VUB) ~ predictions, data = predicted_V)
        g <- ggplotRegression(fit.mod)
        print(mod$finalModel)
        g
    }
    else if (model == "mod5") {
        set.seed(10)
        mod <- train(log(VUB) ~ I(log(H)), data = TrainSet, method = "lm",
            trControl = ctrl, metric = "Rsquared")
        predictions <- predict(mod, ValidSet)
        predicted_V <- data.frame(predictions, log(ValidSet$VUB))
        fit.mod <- lm(log(ValidSet$VUB) ~ predictions, data = predicted_V)
        g <- ggplotRegression(fit.mod)
        print(mod$finalModel)
        g
    }
    else if (model == "modRF") {
        set.seed(10)
        mod <- train(log(VUB) ~ I(log(H)), data = TrainSet, method = "rf",
            ntree = 1000, trControl = ctrl, metric = "Rsquared",
            importance = T, proximity = T)
        predictions <- predict(mod, ValidSet)
        predicted <- data.frame(predictions, log(ValidSet$VUB))
        fit.mod <- lm(log(ValidSet$VUB) ~ predictions, data = predicted)
        g <- ggplotRegression(fit.mod)
        fit.mod$coefficients
        print(mod$finalModel)
        g
    }
    else if (model == "modRF2") {
        rf <- randomForest(VUB ~ H, TrainSet, ntree = 1000, importance = T,
            proximity = TRUE, keep.forest = TRUE)
        print(rf)
        plot(rf, main = "Random Forest")
        predictions <- predict(rf, ValidSet)
        predicted <- data.frame(log(predictions), log(ValidSet$VUB))
        fit.mod <- lm(log(ValidSet$VUB) ~ log(predictions), data = predicted)
        plot(predicted) + abline(fit.mod)
    }
    else {
        stop("No Model selected. Please select model eg., 'mod1','mod2','mod3','mod4','mod5' or 'modRF'")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
