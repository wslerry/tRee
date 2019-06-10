\name{DBH}
\alias{DBH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Diameter at breast height prediction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DBH(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function (data)
{
    library(randomForest)
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
        options(repr.plot.width = 5, repr.plot.height = 5)
        ggplot(fit$model, aes_string(x = names(fit$model)[1],
            y = names(fit$model)[2])) + geom_point() + geom_smooth(method = "lm",
            col = "red", size = 0.5, se = TRUE) + labs(x = "Observations",
            y = "Predictions", title = paste("Adj. R2 = ", signif(summary(fit)$adj.r.squared,
                5), " | ", textlab)) + theme(plot.title = element_text(size = 8,
            face = "bold"))
    }
    res <- cor(data)
    traindata <- base::sample(nrow(data), size = 0.6 * nrow(data),
        replace = FALSE)
    TrainSet <- data[traindata, ]
    ValidSet <- data[-traindata, ]
    ctrl <- caret::trainControl(method = "cv", number = 10, savePredictions = TRUE)
    mod <- train(log(DBH) ~ I(log(H)), data = TrainSet, method = "lm",
        trControl = ctrl, metric = "Rsquared")
    predictions <- predict(mod, ValidSet)
    predicted <- data.frame(log(ValidSet$DBH), predictions)
    fit.mod <- lm(log(ValidSet$DBH) ~ predictions, data = predicted)
    g <- ggplotRegression(fit.mod)
    print(round(res, 3))
    print(mod$finalModel)
    g
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
