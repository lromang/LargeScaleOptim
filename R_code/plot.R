library(ggplot2)


## Read in dataset
ndata        <- read.csv("../tests/lbfgsT02.csv", header = FALSE)
names(ndata) <- c("data_points", "precision")

## Additional variables LBFGS
ndata$T <- rep("02", nrow(ndata))


## Additional variables TN
data$sampleSize  <- rep(.05, nrow(data))
data$CG          <- rep(5, nrow(data))
data$SGD_iters   <- rep(100, nrow(data))
data$SGD_batch   <- rep(1000, nrow(data))

all_data <- rbind(all_data, ndata)

## Plot
ggplot(data = ndata,
       aes(x = data_points,
           y = precision,
           col = as.factor(T))) +
    geom_line() +
    theme(panel.background = element_blank())
