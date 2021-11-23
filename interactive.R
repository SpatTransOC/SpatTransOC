
SpatialDimPlot(e[[1]], interactive = TRUE)
SpatialDimPlot(e[[2]], interactive = TRUE)

SpatialDimPlot(p1, interactive = TRUE)


# subset 030711

outliers <- list()

p1 <- ggplot(GetTissueCoordinates(p[[6]]),
       aes(x = imagecol, y = imagerow)) +
    geom_point() +
    coord_equal()

outliers[[1]] <- CellSelector(p1)
outliers[[1]] <- c(outliers[[1]], CellSelector(p1))

outliers[[2]] <- CellSelector(p1)
outliers[[2]] <- c(outliers[[2]], CellSelector(p1))

outliers[[3]] <- CellSelector(p1)
outliers[[3]] <- c(outliers[[3]], CellSelector(p1))

outliers[[4]] <- CellSelector(p1)
outliers[[4]] <- c(outliers[[4]], CellSelector(p1))

outliers[[5]] <- CellSelector(p1)
outliers[[5]] <- c(outliers[[5]], CellSelector(p1))

outliers[[6]] <- CellSelector(p1)
outliers[[6]] <- c(outliers[[6]], CellSelector(p1))

outliers[[7]] <- CellSelector(p1)
outliers[[7]] <- c(outliers[[7]], CellSelector(p1))

outliers[[8]] <- CellSelector(p1)
outliers[[8]] <- c(outliers[[8]], CellSelector(p1))

outliers[[9]] <- CellSelector(p1)
outliers[[9]] <- c(outliers[[9]], CellSelector(p1))

outliers[[10]] <- CellSelector(p1)
outliers[[10]] <- c(outliers[[10]], CellSelector(p1))

outliers[[11]] <- CellSelector(p1)
outliers[[11]] <- c(outliers[[11]], CellSelector(p1))

outliers[[12]] <- CellSelector(p1)
outliers[[12]] <- c(outliers[[12]], CellSelector(p1))




SpatialDimPlot(a)

library(plotly)

ggplotly(p1)

plotly::plotly_data(ggplotly(p1))

