PointsWithLine <- function(data, mapping, ..., method = "pearson") {
  if ("colour" %in% colnames(mapping)) {
    df <- data.frame(x = data[,eval(as.character(mapping$x)[2], data)],
                     y = data[,eval(as.character(mapping$y)[2], data)],
                     c = data[,eval(as.character(mapping$colour)[2], data)])
  }else{
    df <- data.frame(x = data[,eval(as.character(mapping$x)[2], data)],
                     y = data[,eval(as.character(mapping$y)[2], data)])
  }
  
  xMax <- max(df$x, na.rm = TRUE)
  yMax <- max(df$y, na.rm = TRUE)
  
  if ("colour" %in% colnames(mapping)) {
    graph = ggplot(data, mapping) + geom_point(aes(shape = df$c), stroke = 0.7) +
      scale_shape_manual(values = 1:12) +
      geom_smooth(
        method = "lm", se=FALSE ,
        aes(x = data[,eval(as.character(mapping$x)[2], data)],
            y = data[,eval(as.character(mapping$y)[2], data)], color = NA)
      )
  }else{
    graph = ggplot(data, mapping) + geom_point(stroke = 0, alpha = 0.01, color = "#BB933E") +
      geom_smooth(
        method = "lm", se=FALSE, color = "#255B87"
      )
  }
}