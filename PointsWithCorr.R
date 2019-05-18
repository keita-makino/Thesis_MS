PointsWithCorr <- function(data, mapping, ..., method = "pearson") {
  if ("colour" %in% colnames(mapping)) {
    df <- data.frame(x = data[,eval(as.character(mapping$x)[2], data)],
                     y = data[,eval(as.character(mapping$y)[2], data)],
                     c = data[,eval(as.character(mapping$colour)[2], data)])
  }else{
    df <- data.frame(x = data[,eval(as.character(mapping$x)[2], data)],
                     y = data[,eval(as.character(mapping$y)[2], data)])
  }
  
  xPos = min(df$x)
  yPos = max(df$y, (lm(df$y~df$x) %>% .$coef %>% .[1]) + (lm(df$y~df$x) %>% .$coef %>% .[2]) * max(df$x))

  if ("colour" %in% colnames(mapping)) {
    sumdf <- df %>%
      group_by(c) %>%
      summarise(
        lab = round(cor(x, y),3),
        x = xPos,
        y = yPos*min(as.numeric(c))/max(as.numeric(df$c))
      )
    
    ggplot(data, mapping) + geom_point(stroke = 0, alpha = 0.5, color = "#BB933E") +
      ggplot2::geom_label(
        data = sumdf,
        mapping = ggplot2::aes(x = x, y = y, label = lab, color = c),
        hjust = 0, vjust = 1,
        size = 3,
        inherit.aes = FALSE # do not inherit anything from the ...
      )
  }else{
    sumdf <- df %>%
      summarise(
        lab = round(cor(x, y),3),
        x = xPos,
        y = yPos,
        col = "a"
      )
    
    if (sumdf$lab %>% abs < 0.15) {
      sumdf$col = "#bbbbbb"
    } else if(sumdf$lab < -0.5){
      sumdf$col = "#540829"
    }else if(sumdf$lab < -0.15){
      sumdf$col = "#CC5287"
    }else if(sumdf$lab < 0.5){
      sumdf$col = "#255B87"
    }else{
      sumdf$col = "#072B48"
    }
    
    graph = ggplot(data, mapping) + geom_point(shape = 21, alpha = 0.01, fill = "#BB933E", color = "#BB933E") +
      ggplot2::geom_label(
        data = sumdf,
        mapping = ggplot2::aes(x = x, y = y, label = lab),
        hjust = 0, vjust = 1,
        size = 3,
        inherit.aes = FALSE,
        fontface = ifelse(sumdf$lab %>% abs < 0.5,"plain","bold.italic"),
        alpha = 0.75,
        color = sumdf$col # do not inherit anything from the ...
      ) +
      geom_smooth(
        method = "lm", se=FALSE, color = "#255B87"
      )
  }
}