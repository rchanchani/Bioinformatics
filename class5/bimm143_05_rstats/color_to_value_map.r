map.colors <- function (x,
                        high.low = range(x),
                        palette = cm.colors(100)) {
  ## Description: Map the values of the input vector 'x'
  ##  to the input colors vector 'palette'
  
  # percent value of high_low range
  percent <- ((value-high.low[2]) / (high.low[1] - high.low[2]))
  # find index position and if value is 0, adjust to 1%
  index <- round((length(palette) - 1) * percent) + 1
  return (palette[index])
}