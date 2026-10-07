litedown::reactor(warning = FALSE)

library(easybio)
library(data.table)

# Example data: airquality, focusing on May and June
air <- subset(airquality, Month %in% c(5, 6))
setDT(air)

a <- Artist$new(data = air)

a$plot_scatter(x = Wind, y = Temp)

a$get_all_result()[, .(command)]

a$plot_box(x = factor(Month), y = Ozone)
a$test_wilcox(formula = Ozone ~ Month)
a$plot_box(x = factor(Month), y = Ozone)

a$plot_scatter(x = Wind, y = Temp)
a$test_t(formula = Temp ~ Month)
a$plot_scatter(x = Wind, y = Temp)

a$plot_scatter(x = Wind, y = Temp)$
  test_wilcox(formula = Ozone ~ Month)$
  plot_scatter(x = Wind, y = Temp)

a$plot_scatter(
  fun = \(x) x[, z := Wind * Temp],
  x = Wind,
  y = z
)

# Summarise data for dumbbell
air[, mean_temp := mean(Temp), by = Month]
a$plot_dumbbell(x = Month, y = mean_temp, col = factor(Month))

air[, .N, by = .(Month, Day)]
a2 <- Artist$new(data = air[, .N, by = .(Month, Day)])
a2$plot_bubble(x = Month, y = Day, size = N, col = N)

# Create some example data with positive and negative values
example <- data.table(
  group = letters[1:6],
  value = c(1.5, -0.8, 2.1, -1.2, 0.5, -0.3)
)
a3 <- Artist$new(data = example)
a3$plot_barchart_divergence(group = group, y = value)

a3$plot_lollipop(x = group, y = value)

a$plot_contour(x = Wind, y = Temp)

a$plot_scatter_ellipses(x = Wind, y = Temp, col = factor(Month))

month_counts <- air[, .N, by = Month]
a4 <- Artist$new(data = month_counts)
a4$plot_donut(x = Month, y = N, fill = factor(Month))

a4$plot_pie(y = N, fill = factor(Month))

history <- a$get_all_result()
history

# The last plot
a$result[[length(a$result)]]

# The last statistical test result
for (res in a$result) {
  if (inherits(res, "htest")) print(res)
}

