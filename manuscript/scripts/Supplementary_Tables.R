library(kableExtra)
dt <- mtcars[1:5, 1:6]

print(kbl(dt, booktabs = T, ) %>%
  kable_styling(latex_options = "striped"))


