library(readxl)
library(sankeywheel)
library(highcharter)
library(tidyverse)

###input data
df <- read_xlsx('../data/Figure2_data/Figure2C/input_sankey_MX.csv',col_names = F)
names(df) <- c("from","to","weight")

###
highchart()%>%
   hc_add_series(data = df,
                 type = 'sankey',
                 hcaes(from=from,to=to,weight=weight))%>%
   hc_add_theme(hc_theme_google())


