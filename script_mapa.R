install.packages('units')
install.packages('sf')

require(sf)
require(tidyverse)
require(RColorBrewer)

setwd("~/COVID-19")

regionais <- st_read("ERS_SAUDE_MT_FINAL.shp")
ers <- read.csv2("ers_aglomerados.csv",h=T) %>% 
  mutate(Geocodigo = as.character(Geocodigo))

regionais <- regionais %>% left_join(ers, by = "Geocodigo")

colors <- c("#0739ed","#f25cf2","#FFFFB8","#CECFD1",
            "#E5258C","#884E88","#c9f74a","#FDD22B",
            "#8a7460","#139A4B","#F8AC92","#FF8622",
            "#f52040","#C5B76C","#9ADEF5","#FEF104")

#Mapa por regionais
ggplot(regionais) +
  geom_sf(aes(fill = as.factor(Geocodigo))) +
  theme_minimal() + 
  scale_fill_manual(values = colors,
                    guide = guide_legend(title = "Regionais",reverse = TRUE))

#Mapa por variavel selecionada (r_aux_025)
ggplot(regionais) +
  geom_sf(aes(fill = r_aux_025)) +
  theme_minimal() + 
  scale_fill_gradientn(colours=rev(brewer.pal(9, "YlOrRd")),na.value = "#ffffff") 
 