library(tidyverse)
library(readxl)
source("scripts/functions/read_hplc_bouss.R")

pigments <- c("chl_c1_c2", "chl_c3", "peri", "but", "fuco", "neox", "prasi", "viola", "hex", "diad", "allo", "diat", "zea", "lutein", "t_chlb", "dv_chla", "chla", "t_chla")
pigtosum <- c("chl_c1_c2", "chl_c3", "peri", "but", "fuco", "neox", "prasi", "viola", "hex", "diad", "allo", "diat", "zea", "lutein", "t_chlb")
diapig <-  c("peri", "but", "fuco", "hex", "allo", "zea", "t_chlb")

hplc <- read_hplc("data/Darkedge2021_takuvik_pigments_160922.xlsx") %>% select( - chlb, - dv_chlb)

lonlat <- read_excel("data/Darkedge2021_takuvik_pigments_160922.xlsx") %>% janitor::clean_names() %>% 
  select(latitude, longitude) 


hplc[hplc == "LOD"] <- "0"

hplc$date <- as.Date(as.numeric(hplc$date), origin = "1899-12-30")


hplc <- bind_cols(hplc, lonlat)
hplc <- hplc %>% mutate_at(c(4:21, 23, 24), as.numeric)


library(vegan)

AFC <- cca(select(hplc, all_of(pigtosum)))

scores <- data.frame(scores(AFC, choices = c(1,2,3), display = "site"))
data_ca <- bind_cols(hplc, scores)

pigscore <- data.frame(scores(AFC, choices = c(1,2,3), display = "species"))


ggplot(data_ca)+
  geom_point(aes(x = CA1, y = CA2, colour = station_name))+
  geom_segment(aes(x = 0, xend = CA1, y = 0, yend = CA2), data = pigscore)+
  geom_text(aes(x = CA1, y = CA2, label = rownames(pigscore)), data = pigscore)+
  scale_color_viridis_d()+
  xlim(-2,2)+
  ylim(-2,2)

disthplc <- dist(select(data_ca, CA1, CA2, CA3))

plot(hclust(disthplc, method = "ward.D"))

data_ca$group <- as.factor(cutree(hclust(disthplc, method = "ward.D"),  h = 10))

ggplot(data_ca)+
  geom_point(aes(x = CA1, y = CA2, colour = group), size = 1.5)+
  geom_segment(aes(x = 0, xend = CA1, y = 0, yend = CA2), data = pigscore)+
  geom_text(aes(x = CA1, y = CA2, label = rownames(pigscore)), data = pigscore, size = 6)+
  theme_bw(base_size = 18)+
  coord_equal()

ggplot(data_ca)+
  geom_point(aes(x = date, y = -depth, colour = group))

world <- map_data("world")
usa <- map_data("ame")

ggplot(data_ca)+
  geom_point(aes(x = longitude, y = latitude, colour = date))+
  geom_polygon(aes(x =long, y = lat, group= group), data = world, alpha = 0.5)+
  xlim(-100, -50)+
  ylim(60, 92)+
  coord_quickmap()+
  theme_minimal()

liste_pf <- read_excel(paste(here(),"data/liste_pf.xlsx", sep = "/"), 
                       sheet = "Feuil3")

liste_pf <- liste_pf %>% mutate(station_name = str_extract(station, "[0-9]{3}")) %>% 
  select(station_name, ctd_number, depth) %>% 
  mutate(phytofloat = "yes",
         depth = round(depth, 1))

data_ca2 <- left_join(data_ca, liste_pf)
table(data_ca2$phytofloat)


