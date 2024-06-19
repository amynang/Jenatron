library(openxlsx)
library(tidyverse)




################################## Composition #################################

macro.ident = read.xlsx("H:\\Jenatron\\Jenatron_SP6_macrofauna_data_v1.5.xlsx",
                        sheet = "qry_abund")

macro.iden = macro.ident[,5:100] %>% t()
colnames(macro.iden) = macro.ident$taxon_name
rownames(macro.iden) = str_split(rownames(macro.iden), "_", simplify = T)[,2]

mac.tax.dens = macro.iden %>% 
  as.data.frame() %>% 
  add_column(.before = "Aphididae",
             Unit_quarter = rownames(macro.iden)) %>% 
  as_tibble() %>% 
  # removing Hymenoptera and Clitellata (soil cores do not provide accurate density estimates)
  select(-c(Formicidae_large, Formicidae_larva, Formicidae_small,
            Diapriidae, Mymaridae_ad,
            Lumbricidae, Enchytraeidae)) %>% 
  # we estimate densities per mesocosm based on abundances in a core of 15cm diameter
  mutate(across(where(is.numeric), ~ (./(pi*7.5^2))*(pi*25^2)))

mac.taxonomy = macro.ident[,1:4] %>% 
  filter(taxon_name %in% names(mac.tax.dens)) %>% 
  mutate(Abbreviation = c(# Abbreviations used in Potapov et al 2022
                          "He-Ste","Dpod","Dpod","Dplu-D","Dipt-M","He-Ste",
                          "Cpt-H-H","Cpt-O","Ga-Sn","Chi-Ge","Dpod","Ar-La-We",
                          "Ar-Sm-We","Chi-Li","Chi-Li","Dipt","Dpod-L",
                          "Dpod-L","Chi-Ge","Cpt-P-Sta","Ga-Sn")) %>% 
  arrange(class,order,family) %>% 
  mutate(.before = class,
         group = "macrofauna")


################################### Body size ##################################

mac.length = read.xlsx("H:\\Jenatron\\Jenatron_SP6_macrofauna_data_v1.5.xlsx",
                sheet = "qry_length_width") %>% 
              drop_na(length_mm)

# check for outliers! Geophilomorpha have a really long one

mac.lengths = mac.length %>% 
  # is width larger than length?
  mutate(weird = length_mm < width_mm) %>% 
  # drop those weird ones
  filter(weird != TRUE) %>% select(-weird) %>%
  mutate(.before = ind_ID,
         .keep = "unused",
         Unit_quarter = str_split(source_sample, "_", simplify = T)[,2]) %>% 
  select(Unit_quarter, class, order, family, taxon_name, length_mm, width_mm) %>% 
  # removing Hymenoptera and Clitellata (soil cores do not provide accurate density estimates)
  filter(!(taxon_name %in% c("Formicidae_large", "Formicidae_larva", "Formicidae_small",
                           "Diapriidae", "Mymaridae_ad",
                           "Lumbricidae", "Enchytraeidae")))

mac.masses = mac.lengths %>% 
  # calculate individual fresh bodymasses
  mutate(FreshMass.mg = case_when(# based on Sohlström 2018 10.1002/ece3.4702
                                  # group specific coefficients (model 1)
                                  order ==            "Araneae" ~ 10^(- .281 + (1.368 * log10(length_mm)) + (1.480 * log10(width_mm))),
                                  order ==     "Lithobiomorpha" ~ 10^(- .549 + (1.416 * log10(length_mm)) + (1.543 * log10(width_mm))),
                                  order ==     "Geophilomorpha" ~ 10^(- .419 + ( .964 * log10(length_mm)) + (1.766 * log10(width_mm))),
                                  class ==          "Diplopoda" ~ 10^(-1.400 + (2.443 * log10(length_mm)) + ( .215 * log10(width_mm))),
                                  order ==            "Diptera" ~ 10^(- .309 + ( .997 * log10(length_mm)) + (1.595 * log10(width_mm))),
                                  order ==         "Coleoptera" ~ 10^(- .286 + ( .840 * log10(length_mm)) + (1.954 * log10(width_mm))),
                                  order ==          "Hemiptera" ~ 10^(- .420 + (1.177 * log10(length_mm)) + (1.431 * log10(width_mm))),
                                  # The general relationship (model 3)
                                  TRUE ~ 10^(- .285 + (1.040 * log10(length_mm)) + (1.585 * log10(width_mm))))) %>% 
  select(-c(length_mm,width_mm))
  