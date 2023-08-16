library(openxlsx)
library(tidyverse)


soil.1 = read.xlsx(# soil info from nematode extraction
                   "H:\\Jenatron\\Angelos_Jenatron_Nematode Extraction_ 10_22.xlsx",
                   startRow = 2) %>% 
  # brackets cause trouble
  rename(any_of(setNames(c("weight.pot.[g]","fresh.soil.weight.[g]","dry.weight.[g]","X6"),
                         c("weight.pot.g","fresh.soil.weight.g","dry.weight.g","dry.soil.weight.g" )))) %>% 
  mutate(.after = "dry.soil.weight.g",
         # amount of evaporated water
         water.g = fresh.soil.weight.g - dry.soil.weight.g,
         # proportion of water
         water.prop = water.g/fresh.soil.weight.g)


soil.2 = read.xlsx(# soil info from soil density cores
                   "H:\\Jenatron\\Angelos_Jenatron_SoilDensity.xlsx",
                   sheet = "Sheet1") %>% 
  mutate(# dry soil weight in a core of 5cm diameter 10cm length
         soil.dry.weight.g = dry.bagnsoil.g - bag.weight.g)

nem.abun = read.xlsx(# extracted nematodes
                     "H:\\Jenatron\\Nematode genera Jenatron SP6 nematodes 2022 Amyntas KIM corrected.xlsx",
                     sheet = "abundance",
                     startRow = 3) %>% 
  rename(Unit_quarter = X1) %>% 
  mutate(Unit_quarter = str_c("U",
                              str_split(.$Unit_quarter, 
                                        "((?<=[a-z])(?![a-z])|(?<=\\d)(?!\\d))", 
                                        simplify = T)[,1],
                              "q",
                              str_split(.$Unit_quarter, 
                                        "((?<=[a-z])(?![a-z])|(?<=\\d)(?!\\d))", 
                                        simplify = T)[,2]))

nem.dens = nem.abun %>% 
  mutate(# nematode density per gram of soil
         dens.per.g = total.abundance / soil.1$dry.soil.weight.g[match(nem.abun$Unit_quarter, 
                                                                       soil.1$sample.code)],
         # g of soil in a 0-10cm "pie"
         soil.g.per.mesocosm = (soil.2$soil.dry.weight.g[match(nem.abun$Unit_quarter, soil.2$Unit_quarter)]*
                               pi*25^2)/(pi*2.5^2),
         # nematode density in the pie
         density = dens.per.g*soil.g.per.mesocosm)


nem.com = read.xlsx(# identified nematodes
                    "H:\\Jenatron\\Nematode genera Jenatron SP6 nematodes 2022 Amyntas KIM corrected.xlsx",
                     sheet = "identified nematodes",
                     rows = c(6,8:118),
                     cols = 2:290)
nem.comp = t(nem.com[2:289]) 
colnames(nem.comp) = nem.com$X1
nem.comp[is.na(nem.comp)] = 0

nem.comp = nem.comp %>% as.data.frame() %>% 
  add_column(.before = "Achromadora",
             Unit_quarter = names(nem.com)[-1]) %>% 
  as_tibble() %>% 
  pivot_longer(2:112,
               names_to = "Taxon",
               values_to = "Abundance")%>% 
  mutate(Unit_quarter = str_replace(Unit_quarter, "u", "U"),
         Taxon = case_when(Taxon == "Chrysonemoides" ~ "Chrysonema", # Probably a synonym of Chrysonema Thorne, 1929 (http://nemaplex.ucdavis.edu/Taxadata/G695.aspx)
                           Taxon == "Diphterophora" ~ "Diphtherophora",
                           Taxon == "Macroposthonia/Mesocriconema" ~ "Mesocriconema",
                           Taxon == "Nigolaimus" ~ "Nygolaimus",
                           Taxon == "Paramphidellus" ~ "Paramphidelus",
                           Taxon == "Paratrophorus" ~ "Paratrophurus",
                           Taxon == "Protorhabditis d.l." ~ "Protorhabditis",
                           Taxon == "Rhabditis d.l." ~ "Rhabditis",
                           .default = Taxon),
         Taxon = str_replace(Taxon, " ", "")) %>% 
  group_by(Unit_quarter, Taxon) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Taxon,
              values_from = Abundance) %>% 
  mutate(Total = rowSums(.[2:110])) %>% 
  mutate(across(where(is.numeric), ~ ./Total)) %>% 
  select(-Total)

# density of individual taxa based on total nematode density and proportions 
# of taxa in identified nematodes
nem.tax.dens = nem.comp %>% 
  mutate(across(where(is.numeric), ~.*nem.dens$density[match(nem.comp$Unit_quarter, nem.dens$Unit_quarter)]))


taxa = names(nem.comp)[-1]
source("https://raw.githubusercontent.com/amynang/marcel/main/R/functions.R")
nemaplex = query_nemaplex(taxa)
# calculate st deviation from st error
nemaplex$StDevMass = nemaplex$StderrMass * sqrt(nemaplex$N)
# make taxon a column then drop row names
nemaplex = nemaplex %>% mutate(Taxon = rownames(nemaplex),
                               .keep ="all", 
                               .before = cp_value)
rownames(nemaplex) = NULL
# https://www.wur.nl/en/research-results/chair-groups/plant-sciences/laboratory-of-nematology/nematode-pictures/desmoscolex.htm
nemaplex$feeding[nemaplex$Taxon == "Desmoscolex"] = "3"
# http://nemaplex.ucdavis.edu/Taxadata/G916.aspx
nemaplex$feeding[nemaplex$Taxon == "Prodesmodora"] = "3"

#replace feeding codes with their meaning
nemaplex = nemaplex %>% mutate(feeding.type = case_when(feeding == "1" ~ "herbivore",
                                                        feeding == "2" ~ "fungivore",
                                                        feeding == "3" ~ "bacterivore",
                                                        feeding == "5" ~ "predator",
                                                        feeding == "8" ~ "omnivore"),
                                .keep ="all", 
                                .after = feeding)

# write.csv(nemaplex, "nemaplex.csv")
nemaplex = read.csv("nemaplex.csv", row.names = 1, header = TRUE)

#bodymass data from 10.1890/11-0546.1
muldervonk = read.csv("https://raw.githubusercontent.com/amynang/MulderVonk2011/main/Mulder%26Vonk2011_bodymass%26feeding.csv",
                      sep = ";", dec = ",")

#create an empty-ish dataframe
ecophys = data.frame(Taxon = taxa,
                     AvgMass = NA)
# we will rely on Mulder & Vonk (2011) for bodymass information
# if a taxon is not there we resort to nemaplex
ecophys$AvgMass = ifelse(ecophys$Taxon %in% muldervonk$TAX.MORPHON, 
                         muldervonk$AvgMass[match(ecophys$Taxon, muldervonk$TAX.MORPHON)], 
                         nemaplex$AvgMass[match(ecophys$Taxon, nemaplex$Taxon)])
ecophys$StDevMass = ifelse(ecophys$Taxon %in% muldervonk$TAX.MORPHON, 
                           muldervonk$StDevMass[match(ecophys$Taxon, muldervonk$TAX.MORPHON)], 
                           nemaplex$StDevMass[match(ecophys$Taxon, nemaplex$Taxon)])
# ... except when Mulder & Vonk have measured only one individual
ecophys$AvgMass[ecophys$Taxon == "Boleodorus"] = nemaplex$AvgMass[nemaplex$Taxon == "Boleodorus"]
ecophys$StDevMass[ecophys$Taxon == "Boleodorus"] = nemaplex$StDevMass[nemaplex$Taxon == "Boleodorus"]
ecophys$AvgMass[ecophys$Taxon == "Discolaimus"] = nemaplex$AvgMass[nemaplex$Taxon == "Discolaimus"]
ecophys$StDevMass[ecophys$Taxon == "Discolaimus"] = nemaplex$StDevMass[nemaplex$Taxon == "Discolaimus"]
ecophys$AvgMass[ecophys$Taxon == "Dorylaimellus"] = nemaplex$AvgMass[nemaplex$Taxon == "Dorylaimellus"]
ecophys$StDevMass[ecophys$Taxon == "Dorylaimellus"] = nemaplex$StDevMass[nemaplex$Taxon == "Dorylaimellus"]
ecophys$AvgMass[ecophys$Taxon == "Prionchulus"] = nemaplex$AvgMass[nemaplex$Taxon == "Prionchulus"]
ecophys$StDevMass[ecophys$Taxon == "Prionchulus"] = nemaplex$StDevMass[nemaplex$Taxon == "Prionchulus"]

ecophys$feeding.type = nemaplex$feeding.type[match(ecophys$Taxon, nemaplex$Taxon)]

ecophys = ecophys %>% arrange(feeding.type, AvgMass)

# One taxon has 0 sd!!!
