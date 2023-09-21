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
                           #Taxon == "Protorhabditis d.l." ~ "Protorhabditis",
                           #Taxon == "Rhabditis d.l." ~ "Rhabditis",
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
  mutate(across(where(is.numeric), ~.*nem.dens$density[match(nem.comp$Unit_quarter, nem.dens$Unit_quarter)])) %>% 
  select(-c("Protorhabditisd.l.",
            "Rhabditisd.l."))


taxa = names(nem.comp)[-1]
# source("https://raw.githubusercontent.com/amynang/marcel/main/R/functions.R")
# nemaplex = query_nemaplex(taxa)
# # calculate st deviation from st error
# nemaplex$StDevMass = nemaplex$StderrMass * sqrt(nemaplex$N)
# # make taxon a column then drop row names
# nemaplex = nemaplex %>% mutate(Taxon = rownames(nemaplex),
#                                .keep ="all",
#                                .before = cp_value)
# rownames(nemaplex) = NULL
# # https://www.wur.nl/en/research-results/chair-groups/plant-sciences/laboratory-of-nematology/nematode-pictures/desmoscolex.htm
# nemaplex$feeding[nemaplex$Taxon == "Desmoscolex"] = "3"
# # http://nemaplex.ucdavis.edu/Taxadata/G916.aspx
# nemaplex$feeding[nemaplex$Taxon == "Prodesmodora"] = "3"
# 
# #replace feeding codes with their meaning
# nemaplex = nemaplex %>% mutate(feeding.type = case_when(feeding == "1" ~ "herbivore",
#                                                         feeding == "2" ~ "fungivore",
#                                                         feeding == "3" ~ "bacterivore",
#                                                         feeding == "5" ~ "predator",
#                                                         feeding == "8" ~ "omnivore"),
#                                 .keep ="all",
#                                 .after = feeding)

# write.csv(nemaplex, "nemaplex.csv")
nemaplex = read.csv("wrangling/nemaplex.csv", row.names = 1, header = TRUE) %>%
  # convert micrograms to mg
  mutate(AvgMass = AvgMass/1e3,
         StDevMass = StDevMass/1e3)

#create an empty-ish dataframe
ecophys = data.frame(Taxon = taxa,
                     AvgMass = NA)

ecophys$AvgMass = nemaplex$AvgMass[match(ecophys$Taxon, nemaplex$Taxon)]

ecophys$feeding.type = nemaplex$feeding.type[match(ecophys$Taxon, nemaplex$Taxon)]

ecophys = ecophys %>% arrange(feeding.type, AvgMass) %>% 
  rename(taxon_name = Taxon) %>% 
  mutate(Abbreviation = case_when(feeding.type == "herbivore"   ~ "Ne-H",
                                  feeding.type == "bacterivore" ~ "Ne-B",
                                  feeding.type == "fungivore"   ~ "Ne-F",
                                  feeding.type == "omnivore"    ~ "Ne-O",
                                  feeding.type == "predator"    ~ "Ne-P"))


#library(taxize)
#use_entrez()
#a7755d3fd5110bfa9c79a7ceddad1b4c1a09
#usethis::edit_r_environ()
# taxonomy = tax_name(ecophys$Taxon, 
#                     get = c("class", 
#                             "order", 
#                             "family"), 
#                     db = "ncbi")
# taxonomy$class[taxonomy$query=="Acrolobus"] = "Chromadorea"
# taxonomy$order[taxonomy$query=="Acrolobus"] = "Rhabditida"
# taxonomy$family[taxonomy$query=="Acrolobus"] = "Cephalobidae"
# 
# write.csv(taxonomy, "ncbi_nematodes.csv", 
#           row.names = F)

taxonomy = read.csv("ncbi_nematodes.csv", header = TRUE)


nem.taxonomy = taxonomy %>% 
  mutate(.after = query,
         taxon_name = query) %>% 
  select(-c(db,query)) %>% 
  mutate(Abbreviation = ecophys$Abbreviation[match(.$taxon_name, ecophys$taxon_name)]) %>% 
  mutate(.before = class,
         group = "Nematodes")











# 
# boom = nem.abun %>% 
#   mutate(not.sieved = soil.1$not.sieved[match(.$Unit_quarter, soil.1$sample.code)], 
#          sieved = ifelse(is.na(not.sieved), "yes", "no"), 
#          gram.soil = soil.1$dry.soil.weight.g[match(.$Unit_quarter, soil.1$sample.code)]
#          )
# summary(lm(total.abundance/gram.soil ~ sieved, data = boom))
# library(glmmTMB)
# summary(glmmTMB(total.abundance ~ sieved + offset(log(gram.soil)), 
#                 family = poisson(), 
#                 data = boom))
# 
# 
# 
# 
# library(brms)
# m = brm(bf(total.abundance ~ offset(log(gram.soil)) + Unit_quarter:sieved),
#           family = poisson(),
#           chains = 3,
#           cores = 3,
#           iter = 3000,
#           backend = "cmdstanr",
#           seed = 321,
#           data = boom)
# pp_check(m)
# summary(m)
# conditional_effects(m)
