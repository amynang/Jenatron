library(openxlsx)
library(tidyverse)

potapov = read.xlsx("H:/Literature/brv12857-sup-0003-tables2.xlsx",
                    sheet = "GroupList") %>% 
  rename(any_of(setNames(c("P","A","L","W","S","B","F","Fa"),
                         c("Plants","Algae","Litter","deadWood","Detritus","Bacteria","Fungi","Animals"))))

nem.mes.mac.taxonomy = full_join(nem.taxonomy, mes.taxonomy) %>% 
  full_join(., mac.taxonomy) %>% 
  mutate(Agility = potapov$Agility[match(.$Abbreviation, potapov$Abbreviation)],
         PhysicalProtection = potapov$PhysicalProtection[match(.$Abbreviation, potapov$Abbreviation)],
         Metabolites = potapov$Metabolites[match(.$Abbreviation, potapov$Abbreviation)],
         above = potapov$above[match(.$Abbreviation, potapov$Abbreviation)],
         epi = potapov$epi[match(.$Abbreviation, potapov$Abbreviation)],
         hemi = potapov$hemi[match(.$Abbreviation, potapov$Abbreviation)],
         eu = potapov$eu[match(.$Abbreviation, potapov$Abbreviation)])

nematodes = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$group == "Nematodes"]
meso = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$group == "mesofauna"]
macro = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$group == "macrofauna"]
# indexing for the metabolic loss regressions
insecta      = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$class == "Insecta"]
araneae      = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$order == "Araneae"]
chilopoda    = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$class == "Chilopoda"]
progoneata   = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$class %in% c("Pauropoda","Symphyla","Diplopoda")]
mesostigmata = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$order == "Mesostigmata"]
oribatida    = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$order == "Sarcoptiformes"]
prostigmata  = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$order == "Trombidiformes"]


nem.mes.mac.dens = full_join(nem.tax.dens, mes.tax.dens, by = join_by(Unit_quarter)) %>% 
  full_join(., mac.tax.dens, by = join_by(Unit_quarter))


masses = full_join(mes.masses, mac.masses)
mes.mac.mean.masses = masses %>% 
  group_by(Unit_quarter, taxon_name) %>% 
  summarise(Bodymass.mg = mean(FreshMass.mg)) %>% 
  ungroup()
across.unit.avg = mes.mac.mean.masses %>% 
  group_by(taxon_name) %>% 
  summarise(Bodymass.mg = mean(Bodymass.mg)) %>% 
  ungroup()
# split density dataframe to mesocosm-specific dataframes
att = nem.mes.mac.dens %>% 
  add_column(.before = "Achromadora",
             plants = 1, # dummy values so that basals are part of the foodweb
             detritus = 1,
             microbes = 1) %>% 
  pivot_longer(2:length(names(.)),
               names_to = "taxon_name",
               values_to = "density")%>% 
  filter(density>0) %>% 
  # finally, create a list: each element is a mesocosm
  split(., with(.,Unit_quarter))

for (i in 1:length(att)) {
  att[[i]] = att[[i]] %>% 
    mutate(.after = taxon_name,
           Agility = nem.mes.mac.taxonomy$Agility[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
           PhysicalProtection = nem.mes.mac.taxonomy$PhysicalProtection[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
           Metabolites = nem.mes.mac.taxonomy$Metabolites[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
           above = nem.mes.mac.taxonomy$above[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
           epi = nem.mes.mac.taxonomy$epi[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
           hemi = nem.mes.mac.taxonomy$hemi[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
           eu = nem.mes.mac.taxonomy$eu[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)]) %>% 
    mutate(Bodymass.mg = case_when(.$taxon_name %in% nematodes ~ ecophys$AvgMass[match(.$taxon_name, ecophys$taxon_name)],
                                   .$taxon_name %in% c(meso,macro) ~ mes.mac.mean.masses$Bodymass.mg[match(paste(.$Unit_quarter, .$taxon_name), 
                                                                                                           paste(mes.mac.mean.masses$Unit_quarter, mes.mac.mean.masses$taxon_name))],
                                   TRUE ~ 1)) %>% 
    # mean imputation for taxa  that were counted in a mesocosm but not measured (11 cases)
    mutate(Bodymass.mg = ifelse(is.na(Bodymass.mg), across.unit.avg$Bodymass.mg[match(.$taxon_name, 
                                                                                      across.unit.avg$taxon_name)], 
                                Bodymass.mg)) %>% 
    mutate(ind.Loss.Jh = case_when(# based on Ehnes 2011 10.1111/j.1461-0248.2011.01660.x
             # group specific coefficients (phylogenetic model)
             .$taxon_name %in% araneae ~ exp(24.581475 + .5652537*log(Bodymass.mg) - .7093476*(1/(8.62*1e-5*(20+273.15)))),
             .$taxon_name %in% insecta ~ exp(21.972050 + .7588950*log(Bodymass.mg) - .6574038*(1/(8.62*1e-5*(20+273.15)))),
             .$taxon_name %in% chilopoda ~ exp(28.252911 + .5580991*log(Bodymass.mg) - .8030069*(1/(8.62*1e-5*(20+273.15)))),
             .$taxon_name %in% progoneata ~ exp(22.347024 + .5713411*log(Bodymass.mg) - .6700449*(1/(8.62*1e-5*(20+273.15)))),
             .$taxon_name %in% oribatida ~ exp(22.022770 + .6793706*log(Bodymass.mg) - .7060855*(1/(8.62*1e-5*(20+273.15)))),
             .$taxon_name %in% prostigmata ~ exp(10.281495 + .6599399*log(Bodymass.mg) - .4125318*(1/(8.62*1e-5*(20+273.15)))),
             .$taxon_name %in% mesostigmata ~ exp(9.6740230 + .6904864*log(Bodymass.mg) - .3792541*(1/(8.62*1e-5*(20+273.15)))),
             # The general relationship (linear model)
             TRUE ~ exp(23.055335 + .6950710*log(Bodymass.mg) - .6864200*(1/(8.62*1e-5*(20+273.15))))),
           Biomass.mg = density*Bodymass.mg,
           pop.Loss.Jh = if_else(.$taxon_name %in% c("plants","detritus","microbes"), 0, # not necessary, basal losses are ignored
                                 density*ind.Loss.Jh)) %>% 
    #https://www.ecologycenter.us/species-richness/the-importance-of-transfer-efficiencies-in-determining-energy-pathways.html
    mutate(efficiency = case_when(# based on Lang 2017 10.1111/oik.04419
                                  .$taxon_name == "detritus" ~ .158,
                                  .$taxon_name == "plants" ~ .545,
                                  #.$taxon_name == "microbes" ~ .545, # not sure this makes sense, 
                                                                      # if complex structural chemicals 
                                                                      # such as cellulose and lignin lowers assim. ef. of two above,
                                                                      # then this should not apply to microbes
                                  TRUE ~ .906))
}

# View(att[[1]])
# 
# long = do.call(rbind.data.frame, att)
