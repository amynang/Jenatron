library(openxlsx)
library(tidyverse)


d.1 = read.xlsx("H:\\Jenatron\\PlotList_Jenatron.xlsx",
                sheet = "Tabelle4")[,-(8:9)]

d.1.1 = d.1 %>% mutate(.after = subunit,
                       Uq = paste0(unit,"_",subunit)) %>% 
  
                mutate(.after = soil.history,
                       soil_history = case_when(grepl("B_", soil.history) ~ 0,
                                                grepl("E_", soil.history) ~ 1),
                       plant_history = case_when(plant.history == "+" ~ 1,
                                                 plant.history == "-" ~ 0)) %>% 
                mutate(.after = soil.history,
                       treatment = case_when(soil_history == 0 &
                                               plant_history == 0 ~ "00",
                                             soil_history == 0 &
                                               plant_history == 1 ~ "01",
                                             soil_history == 1 &
                                               plant_history == 0 ~ "10",
                                             soil_history == 1 &
                                               plant_history == 1 ~ "11") %>% as.factor())

d.2 = read.xlsx("H:\\Jenatron\\Jenatron_SP6_macrofauna_data_v1.5.xlsx",
                sheet = "qry_abund")

View(t(d.2[d.2$order == "Coleoptera",]))

d.3 = read.xlsx("H:\\Jenatron\\Jenatron_SP6_macrofauna_data_v1.5.xlsx",
                sheet = "qry_length_width")

d.3.1 = d.3 %>% drop_na(length_mm)

# check for outliers! Geophilomorpha have a really long one

d.3.1[d.3.1$order == "Diptera",] %>% 
  mutate(.after = width_mm,
         FreshMass.mg = 10^(- .309 + ( .997 * log10(length_mm)) + (1.595 * log10(width_mm)))) %>% 
  group_by(order) %>% 
  summarise( MeanMass.mg = mean(FreshMass.mg),
             StDMass.mg = sd(FreshMass.mg))





d.4 = d.2[,5:100] %>% t()
colnames(d.4) = d.2$taxon_name
rownames(d.4) = str_split(rownames(d.4), "_", simplify = T)[,2]

d.5 = d.4 %>% 
  as.data.frame() %>% 
  add_column(.before = "Aphididae",
             Unit_quarter = rownames(d.4)) %>% 
  as_tibble() %>% 
  # removing Hymenoptera and Clitellata (soil cores do not provide accurate density estimates)
  select(-c(Formicidae_large, Formicidae_larva, Formicidae_small,
            Diapriidae, Mymaridae_ad,
            Lumbricidae, Enchytraeidae)) %>% 
  # we estimate densities per mesocosm based on abundances in a core of 15cm diameter
  mutate(across(where(is.numeric), ~ (./(pi*7.5^2))*(pi*25^2)))

d.3.2 = d.3.1 %>% 
  # is width larger than length?
  mutate(weird = length_mm < width_mm) %>% 
  # drop those weird ones
  filter(weird != TRUE) %>% select(-weird) %>% 
  # calculate individual fresh bodymasses
  mutate(FreshMass.mg = case_when(
    class ==     "Chilopoda" ~ 10^(- .549 + (1.416 * log10(length_mm)) + (1.543 * log10(width_mm))),
    class ==     "Diplopoda" ~ 10^(-1.400 + (2.443 * log10(length_mm)) + (0.215 * log10(width_mm))))
  )
  
  
  
  
  
  
  
library(brms)

m = brm(bf(FreshMass.mg ~ 0 + taxon_name + (taxon_name|source_sample),
           sigma ~ 0 + taxon_name,
           alpha ~ 0 + taxon_name 
           ),
        family = skew_normal(),
        chains = 4,
        iter = 4000,
        cores = 4,
        #control = list(adapt_delta = 0.99),
        backend = "cmdstanr",
        seed = 2023,
        data = d.3.2 %>% filter(class %in% c("Chilopoda","Diplopoda")))


m = brm(bf(FreshMass.mg ~ 0 + taxon_name 
                       + (0 + taxon_name|source_sample)
                       + (0 + taxon_name|order),
           sigma ~ 0 + taxon_name
           ),
        family = lognormal(),
        chains = 4,
        iter = 4000,
        cores = 4,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr",
        seed = 2023,
        data = d.3.2 %>% 
          filter(class %in% c("Chilopoda",
                              "Diplopoda")))


m = brm(bf(FreshMass.mg ~ 0 + (1|source_sample)),
        #bf(FreshMass.mg ~ 0 + source_sample),
        family = lognormal(),
        chains = 4,
        iter = 4000,
        cores = 4,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr",
        seed = 2023,
        data = d.3.2 %>% 
          filter(taxon_name == "Polydesmidae_large"))

summary(m)
pp_check(m, ndraws = 100)
ranef(m)
View(ranef(m)[[1]][,,2])

plot(conditional_effects(m, #prob = .99,
                         #method = "posterior_predict",
                         effects = "source_sample:taxon_name", 
                         re_formula = NULL))

m = brm(bf(FreshMass.mg ~ 0 + taxon_name 
           + (0 + taxon_name|source_sample),
           sigma ~ 0 + taxon_name
           ),
        family = lognormal(),
        chains = 4,
        iter = 4000,
        cores = 4,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr",
        seed = 2023,
        data = d.3.2)
