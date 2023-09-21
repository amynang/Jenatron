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
  
# write.csv(mac.masses, "macromasses.csv")  
 
# mac.mean.masses = mac.masses %>%
#   group_by(Unit_quarter, taxon_name) %>%
#   summarise(N = n(),
#             mean.mass.mg = mean(FreshMass.mg))
#   
#   
#   
# library(brms)
# 
# #### partial pooling ####
# mpp = brm(bf(FreshMass.mg ~ 0 + taxon_name# + (0 + taxon_name|Unit_quarter)
#              #,sigma ~ 0 + taxon_name
#              ),
#           family = Gamma(link = "log"),
#           chains = 4,
#           cores = 4,
#           iter = 4000,
#           backend = "cmdstanr",
#           seed = 123,
#           data = mac.masses)
# 
# 
# 
# #### partial pooling within groups, group means varying by unit ####
# mppgp = brm(bf(FreshMass.mg ~ 0 + (1|gr(taxon_name, by = order)) + (0 + taxon_name|Unit_quarter)
#                #,sigma ~ 0 + (1|gr(taxon_name, by = order))
#                ),
#             family = lognormal(),
#             chains = 4,
#             cores = 4,
#             iter = 4000,
#             control = list(adapt_delta = 0.95),
#             backend = "cmdstanr",
#             seed = 321,
#             data = mac.masses)
# 
# summary(mppgp)
# 
# # "pop-level" mu and sigma for each taxon, 
# # followed by unit differences from pop-level mu
# ranef(mppgp) 
# plot(conditional_effects(mppgp,
#                          effects = "u:t", 
#                          re_formula = NULL))
# ce = conditional_effects(mppgp,
#                          effects = "u:t", 
#                          re_formula = NULL)
# ce[[1]] # estimate__ gives unit specific estimates for mean body-mass
# 
# 
# 
# m = brm(bf(FreshMass.mg ~ 0 + taxon_name + (taxon_name|source_sample),
#            sigma ~ 0 + taxon_name,
#            alpha ~ 0 + taxon_name 
#            ),
#         family = skew_normal(),
#         chains = 4,
#         iter = 4000,
#         cores = 4,
#         #control = list(adapt_delta = 0.99),
#         backend = "cmdstanr",
#         seed = 2023,
#         data = d.3.2 %>% filter(class %in% c("Chilopoda","Diplopoda")))
# 
# 
# m = brm(bf(FreshMass.mg ~ 0 + taxon_name 
#                        + (0 + taxon_name|source_sample)
#                        + (0 + taxon_name|order),
#            sigma ~ 0 + taxon_name
#            ),
#         family = lognormal(),
#         chains = 4,
#         iter = 4000,
#         cores = 4,
#         control = list(adapt_delta = 0.99),
#         backend = "cmdstanr",
#         seed = 2023,
#         data = d.3.2 %>% 
#           filter(class %in% c("Chilopoda",
#                               "Diplopoda")))
# 
# 
# m = brm(bf(FreshMass.mg ~ 0 + (1|source_sample)),
#         #bf(FreshMass.mg ~ 0 + source_sample),
#         family = lognormal(),
#         chains = 4,
#         iter = 4000,
#         cores = 4,
#         control = list(adapt_delta = 0.99),
#         backend = "cmdstanr",
#         seed = 2023,
#         data = d.3.2 %>% 
#           filter(taxon_name == "Polydesmidae_large"))
# 
# summary(m)
# pp_check(m, ndraws = 100)
# ranef(m)
# View(ranef(m)[[1]][,,2])
# 
# plot(conditional_effects(m, #prob = .99,
#                          #method = "posterior_predict",
#                          effects = "source_sample:taxon_name", 
#                          re_formula = NULL))
# 
# m = brm(bf(FreshMass.mg ~ 0 + taxon_name 
#            + (0 + taxon_name|source_sample),
#            sigma ~ 0 + taxon_name
#            ),
#         family = lognormal(),
#         chains = 4,
#         iter = 4000,
#         cores = 4,
#         control = list(adapt_delta = 0.99),
#         backend = "cmdstanr",
#         seed = 2023,
#         data = d.3.2)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# mat = as.data.frame(t(structure(c(
# 1.45,	2.05	,1.17,
# 3.7	,  2.84	,1.56,
# 1.3	,  1.97	,1.07,
# 3.73,	3.71	,1.63,
# 4.18,	3.95	,1.68,
# 1.79,	2.97	,1.26,
# 4.05,	3.63	,1.51,
# 2.78,	3.03	,1.26,
# 1.88,	2.77	,1.14,
# 5.67,	4.34	,1.78,
# 4.8	,  4.26	,1.74,
# 3.08,	3.19	,1.29,
# 4.9	,  4.09	,1.63,
# 1.45,	2.59	,1.03,
# 1.82,	2.98	,1.18,
# 3.78,	3.5	  ,1.34,
# 3.6	,  3.35	,1.26,
# 6.28,	4.68	,1.75,
# 2.61,	3.53	,1.31,
# 2.59,	3.07	,1.13,
# 4.58,	4.1	  ,1.49,
# 2.84,	3.9	  ,1.4,
# 5.47,	4.57	,1.63,
# 4.22,	3.99	,1.42,
# 7.52,	4.92	,1.75,
# 1.47,	2.91	,1.02,
# 4.69,	4.27	,1.47,
# 0.84,	2.68	,0.92,
# 0.99,	2.66	,0.9,
# 0.73,	2.42	,0.81,
# 4.75,	4.27	,1.41,
# 0.93,	2.7	  ,0.89,
# 1.1	,  2.93	,0.96,
# 1.06,	3.09	,1.01,
# 5.05,	4.59	,1.49,
# 1.24,	2.67	,0.86,
# 0.88,	2.47	,0.79,
# 1.04,	2.44	,0.78,
# 1.08,	3.05	,0.97,
# 1.14,	3.18	,1.01,
# 4.58,	4.05	,1.28,
# 1.51,	3.42	,1.08,
# 1.08,	2.93	,0.9,
# 0.93,	2.81	,0.86,
# 0.83,	2.29	,0.7,
# 0.68,	2.49	,0.76,
# 6.86,	5.06	,1.54,
# 1	  ,  2.68	,0.81,
# 0.93,	2.61	,0.78,
# 1.31,	3.18	,0.91,
# 0.96,	2.68	,0.76,
# 0.9	,   2.6	,0.72,
# 0.93,	2.77	,0.72,
# 0.82,	2.55	,0.65), .Dim = c(3L,54L))))
# colnames(mat) = c("m","l","w")
# 
# log10(m) ~ log10(l) + log10(w)
# m = brm(bf(log10(m) ~ log10(l) + log10(w)),
#         family = gaussian(),
#         chains = 4,
#         iter = 4000,
#         cores = 4,
#         control = list(adapt_delta = 0.9),
#         backend = "cmdstanr",
#         seed = 2023,
#         data = mat)
# summary(m)
