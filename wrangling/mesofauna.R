library(openxlsx)
library(tidyverse)

################################## Composition #################################

# counted mesofauna (in core abundances)
mes.counts = read.xlsx("H:\\Jenatron\\Angelos_Jenatron_mesofauna_in_Core.xlsx",
                       sheet = "counts")
# identified Collembola & Acari
mes.ident = read.xlsx("H:\\Jenatron\\Jenatron_SP6_mesofauna_data_v3.10.xlsx",
                sheet = "qry_abund_mesof")
# correct name
mes.ident$taxon_name[mes.ident$taxon_name == "Blattisociidae_x000D_\nBlattisociidae"] = "Blattisociidae"
mes.ident$taxon_name[mes.ident$taxon_name == "Eosentomidae"] = "Protura"
mes.ident$taxon_name[mes.ident$taxon_name == "Pauropodidae"] = "Pauropoda"

mes.iden = mes.ident[,5:244] %>% t()
mes.iden[is.na(mes.iden)] = 0
colnames(mes.iden) = mes.ident$taxon_name

mes.id = mes.iden %>% as.data.frame() %>% 
  add_column(.before = "Acaridae",
             Unit_quarter = str_split(rownames(mes.iden), "_", simplify = T)[,2]) %>% 
  as_tibble() %>% 
  pivot_longer(2:43,
               names_to = "Taxon",
               values_to = "Abundance") %>% 
  group_by(Unit_quarter, Taxon) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Taxon,
              values_from = Abundance)


collembola = mes.ident$taxon_name[mes.ident$class == "Collembola"]
mites = mes.ident$taxon_name[mes.ident$class == "Arachnida"] 
pauropoda = mes.ident$taxon_name[mes.ident$class == "Pauropoda"]
protura = mes.ident$taxon_name[mes.ident$class == "Protura"]

mes.tax.dens = mes.id %>% 
  # *.out is the number that were identified
  mutate(Collembola.out = rowSums(.[collembola]),
              Acari.out = rowSums(.[mites]),
          Pauropoda.out = rowSums(.[pauropoda]),
            Protura.out = rowSums(.[protura]),
         .before = "Acaridae") %>% 
  # *.in is the number we counted originally
  add_column(.before = "Collembola.out",
             Collembola.in = mes.counts$Collembola[match(.$Unit_quarter, mes.counts$Unit_quarter)]) %>% 
  add_column(.before = "Acari.out",
             Acari.in = mes.counts$Acari[match(.$Unit_quarter, mes.counts$Unit_quarter)]) %>% 
  add_column(.before = "Pauropoda.out",
             Pauropoda.in = mes.counts$Pauropoda[match(.$Unit_quarter, mes.counts$Unit_quarter)]) %>% 
  add_column(.before = "Protura.out",
             Protura.in = mes.counts$Protura[match(.$Unit_quarter, mes.counts$Unit_quarter)]) %>% 
  add_column(.after = "Protura.out",
             Symphyla = mes.counts$Symphyla[match(.$Unit_quarter, mes.counts$Unit_quarter)]) %>% 
  rowwise() %>% 
  mutate(.before = "Collembola.in",
         # we subtract pauropoda and protura that were originally mistaken for collembola
         Collembola.in = Collembola.in - (Pauropoda.out + Protura.out),
         Acari = max(Acari.in,Acari.out),
         Collembola = max(Collembola.in,Collembola.out),
         Pauropoda = Pauropoda.in + Pauropoda.out,
         Protura = Protura.in + Protura.out) %>% 
  ungroup() %>% 
  # to estimate in-core densities of the different taxa
  # we divide each taxon by the total identified, then we multiply by the number counted
  mutate(across(all_of(collembola), ~(./Collembola.out)*Collembola)) %>% 
  mutate(across(all_of(mites), ~(./Acari.out)*Acari)) %>% 
  select(-c(Collembola.in,Collembola.out,Collembola,
            Acari.in,Acari.out,Acari,
            Pauropoda.in,Pauropoda.out,
            Protura.in,Protura.out)) %>% 
  # we estimate densities per mesocosm based on densities in a core of 5cm diameter
  mutate(across(where(is.numeric), ~ (./(pi*2.5^2))*(pi*25^2)))


mes.taxonomy = mes.ident[,1:4] %>% 
  arrange(class,order,family) %>% 
  mutate(.before = class,
         group = "mesofauna") %>% 
  add_row(taxon_name = "Symphyla",
          group = "mesofauna",
          class = "Symphyla",
          order = "Symphyla") %>% 
  mutate(Abbreviation = c(# Abbreviations used in Potapov et al 2022
                          rep("Me", 7),"Ast",
                          rep("Ori", 12),
                          rep("Prost", 8),
                          rep("Cla-A", 3),
                          rep("Cla-D", 3),
                          "Cla",
                          rep("Cla-On", 2),
                          "Cla-D",rep("Cla-A", 2),
                          "Pau","Protu","Sym"))
mes.taxonomy$Abbreviation[mes.taxonomy$taxon_name %in% c("Astegistidae","Autognethidae")] = "Ori-M" # fungivores Maraun 2023
mes.taxonomy$Abbreviation[mes.taxonomy$taxon_name %in% c("Euphthiracaridae")] = "Ori-M"             # and decomp
mes.taxonomy$Abbreviation[mes.taxonomy$taxon_name %in% c("Hypochthoniidae")] = "Ori-C"
mes.taxonomy$Abbreviation[mes.taxonomy$taxon_name %in% c("Oppiidae")] = "Ori-C"                     # and fungi
mes.taxonomy$Abbreviation[mes.taxonomy$taxon_name %in% c("Oribatulidae")] = "Ori-M"                 # and decomp
mes.taxonomy$Abbreviation[mes.taxonomy$taxon_name %in% c("Protoribatidae")] = "Ori-M"               # fungivores


################################### Body size ##################################

# length measurements for protura, pauropoda, symphyla
other.lengths = read.xlsx("H:\\Jenatron\\Angelos_Jenatron_mesofauna_in_Core.xlsx",
                       sheet = "lengths") %>% 
  pivot_longer(2:4,
               names_to = "taxon_name",
               values_to = "length_mm") %>% 
  arrange(Unit_quarter, taxon_name) %>% 
  drop_na(length_mm) %>% 
  mutate(.after = "Unit_quarter",
         order = case_when(taxon_name == "Pauropoda" ~ "Tetramerocerata", 
                           taxon_name == "Protura" ~ "Protura",
                           taxon_name == "Symphyla" ~ "Symphyla")) %>% 
  mutate(.after = "taxon_name",
         .keep = "unused",
         length_micro = length_mm*1e3)

# length measurements for mites & Collembola
mes.lengths = read.xlsx("H:\\Jenatron\\Jenatron_SP6_mesofauna_data_v3.10.xlsx",
                        sheet = "qry_length_mesof")
mes.lengths$taxon_name[mes.lengths$taxon_name == "Blattisociidae_x000D_\nBlattisociidae"] = "Blattisociidae"
mes.lengths$taxon_name[mes.lengths$taxon_name == "Eosentomidae"] = "Protura"
mes.lengths$taxon_name[mes.lengths$taxon_name == "Pauropodidae"] = "Pauropoda"

mes.lengths = mes.lengths %>% 
  rename("length_micro" = "length_µm") %>% 
  drop_na(length_micro) %>% 
  
  mutate(.before = ind_ID,
         .keep = "unused",
         Unit_quarter = str_c(U,q)) %>% 
  select(Unit_quarter, class, order, family, taxon_name, length_micro) %>%
  # now we bring in protura, pauropoda, symphyla
  full_join(., other.lengths, by = join_by(Unit_quarter, order, taxon_name, length_micro))


mes.masses = mes.lengths %>% 
  # calculate individual fresh bodymasses
  mutate(FreshMass.mg = case_when(# based on Mercer 2001 10.1017/S0954102001000219
                                  # Mites
                                  order == "Trombidiformes" ~ (10^(2.124 + 2.808*log10(length_micro/1e3)))/1e3,
                                  order ==   "Mesostigmata" ~ (10^(2.064 + 2.857*log10(length_micro/1e3)))/1e3,
                                  order == "Sarcoptiformes" ~ (10^(2.146 + 2.770*log10(length_micro/1e3)))/1e3,
                                  # Collembola, Protura, Pauropoda, Symphyla
                                  TRUE ~ (10^(1.339 + 1.992*log10(length_micro/1e3)))/1e3)) %>% 
  select(-length_micro)


mes.mean.masses = mes.masses %>% 
  group_by(taxon_name) %>% 
  summarise(N = n(),
            mean.mass.mg = mean(FreshMass.mg))

# write.csv(mes.masses, "mesomasses.csv") 



# library(tidyverse)
# library(brms)
# 
# mes.masses = read.csv("mesomasses.csv", row.names = 1, header = TRUE)
# mes.mean.masses = mes.masses %>%
#   group_by(Unit_quarter, taxon_name) %>%
#   summarise(N = n(),
#             mean.mass.mg = mean(FreshMass.mg))
# 
# 
# #### partial pooling within groups, group means varying by unit ####
# mpp = brm(bf(FreshMass.mg ~ 0 + taxon_name + (0 + taxon_name|Unit_quarter)),
#           family = Gamma(link = "log"),
#           init = 0,
#           chains = 4,
#           cores = 4,
#           iter = 4000,
#           control = list(adapt_delta = 0.95),
#           backend = "cmdstanr",
#           seed = 111,
#           data = mes.masses)


