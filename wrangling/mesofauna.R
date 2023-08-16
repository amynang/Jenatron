library(openxlsx)
library(tidyverse)

# counted mesofauna (in core abundances)
mes.counts = read.xlsx("H:\\Jenatron\\Angelos_Jenatron_mesofauna_in_Core.xlsx",
                       sheet = "counts")
# identified Collembola & Acari
mes.ident = read.xlsx("H:\\Jenatron\\Jenatron_SP6_mesofauna_data_v3.10.xlsx",
                sheet = "qry_abund_mesof")
# correct name
mes.ident$taxon_name[mes.ident$taxon_name == "Blattisociidae_x000D_\nBlattisociidae"] = "Blattisociidae"

# length measurements for mites & Collembola
mes.lengths = read.xlsx("H:\\Jenatron\\Jenatron_SP6_mesofauna_data_v3.10.xlsx",
                sheet = "qry_length_mesof")
mes.lengths$taxon_name[mes.lengths$taxon_name == "Blattisociidae_x000D_\nBlattisociidae"] = "Blattisociidae"
mes.lengths = mes.lengths %>% 
  rename("length_micro" = "length_µm") %>% 
  drop_na(length_micro) %>% 
  mutate(.before = ind_ID,
         .keep = "unused",
         Unit_quarter = str_c(U,q)) %>% 
  select(Unit_quarter, class, order, family, taxon_name, length_micro)

d.4 = mes.ident[,5:244] %>% t()
d.4[is.na(d.4)] = 0
colnames(d.4) = mes.ident$taxon_name

d.5 = d.4 %>% as.data.frame() %>% 
  add_column(.before = "Acaridae",
             Unit_quarter = str_split(rownames(d.4), "_", simplify = T)[,2]) %>% 
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

d.6 = d.5 %>% 
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


rownames(d.4) = str_split(rownames(d.4), "_", simplify = T)[,2]
