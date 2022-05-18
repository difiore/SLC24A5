seq <- c("CCCTTGGATTGTCTCAGGATGTTGCAGGCACAACTTTCATGGCAGCGGGCAGTTCAGCTCCTGAATTAGTTACTGCTTTCCTAG", "CCCTTGGATTGTCTCAGGATGTTGCAGGCGCAACTTTCATGGCAGCGGGCAGTTCAGCTCCTGAATTAGTTACTGCTTTCCTAG")
allele_names <- c("Ref", "Alt")

allele_table <- tibble(allele = allele_names, seq=seq)

library(tidyverse)
f <- read_tsv("~/Development/Repos/SLC24A5/rs1426654_frequency.tsv",skip = 12)
f <- f %>% filter(`#Study` == "1000Genomes") %>%
  rowwise() %>%
  mutate(Ref = as.numeric(str_split(`Ref Allele`, "=")[[1]][2]),
         Alt = as.numeric(str_split(`Alt Allele`, "=")[[1]][2])) %>%
  filter(Group == "Sub") %>%
  select(-c("Ref Allele", "Alt Allele", "BioProject ID","BioSample ID"))

t <- tibble(name=character(), pop=character(), allele = character())
pop_size <- 1000
for (p in f$Population) {
  pop <- f %>% filter(Population == p)
  for (i in 1:pop_size){
    allele <- if_else(runif(1, min = 0, max = 1) < pop$Ref, "Ref","Alt")
    r <- tibble(name = paste0(pop$Population, i), pop = pop$Population, allele = allele)
    t <- bind_rows(t,r)
  }
}
t <- left_join(t, allele_table, by = c("allele" = "allele"))
head(t)

# shuffle them up
t <- t %>%
  group_by(pop) %>%
  sample_n(100) %>%
  arrange(name) %>%
  select(name, seq)

write_csv(t, "~/Desktop/SLC24A5_exon3.csv")



# create fasta
pops <- c("Europe", "African")

t2 <- tibble(name=character(), pop=character(), allele = character())
for(i in pops){
  r <- t %>% filter(pop == i) %>% slice_sample(n=8)
  t2 <- bind_rows(t2, r)
}

for (i in 1:nrow(t2)){
  header <- paste0(">", t2[i,]$name)
  write_lines(header, "~/Desktop/temp.txt", append = TRUE)
  write_lines(t2[i,]$seq, "~/Desktop/temp.txt", append = TRUE)
}
