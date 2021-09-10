########  CasTn analysis
## Auth: Lidimarie Trujillo Rodriguez
## email: Ltrujillo@ufl.edu
## Reisch Lab
## Microbiology and Cell Science
## University of Florida


###  Import dataset
###############################
data <- read.csv("raw.csv", header = TRUE, skip = 6)
library(ggplot2)
library(dplyr)
library(formattable)
library(table1)

data <- data[data$Grow_LB > 0,] %>%
  rowwise() %>%
  mutate(
    Phenotype = if_else(is.na(Fluoresence), sum(c(
      Grow_LB, Grow_M9, Grow_M9_AA
    )), Fluoresence),
    Genotype = if_else(is.na(RE.LE_PCR -
                               WT_PCR), RE.LE_PCR, RE.LE_PCR - WT_PCR)
  ) # Phenotype/Genotype score, only keep mutants that grew in LB

#Statistics by biological replicate
#############################
table <- data %>%
  group_by(Bacteria, Target.oligo, CasTn, Target.Gene, Rep) %>%
  count()

table <-
  cbind(table, (
    data %>% group_by(Bacteria, Target.oligo, CasTn, Target.Gene, Rep) %>% tally(Phenotype == 0 |
                                                                                   Phenotype == 2, name = "On-target")
  )$`On-target`)

colnames(table)[length(colnames(table))] <- "On-target"
table$"On-target%" <- table$`On-target` / table$n
#### Phenotypic On-target
## 0 is on-target fluoresence phenotype, 2 is on-target auxotrophy phenotype

table = cbind(table, (
  data %>% group_by(Bacteria, Target.oligo, CasTn, Target.Gene, Rep) %>% tally(!is.na(Genotype) ==
                                                                                 TRUE, name = "On-target")
)$`On-target`)
colnames(table)[length(colnames(table))] <- "PCR n"
table = cbind(table, (
  data %>% group_by(Bacteria, Target.oligo, CasTn, Target.Gene, Rep) %>% tally(Genotype >
                                                                                 0, name = "On-target")
)$`On-target`)
colnames(table)[length(colnames(table))] <- "PCR On-target"
table$"PCR On-target%" = table$`PCR On-target` / table$`PCR n`
#### Genotype On-target
## 1 is on-target , 0 is off-target


table = cbind(table, (
  data %>% group_by(Bacteria, Target.oligo, CasTn, Target.Gene, Rep) %>% tally(!is.na(pUC_PCR) ==
                                                                                 TRUE, name = "On-target")
)$`On-target`)
colnames(table)[length(colnames(table))] <- "pUC n"
table = cbind(table, (
  data %>% group_by(Bacteria, Target.oligo, CasTn, Target.Gene, Rep) %>% tally(pUC_PCR >
                                                                                 0, name = "On-target")
)$`On-target`)
colnames(table)[length(colnames(table))] <- "pUC PCR"
table$"plasmid insert" = if_else(is.nan(table$`pUC PCR` / table$`pUC n`),
                                 0,
                                 table$`pUC PCR` / table$`pUC n`)
#### cointegrate frequency
## 1 is cointegrate , 0 is no cointegrate

#### Summary statistics by biological replicate
table_stat = table %>% group_by(Bacteria, Target.oligo, CasTn, Target.Gene) %>% summarise(
  mean.pheno = mean(`On-target%`) * 100,
  #Phenotype on-target %
  sd.pheno = sd(`On-target%`) * 100,
  mean.pcr = mean(`PCR On-target%`) * 100,
  #Genotype on-target %
  sd.pcr = sd(`PCR On-target%`) * 100,
  mean.plasmid = mean(`plasmid insert`) * 100,
  #Cointegrate %
  sd.plasmid = sd(`plasmid insert`) * 100
)

table_stat$mean.pheno[grep("Vio", table_stat$Target.oligo)] <- c(100)
# Violacein all mutants presented violet
table = table[!table$Target.oligo == "before-2517",]
# Remove 2517 before arginine addition

# Table with percentages
formattable(
  table_stat,
  list(
    `Bacteria` = formatter("span", style = ~ style(color = "Black", font.style = "italic")),
    `mean.pheno` = color_tile("transparent", "springgreen4"),
    `mean.pcr` = color_tile("transparent", "chartreuse4"),
    `mean.plasmid` = color_tile("transparent", "coral")
  )#List
)#formattable
########################################################

#Plotting and Statistics
#Overall on-target efficiency and sgRNA length
############################
library(ggpubr)
my_comparisons <-
  list(
    c("Agrobacterium fabrum", "Pseudomonas putida"),
    c("Pseudomonas putida", "Burkholderia thailandensis"),
    c("Agrobacterium fabrum", "Burkholderia thailandensis")
  )

##Genotype
ggplot(table_stat, aes(x = Bacteria, y = mean.pcr)) +
  geom_boxplot() +
  geom_jitter(data = table_stat) +
  theme_classic() +
  xlab("Bacteria") +
  ylab("On-target genotype (%)") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  stat_compare_means(method = "kruskal.test", label.y = 150) +
  labs(title = "On-target genotype")
ggsave(
  "on-target genotype.tiff",
  width = 6,
  height = 4,
  units = "in",
  device = "tiff",
  dpi = 320
)

##Phenotype
ggplot(table_stat, aes(x = Bacteria, y = mean.pheno)) +
  geom_boxplot() +
  geom_jitter(data = table_stat) +
  theme_classic() +
  xlab("Bacteria") +
  ylab("On-target phenotype (%)") +
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
  #stat_compare_means(method = "kruskal.test", label.y = 150)+
  labs(title = "On-target phenotype")
ggsave(
  "on-target phenotype.tiff",
  width = 6,
  height = 4,
  units = "in",
  device = "tiff",
  dpi = 320
)

############ Target length
length = read.csv("target.length.csv", header = TRUE)
table_stat = merge(
  table_stat,
  length,
  by.x = 'Target.oligo',
  by.y = 'Target',
  all.x = TRUE,
  all.y = FALSE
)
table_stat$Comparison[1:4] <- "2300"
table_stat$Comparison[12] <- "2517"
table_stat$Comparison[10] <- "2511"
table_stat$Comparison[15:16] <- "2523"
table_stat$Comparison[19] <- "2542"
table_stat$Comparison[40:42] <- "2708"
table_stat3 = table_stat[!is.na(table_stat$Comparison),]

## Statistics different lenth sgRNA
shapiro.test(table_stat3$mean.pcr) ## Non-parametric
kruskal.test(table_stat3$mean.pcr ~ table_stat3$length.bp.) ##global Geno
kruskal.test(table_stat3$mean.pheno ~ table_stat3$length.bp.) ##global Pheno

### Statistics sgRNA length vs genotype by bacteria
bt_sglength = table_stat3[table_stat3$Bacteria == "Burkholderia thailandensis" , ]
kruskal.test(bt_sglength$mean.pcr ~ bt_sglength$length.bp.)

af_sglength = table_stat3[table_stat3$Bacteria == "Agrobacterium fabrum" , ]
kruskal.test(af_sglength$mean.pcr ~ af_sglength$length.bp.)

pp_sglength = table_stat3[table_stat3$Bacteria == "Pseudomonas putida" , ]
kruskal.test(pp_sglength$mean.pcr ~ pp_sglength$length.bp.)

########################################################


## Pearson correlation
########################
cor_all = cor.test(table_stat$mean.pheno,
                   table_stat$mean.pcr,
                   use = "complete.obs",
                   method = "pearson")

cor_all

table_stat2 <-
  table_stat[table_stat$Bacteria == "Burkholderia thailandensis", ]
cor_bt = cor.test(table_stat2$mean.pheno,
                  table_stat2$mean.pcr,
                  use = "complete.obs",
                  method = "pearson")
cor_bt

table_stat2 <-
  table_stat[table_stat$Bacteria == "Agrobacterium fabrum", ]
cor_af = cor.test(table_stat2$mean.pheno,
                  table_stat2$mean.pcr,
                  use = "complete.obs",
                  method = "pearson")
cor_af

table_stat2 <-
  table_stat[table_stat$Bacteria == "Pseudomonas putida", ]
cor_pp = cor.test(table_stat2$mean.pheno,
                  table_stat2$mean.pcr,
                  use = "complete.obs",
                  method = "pearson")
cor_pp
########################################################

## Cointegrate statistics
######################
#Overall
plasmid =  data %>%   group_by() %>%
  summarise(
    N = n(),
    mean = mean(pUC_PCR, na.rm = TRUE),
    median = median(pUC_PCR, na.rm = TRUE),
    sd = sd(pUC_PCR, na.rm = TRUE),
    std = (sd(pUC_PCR, na.rm = TRUE) / sqrt(length(
      sd(pUC_PCR, na.rm = TRUE)
    )))
  )
plasmid

#By bacteria and CasTn
plasmid =  data %>%   group_by(Bacteria, CasTn) %>%
  summarise(
    N = n(),
    mean = mean(pUC_PCR, na.rm = TRUE),
    median = median(pUC_PCR, na.rm = TRUE),
    sd = sd(pUC_PCR, na.rm = TRUE),
    std = (sd(pUC_PCR, na.rm = TRUE) / sqrt(length(
      sd(pUC_PCR, na.rm = TRUE)
    )))
  )
plasmid

# By bacteria, biological replicate, target and CasTn
plasmid =  data %>%   group_by(Bacteria, Target.oligo, CasTn, Target.Gene, Rep) %>%
  summarise(
    N = n(),
    mean = mean(pUC_PCR, na.rm = TRUE),
    median = median(pUC_PCR, na.rm = TRUE),
    sd = sd(pUC_PCR, na.rm = TRUE),
    std = (sd(pUC_PCR, na.rm = TRUE) / sqrt(length(
      sd(pUC_PCR, na.rm = TRUE)
    )))
  )
plasmid
## sd = standard deviation, std = standard error

### Table of cointegrate formation frequenncy by bacteria
plasmid = plasmid[-which(plasmid$Bacteria == "Agrobacterium fabrum" &
                           plasmid$CasTn == "pBLLct2"), ]
plasmid$mean = as.numeric(plasmid$mean)
table1(~ mean | Bacteria, data = plasmid)

## Graph and comparion cointegrate
ggplot(plasmid, aes(x = Bacteria, y = mean)) +
  geom_boxplot(geom = 'errorbar',
               linetype = 1,
               width = 0.5) +
  geom_jitter(data = plasmid, width = 0.1) +
  ylab("Plasmid insertion frequency") +
  xlab("") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  stat_compare_means(method = "kruskal.test", label.y = 1.5) +
  theme_bw() +
  ggtitle("Whole plasmid insertion")
########################################################

# Transformation efficiency
#######################################################
te = read.csv("CFU.csv", header = FALSE, skip = 1)
colnames(te) = c("Bacteria", "CasTn", "Target", "CFU")
te$CFU = as.numeric(as.character(te$CFU))
te$bac_cas = paste(te$Bacteria, te$CasTn, sep = "_")
te$bac_target = paste(te$Bacteria, te$Target, sep = "_")

### CFU by bacteri and CasTn
table1(
  ~ CFU |
    bac_cas ,
  data = te,
  transpose = FALSE,
  overall = FALSE,
  footnote = "Transformation efficiency by CasTn.",
  rowlabelhead = "CasTn"
)

### CFU by bacteria and target
te1 = te[te$Bacteria == unique(te$Bacteria)[1], ] #BT
table1(
  ~ CFU |
    bac_target ,
  data = te1,
  transpose = TRUE,
  overall = FALSE,
  footnote = "BT Transformation efficiency by sgRNA.",
  rowlabelhead = "sgRNA Target"
)

te1 = te[te$Bacteria == unique(te$Bacteria)[2], ] #PP
table1(
  ~ CFU |
    bac_target ,
  data = te1,
  transpose = TRUE,
  overall = FALSE,
  footnote = "PP Transformation efficiency by sgRNA.",
  rowlabelhead = "sgRNA Target"
)

te1 = te[te$Bacteria == unique(te$Bacteria)[3], ] #AF
table1(
  ~ CFU |
    bac_target ,
  data = te1,
  transpose = TRUE,
  overall = FALSE,
  footnote = "AF Transformation efficiency by sgRNA.",
  rowlabelhead = "sgRNA Target"
)
########################################################

# TnSeq data graphing
######################
library(scales)
bedgraph = list.files(pattern = "*\\.bedgraph",
                      full.names = TRUE,
                      recursive = TRUE)


bedgraph <-
  bedgraph[grep("bin1000", bedgraph)]##### Off-target analysis



#Import .bedgraph files
files = list()
for (i in seq_along(bedgraph)) {
  files[i] = lapply(bedgraph[i], function(x)
    read.delim(x, header = FALSE)) ##Import data
  names(files)[i] <-
    gsub(
      ".*/|[.bedgraph]",
      "",
      list.files(
        pattern = "*\\.bedgraph",
        full.names = FALSE,
        recursive = TRUE
      )
    )[i]
  files[i] = lapply(files[i], setNames, c("Chrom", "TnStart", "TnStop", "depth")) ## Setting DF names and Colnames
}

#names(files) <- gsub(".*/|[.bedgraph]","",bedgraph) ##### Off-target analysis

##Make a new column of percentages
files = lapply(files, function(x) {
  cbind(x, log10 = with(x, ifelse(depth == 0, 0, c(log10(
    depth
  )))))
}) ### Transform log10
files = lapply(files, function(x) {
  cbind(x, percentage = with(x, ifelse(depth < 1, 0, c(depth / sum(
    depth
  )))))
})
name <- as.list(names(files)) ##Make list of names
for (i in seq_along(files[grep("P.*nomic", names(files))])) {
  files[grep("P.*nomic", names(files))][[i]][, 1] <- "chr1"
} #Pputida first columns to chr1

## Separate each DF in list by chromosome/plasmid
files2 <-
  lapply(files, function(x)
    split(x, x[, 1])) ##Split into unique chromosomes
files2 <-
  sapply(files2, `[`, 1:6, simplify = TRUE) ##Iterate over every 6 elements in list
files2 <- files2[lapply(files2, length) > 0] ##Remove empty DF

## Give unique names to df
name <-
  sapply(seq_along(files), function(i)
    paste(names(files[i]), unique(files[[i]][, 1]), sep = " "), simplify = TRUE)
name <- unlist(name)
names(files2) <- name
names(files2)


#### Off-targets across genome ##################
off_targets = lapply(files, function(x) {
  x[x$percentage > 0.001, ]
})
a = files3$At2511o2nomic
############################################


# Plotting TnSeq Functions
circular.graph <- function(dataset, x, y, name, angle) {
  require(ggplot2)
  require(scales)
  
  mindepth <- dataset[dataset[, y] == min(dataset[, y]), x]
  medepth <- dataset[dataset[, y] == median(dataset[, y]), x]
  maxdepth <-
    dataset[dataset[, y] == max(dataset[, y]), x] ## Calculate higher limit
  sum_y <- sum(dataset[, y])
  dataset[, ncol(dataset) + 1] <- dataset[, y] / sum_y
  
  colnames(dataset)[x] <- "coord"
  colnames(dataset)[y] <- "depth"
  #colnames(dataset)[ncol(dataset)] <- "percentage"
  
  ggplot(data = dataset, aes(x = dataset[, x], y = dataset[, y])) + ## Import data
    geom_bar(stat = "identity", width = 2, aes(color = c(dataset[, x] == maxdepth))) + ##Bargraph depth, color highest
    scale_color_manual(values = c('black', 'red')) +
    geom_text(
      data = dataset[dataset$percentage > 0.2,],
      aes(
        label = scales::percent(percentage),
        x = coord,
        ##Label all above x percent
        y = depth
      ),
      stat = "identity",
      size = 3,
      vjust = -0.5,
      hjust = -0.3
    ) +
    ## adds percentage
    geom_segment(
      data = dataset,
      aes(
        x = min(dataset[, x]),
        y = 50000,
        xend = max(dataset[, x]),
        yend = 50000
      ),
      colour = "grey36",
      linetype = "dashed",
      alpha = 1,
      size = 0.3 ,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = dataset,
      aes(
        x = min(dataset[, x]),
        y = 150000,
        xend = max(dataset[, x]),
        yend = 150000
      ),
      colour = "grey36",
      linetype = "dashed",
      alpha = 1,
      size = 0.3 ,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = dataset,
      aes(
        x = min(dataset[, x]),
        y = 250000,
        xend = max(dataset[, x]),
        yend = 250000
      ),
      colour = "grey36",
      linetype = "dashed",
      alpha = 1,
      size = 0.3 ,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = dataset,
      aes(
        x = min(dataset[, x]),
        y = 500000,
        xend = max(dataset[, x]),
        yend = 500000
      ),
      colour = "grey36",
      linetype = "dashed",
      alpha = 1,
      size = 0.3 ,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = dataset,
      aes(
        x = min(dataset[, x]),
        y = 375000,
        xend = max(dataset[, x]),
        yend = 375000
      ),
      colour = "grey36",
      linetype = "dashed",
      alpha = 1,
      size = 0.3 ,
      inherit.aes = FALSE
    ) +
    ## add y scale
    annotate(
      "text",
      x = rep(max(dataset$coord) - 150000, 5),
      y = c(75000, 175000, 275000, 400000, 525000),
      #75000
      label = c("50K RPM", "150K RPM", "250K RPM", "375K RPM", "500K RPM") ,
      color = "black",
      size = 3 ,
      angle = angle,
      hjust = 1
    ) +
    ylim(-500000, c(max(dataset[, y]) + 100000)) +
    theme_bw() +
    scale_x_continuous(
      name = c(name),
      breaks = c(0e+00, 1e+06, 2e+06, 3e+06, 4e+06, 5e+06),
      labels = c("0Mbp", "1Mbp", "2Mbp", "3Mbp", "4Mbp", "5Mbp")
    ) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    coord_polar()
} #Circulat plot

linear.graph.plasmids <- function(dataset, x, y, name, z) {
  require(ggplot2)
  require(scales)
  
  colnames(dataset)[x] <- "coord"
  colnames(dataset)[y] <- "depth"
  
  dataset = dataset[dataset[, y] > 0, ]
  
  ggplot(data = dataset,
         aes(x = dataset[, x], y = dataset[, y]),
         color = 'balck') + ## Import data
    geom_bar(stat = "identity", width = 2, aes(color = 'black')) +
    scale_color_manual(values = c('black')) +
    theme_bw() +
    ylab("Reads per Million") +
    scale_x_continuous(name = c("coordinate")) +
    scale_y_log10() +
    ggtitle(c(name)) +
    theme(legend.position = "none")
}##Basic plot

linear.graph.simple <- function(dataset, x, y, name, z) {
  require(ggplot2)
  require(scales)
  
  colnames(dataset)[x] <- "coord"
  colnames(dataset)[y] <- "depth"
  
  dataset = dataset[dataset[, y] > 0, ]
  
  ggplot(data = dataset,
         aes(x = dataset[, x], y = dataset[, y]),
         color = 'balck') + ## Import data
    geom_bar(stat = "identity", width = 2, aes(color = 'black')) +
    scale_color_manual(values = c('black')) +
    theme_bw() +
    ylab("Reads per Million") +
    scale_y_continuous(
      labels = function(x)
        format(x, scientific = TRUE),
      limits = c(0, z)
    ) +
    scale_x_continuous(
      name = c("Genomic coordinate"),
      breaks = c(0e+00, 1e+06, 2e+06, 3e+06, 4e+06, 5e+06, 6e+06),
      labels = c("0Mbp", "1Mbp", "2Mbp", "3Mbp", "4Mbp", "5Mbp", "6Mbp")
    ) +
    ggtitle(c(name)) +
    theme(legend.position = "none")
}##Basic plot

linear.graph <- function(dataset, x, y, name, coord1, coord2, z) {
  require(ggplot2)
  require(scales)
  dataset <- dataset[dataset[, y] > 0,]
  mindepth <- dataset[dataset[, y] == min(dataset[, y]), x]
  medepth <- dataset[dataset[, y] == median(dataset[, y]), x]
  maxdepth <-
    dataset[dataset[, y] == max(dataset[, y]), x] ## Calculate higher limit
  sum_y <- sum(dataset[, y])
  dataset[, ncol(dataset) + 1] <- dataset[, y] / sum_y
  
  
  colnames(dataset)[x] <- "coord"
  colnames(dataset)[y] <- "depth"
  colnames(dataset)[6] <- "proportion"
  
  
  ggplot(data = dataset, aes(x = dataset[, x], y = dataset[, y])) + ## Import data
    geom_bar(
      stat = "identity",
      width = 2,
      aes(color =  dataset[, x] >= coord1 &
            dataset[, x] <= coord2)
    ) + ##Bargraph depth, color highest
    scale_color_manual(values = c('black', 'red')) +
    geom_text(
      data = dataset[dataset[, x] >= coord1 &
                       dataset[, x] <= coord2, ],
      aes(
        label = scales::percent(sum(proportion)),
        x = coord2,
        ##Label all above x percent
        y = max(depth)
      ),
      stat = "identity",
      size = 3 ,
      vjust = 0.7,
      hjust = 1.3
    ) +
    
    ## adds percentage
    theme_bw() +
    ylab("Reads per Million") +
    scale_y_continuous(
      labels = function(x)
        format(x, scientific = TRUE),
      limits = c(0, z)
    ) +
    scale_x_continuous(
      name = c("Genomic coordinate"),
      breaks = c(0e+00, 1e+06, 2e+06, 3e+06, 4e+06, 5e+06, 6e+06),
      labels = c("0Mbp", "1Mbp", "2Mbp", "3Mbp", "4Mbp", "5Mbp", "6Mbp")
    ) +
    ggtitle(c(name)) +
    theme(legend.position = "none")
} #Colro/Percent plot
##Use original proportion and define color coordinates
#linear.graph(bed,2,4,"Name") ##Test
## Function to graph circular chromosome TnSeq hits
## data = dataset
## x= col number with Chr coordinate
## y = col number for depth (TnSeq read count)
## z = upper limit y coordinate
## name = name to parse on to ggtitle


plots.circ <- mapply(function(x, y)
  circular.graph(x, 2, 4, y, angle = 25)
  ,
  files2,
  name,
  SIMPLIFY = FALSE) #Circular ploting

plots.lin <-
  mapply(function(x, y)
    linear.graph.plasmids(x, 2, 4, y, 1e+06),
    files2,
    name,
    SIMPLIFY = FALSE) #Basic plotting for plasmids

plots.lin <- mapply(function(x, y)
  linear.graph.simple(x, 2, 4, y, 5e+04),
  files2,
  name,
  SIMPLIFY = FALSE)


plots.lin <-
  mapply(function(x, y)
    linear.graph(x, 2, 4, y, 323000, 325000, 1e+06),
    files2,
    name,
    SIMPLIFY = FALSE) #Linear plotting with color and percentage

## Select plot for viewing (examples)
names(files) # 1 Af CysD-1 target (At2508)
plots.circ[1]
plots.lin[1]

## Target coordinates 1Kbp bin for plotting color and percentage.
# argC Bt - 3554000,3555000
# trpF Bt - 797000, 798000 #BtS2.* chr2
# eyfp Bt - 323000,325000 #Bt2686.*|Bt2708.*
# eyfp Pp - 6171000, 6172000 #P2686.*|P2708.*
# Pp 2544 - 2263000, 2264000 #P2544.*
# Pp 2542 - 5881000, 5882000 #P2542
# AF eyfp - 1771000, 1772000 #At2686.*
# AF CysD - 816000,817000 #At25.*|At27.*

# linear plots bulks ggave
mapply(
  function(x, y)
    ggsave(
      paste0(y, ".png"),
      x,
      device = "png",
      units = "in",
      height = 1.5,
      width = 6,
      dpi = 620
    ),
  plots.lin,
  name
)
##save

#Saving individual plots
#plots2 <-plots.lin[c(2,6)]
#name2 <- name[c(2,6)]
#plots2 <- plots.lin[grep(".*omic.*",names(plots.lin))] ##Select plots to save
#name2 <- name[grep(".*omic.*",name)]; name2 ##Name of select plots to save
#mapply(function(x,y)ggsave(paste0(y,".png"),x,path = ., device = "png",units = "in",height = 1.5,width = 6, dpi =620),plots2,name2)

## Save circular plot
#plots2 <-plots.circ[c(5)]
#name2 <- name[c(5)]
#mapply(function(x,y)ggsave(paste0(y,".png"),x,path = ., device = "png",units = "in",height = 6,width = 6, dpi =620),plots2,name2)
########################################################

# Plotting zoomed in graphs
##############################################
data <- read.delim("Pp2544-bin1.bedgraph", header = TRUE)
colnames(data) <- c("Chrom", "TnStart", "TnStop", "depth")
## bedgraphs with bin= 1 bp for zoomed-in graphs

data <-
  data[data$TnStart > 2263000 & data$TnStop < 2264000, ] ##Pp trpF

## Ranges used to plot
## 322300 323000  Bt eyfp
## 3554000 3555000 Bt argC
## 5881000 5882000 Pp serA
## 2263000 data$TnStop Pp trpF (example above)
## 6171000 6172000 Pp eyfp
## 816000 816400 Af cysD
## 2795600 2796600 Af leuB
## 1771250 1771500 Af eyfp

df2 <- expand.grid(
  lineend = c('round', 'butt', 'square'),
  linejoin = c('round', 'mitre', 'bevel'),
  stringsAsFactors = FALSE
)
inst_range <-
  max(data[data$depth > 1e5, ]$TnStart) - min(data[data$depth > 1e5, ]$TnStart)
inst_dist <- max(data[data$depth > 1e5, ]$TnStart) - 2263200

zoom_plot <- ggplot(data = data, aes(x = TnStart, y = depth)) +
  geom_bar(stat = "identity") +
  xlab("Chromosome coordinate") +
  ylab("Read counts RPM") +
  theme_minimal() +
  scale_y_continuous(
    labels = function(x)
      format(x, scientific = TRUE)
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )
zoom_plot
########################################################