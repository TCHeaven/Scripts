```R
# Install necessary packages
install.packages("rlang", update = TRUE)
install.packages("BiocManager", update = TRUE)
library(BiocManager)
BiocManager::install("ggtree", update = TRUE)
install.packages("ape", dependencies = TRUE, update = TRUE)
install.packages("phangorn", dependencies = TRUE, update = TRUE)

# Load required libraries
library(ggplot2)
library(ggtree)
library(ape)
library(phangorn)
library(dplyr)
library(heatmaply)

# Set working directory
setwd("C:/Users/did23faz/OneDrive - Norwich Bioscience Institutes/Desktop/R")

# Read the distance matrix from CSV file
dist_matrix <- read.csv("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193+ligustri_mperc.csv", row.names = 1)
dist_matrix <- read.csv("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193+ligustri_genic.csv", row.names = 1)

# Convert the matrix to a distance object
diss_matrix <- as.dist(dist_matrix)

# Build the phylogenetic tree using Neighbor-Joining
phy_tree <- nj(diss_matrix)
plot(phy_tree)

# Root the tree with "ligustri"
phy_tree_rooted <- root(phy_tree, outgroup = "ligustri")
plot(phy_tree_rooted)

# Set plot parameters with smaller font size
p <- ggtree(phy_tree_rooted) + 
  geom_tiplab(size = 1.5) +
  geom_tippoint(size = 0) +
  theme(legend.text = element_text(size = 5), text = element_text(size = 5))

# Shorten the "ligustri" branch
ggtree(phy_tree_rooted) + geom_text(aes(label=node), hjust=-.3)
p$data[p$data$node %in% c(194), "x"] = mean(p$data$x)

# Plot the tree
print(p)

#Get sample info
sample_info <- read.csv("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/PCA_file_host.csv", row.names = 1)
row.names(sample_info) <- sample_info[, 1]
sample_info <- sample_info[, -1]
tip_names <- phy_tree$tip.label
row_names <- rownames(sample_info)
matching_rows <- row_names %in% tip_names
sample_info <- sample_info[matching_rows, ]
sample_info <- sample_info[, -c(1, 2, 3, 4, 6, 11)]
sample_info <- sample_info[, c(5, 3, 2, 4, 1)]

#Make a dataframe for info about who sampled each sample
Sampler <- data.frame(sample_info[, 1])
row.names(Sampler) <- row.names(sample_info)
#Plot sample collector info with the phylogeny
q <- gheatmap(p, Sampler, offset=0, width=0.5, font.size=3, 
              colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("Bass", "JIC"), 
                    values=c("steelblue", "firebrick"), name="Sampler")
# Adjust plot margins and size
q <- q + theme(plot.background = element_rect(fill = "white"),
               panel.grid.major = element_blank(),
               plot.margin = margin(0, 0, 30, 0, "pt"),
               panel.grid.minor = element_blank(),
               plot.title = element_text(hjust = 0.5),
               axis.text.x = element_text(angle = -45, hjust = 0))

# Increase the size of the plotting area
q <- q + coord_cartesian(clip = "off")
print(q)

#Make a dataframe for info about geography of each sample
Geo <- data.frame(sample_info[, c(2,3)])
row.names(Geo) <- row.names(sample_info)
#Plot geography info with the phylogeny
r <- gheatmap(p, Geo, offset=0, width=0.5, font.size=3, 
              colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("Europe", "Australia","Asia","America","Africa","Italy","United Kingdom", "Belgium","France", "Spain", "US", "Israel", "Kenya", "Netherlands", "Germany", "Greece", "Zimbabwe", "Algeria", "Belarus", "Armenia", "Japan", "Argentina", "Switzerland", "South Korea","Hungary", "China", "Chile", "Tunisia","NA"), 
                    values=c("steelblue", "gold", "forestgreen", "firebrick", "grey86","blue", "dodgerblue", "cornflowerblue", "deepskyblue", "Blue4", "darkred", "darkslategray1", "gray56", "navyblue", "royalblue", "cadetblue1", "grey21", "ivory2", "turquoise1", "turquoise4", "green", "orangered", "slateblue", "springgreen","mediumpurple1", "chartreuse4", "red4", "grey0","snow" ), name="Geography")
# Adjust plot margins and size
r <- r + theme(plot.background = element_rect(fill = "white"),
               panel.grid.major = element_blank(),
               plot.margin = margin(0, 0, 30, 0, "pt"),
               panel.grid.minor = element_blank(),
               plot.title = element_text(hjust = 0.5),
               axis.text.x = element_text(angle = -45, hjust = 0))

# Increase the size of the plotting area
r <- r + coord_cartesian(clip = "off")
print(r)

#Make a dataframe for info about host of each sample
Host <- data.frame(sample_info[, c(4,5)])
row.names(Host) <- row.names(sample_info)
#Plot host info with the phylogeny
s <- gheatmap(p, Host, offset=0, width=0.5, font.size=3, 
              colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("Tobacco", "Solanacae", "Brassica", "Lab", "Others", "Sugar beet", "Prunus", "Potato", "Pepper", "Tomato", "Canola", "Chinese cabbage", "Kales", "Cabbage", "Brussel sprouts", "Broccoli", "Arabidopsis", "Peas", "Chrysantemum", "unknown", "Aubergine", "Weed", "Almond", "Peach", "Nectarine" ), 
                    values=c("gold", "firebrick", "forestgreen", "grey0", "steelblue", "grey86", "hotpink", "firebrick4", "red", "tomato", "springgreen", "springgreen4", "turquoise", "olivedrab", "darkseagreen1", "chartreuse", "forestgreen", "slateblue", "blue", "steelblue", "royalblue", "darkblue", "mediumorchid1", "darkorchid4", "magenta"), name="Host")
# Adjust plot margins and size
s <- s + theme(plot.background = element_rect(fill = "white"),
               panel.grid.major = element_blank(),
               plot.margin = margin(0, 0, 30, 0, "pt"),
               panel.grid.minor = element_blank(),
               plot.title = element_text(hjust = 0.5),
               axis.text.x = element_text(angle = -45, hjust = 0))

# Increase the size of the plotting area
s <- s + coord_cartesian(clip = "off")
print(s)

#All
t <- gheatmap(p, sample_info, offset=0, width=0.5, font.size=3, 
              colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("Bass", "JIC", "Europe", "Australia","Asia","America","Africa","Italy","United Kingdom", "Belgium","France", "Spain", "US", "Israel", "Kenya", "Netherlands", "Germany", "Greece", "Zimbabwe", "Algeria", "Belarus", "Armenia", "Japan", "Argentina", "Switzerland", "South Korea","Hungary", "China", "Chile", "Tunisia", "Tobacco", "Solanacae", "Brassica", "Lab", "Others", "Sugar beet", "Prunus", "Potato", "Pepper", "Tomato", "Canola", "Chinese cabbage", "Kales", "Cabbage", "Brussel sprouts", "Broccoli", "Arabidopsis", "Peas", "Chrysantemum", "unknown", "Aubergine", "Weed", "Almond", "Peach", "Nectarine","NA" ), 
                    values=c("steelblue", "firebrick", "steelblue", "gold", "forestgreen", "firebrick", "grey86","blue", "dodgerblue", "cornflowerblue", "deepskyblue", "Blue4", "darkred", "darkslategray1", "gray56", "navyblue", "royalblue", "cadetblue1", "grey21", "ivory2", "turquoise1", "turquoise4", "green", "orangered", "slateblue", "springgreen","mediumpurple1", "chartreuse4", "red4", "grey0", "gold", "firebrick", "forestgreen", "grey0", "steelblue", "grey86", "hotpink", "firebrick4", "red", "tomato", "springgreen", "springgreen4", "turquoise", "olivedrab", "darkseagreen1", "chartreuse", "forestgreen", "slateblue", "blue", "steelblue", "royalblue", "darkblue", "mediumorchid1", "darkorchid4", "magenta","snow"), name="All")

# Adjust plot margins and size
t <- t + theme(plot.background = element_rect(fill = "white"),
               panel.grid.major = element_blank(),
               plot.margin = margin(0, 0, 30, 0, "pt"),
               panel.grid.minor = element_blank(),
               plot.title = element_text(hjust = 0.5),
               axis.text.x = element_text(angle = -45, hjust = 0))

# Increase the size of the plotting area
t <- t + coord_cartesian(clip = "off")

# Print the plot
print(t)

```
mv * /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/.
echo DONE
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/tmp_54396579
exit
exit

