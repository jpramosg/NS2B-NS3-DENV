library(foreach)
library(doParallel)
library(bio3d)
library(ggplot2)
#demo("pdb")
#demo("pca")
#demo("md")

#### Calculate PCA for molecular dynamics simulations #### 

#### Author: Jean Pierre Ramos

#pc$au[,1] # relative PC1
#pc$au[,2] # relative PC2
#pc$au[,3] # relative PC3

path <- "/home/jpramosg/Desktop/MD/dengue/data/APO/"
setwd(path)

pdb_apo <- read.pdb("protein_only.pdb") 
dcd_apo <- read.dcd("protein_only_r1.dcd")

path <- "/home/jpramosg/Desktop/MD/dengue/data/HOLO/"
setwd(path)


pdb_holo <- read.pdb("only_protein.pdb") 
dcd_holo <- read.dcd("only_protein_1.dcd")

#print(pdb$xyz)
#print(dcd)

#######################################################################
################### NS2B ##############################################
#######################################################################


NS2B_inds <- atom.select(pdb_apo, segid = "BP1", elety = "CA")
NS2BHo_inds <- atom.select(pdb_holo, segid = "BP1", elety = "CA")

xyz_NS2B <- fit.xyz(fixed = pdb_apo$xyz, mobile = dcd_apo, fixed.inds = NS2B_inds$xyz, mobile.inds = NS2B_inds$xyz)
xyz_NS2BHo <- fit.xyz(fixed = pdb_holo$xyz, mobile = dcd_holo, fixed.inds = NS2BHo_inds$xyz, mobile.inds = NS2BHo_inds$xyz)
gc()

dim(xyz_NS2B) == dim(dcd_apo)

dim(xyz_NS2BHo) == dim(dcd_holo)

pc_NS2B <- pca.xyz(xyz_NS2B[,NS2B_inds$xyz])
pc_NS2BHo <- pca.xyz(xyz_NS2BHo[,NS2BHo_inds$xyz])

NS2B_pca <- as.data.frame(pc_NS2B$z)
NS2BHo_pca <- as.data.frame(pc_NS2BHo$z)

## select the PCA 1 and 2 

NS2B_pca <- NS2B_pca[, c("V1", "V2")]
NS2BHo_pca <- NS2BHo_pca[, c("V1", "V2")]

NS2B_pca["Protein"] <- rep(c("NS2B_apo"),c(length(NS2B_pca$V1)))
NS2BHo_pca ["Protein"] <- rep(c("NS2B_holo"), c((length(NS2BHo_pca$V1))))

all_protein_pca <- rbind(NS2B_pca, NS2BHo_pca)

pc_NS2B 
pc_NS2BHo

color <- c("NS2B_apo" = "#000000", "NS2B_holo" = "#07D511")

p_PCA_NS2B <- ggplot(all_protein_pca,aes(V1,V2, colour=Protein)) + theme_bw() + geom_point(size=.9) +
  theme(axis.text = element_text(family = "Roboto",size = 32), 
        axis.title = element_text( family = "Roboto",size = 34),
        legend.position = "none",
        plot.margin = margin(t = 10,  # Top margin
                             r = 40,  # Right margin
                             b = 10,  # Bottom margin
                             l = 10))+
  scale_colour_manual(values= color)  + 
  geom_vline(xintercept = 0,
             linetype="dashed", 
             color = "#C6C6C6", size = 1.5) +
  geom_hline(yintercept = 0,
             linetype="dashed",
             color = "#C6C6C6", size = 1.5) + 
  ylab("PC2 [13.45%, 19.82%]") +
  xlab("PC1 [36.52%, 25.21%]")

## Acumulative 

NS2B_relative <- as.data.frame(pc_NS2B$au[,1:2])
NS2BHo_relative <- as.data.frame(pc_NS2BHo$au[,1:2])

## getting residue 
## pdb_apo$atom$resno select all residues and NS2B_inds$atom select only the select previously  
NS2B_residues <- unique(pdb_apo$atom$resno[NS2B_inds$atom]) 
NS2B_relative["residue"] = NS2B_residues
NS2B_relative["Protein"] = rep(c("NS2B_apo"), c((length(NS2B_relative$V1))))

NS2BHo_residues <- unique(pdb_holo$atom$resno[NS2BHo_inds$atom]) 
NS2BHo_relative["residue"] = NS2BHo_residues
NS2BHo_relative["Protein"] = rep(c("NS2B_holo"), c((length(NS2BHo_relative$V1))))

all_relative_pca <- rbind(NS2B_relative, NS2BHo_relative)

p_acumulative_NS2B <- ggplot(all_relative_pca, aes(residue, V1, fill=Protein,colour=Protein)) + theme_bw() +
  ylab("Relative contribution to PC1") +
  xlab("Residue index of protein") +
  theme(axis.text = element_text(family = "Roboto",size = 18),
        axis.title = element_text( family = "Roboto",size = 24),
        legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 5,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5)) + 
  geom_line(size=0.8) + 
  scale_colour_manual(values= color) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(linewidth =7)))

ggsave("~/Desktop/MD/dengue/plots/PCA_NS2B.png", plot = p_PCA_NS2B, dpi = 300, width = 10, height = 8, units = "in") 

ggsave("~/Desktop/MD/dengue/plots/acumulative_NS2B.png", plot = p_acumulative_NS2B, dpi = 300, width = 10, height = 8, units = "in") 

#######################################################################
################### NS3 ##############################################
#######################################################################


NS3_inds <- atom.select(pdb_apo, segid = "AP1", elety = "CA")
NS3Ho_inds <- atom.select(pdb_holo, segid = "AP1", elety = "CA")

xyz_NS3 <- fit.xyz(fixed = pdb_apo$xyz, mobile = dcd_apo, fixed.inds = NS3_inds$xyz, mobile.inds = NS3_inds$xyz)
xyz_NS3Ho <- fit.xyz(fixed = pdb_holo$xyz, mobile = dcd_holo, fixed.inds = NS3Ho_inds$xyz, mobile.inds = NS3Ho_inds$xyz)
gc()

dim(xyz_NS3) == dim(dcd_apo)

dim(xyz_NS3Ho) == dim(dcd_holo)

pc_NS3 <- pca.xyz(xyz_NS3[,NS3_inds$xyz])
pc_NS3Ho <- pca.xyz(xyz_NS3Ho[,NS3Ho_inds$xyz])

NS3_pca <- as.data.frame(pc_NS3$z)
NS3Ho_pca <- as.data.frame(pc_NS3Ho$z)

## select the PCA 1 and 2 

NS3_pca <- NS3_pca[, c("V1", "V2")]
NS3Ho_pca <- NS3Ho_pca[, c("V1", "V2")]

NS3_pca["Protein"] <- rep(c("NS3_apo"),c(length(NS3_pca$V1)))
NS3Ho_pca ["Protein"] <- rep(c("NS3_holo"), c((length(NS3Ho_pca$V1))))

all_protein_pca_NS3 <- rbind(NS3_pca, NS3Ho_pca)

pc_NS3
pc_NS3Ho

color <- c("NS3_apo" = "#000000", "NS3_holo" = "#07D511")

p_PCA_NS3 <- ggplot(all_protein_pca_NS3,aes(V1,V2, colour=Protein)) + theme_bw() + geom_point(size=.9) +
  theme(axis.text = element_text(family = "Roboto",size = 32), 
        axis.title = element_text( family = "Roboto",size = 34),
        legend.position = "none",
        plot.margin = margin(t = 10,  # Top margin
                             r = 40,  # Right margin
                             b = 10,  # Bottom margin
                             l = 10))+
  scale_colour_manual(values=color)  + 
  geom_vline(xintercept = 0,
             linetype="dashed", 
             color = "#C6C6C6", size = 1.5) +
  geom_hline(yintercept = 0,
             linetype="dashed",
             color = "#C6C6C6", size = 1.5) + 
  ylab("PC2 [31.11%, 51.22%]") +
  xlab("PC1 [90.81%, 96.53%]")


## Acumulative 

NS3_relative <- as.data.frame(pc_NS3$au[,1:2])
NS3Ho_relative <- as.data.frame(pc_NS3Ho$au[,1:2])

## getting residue 
## pdb_apo$atom$resno select all residues and NS2B_inds$atom select only the select previously  
NS3_residues <- unique(pdb_apo$atom$resno[NS3_inds$atom]) 
NS3_relative["residue"] = NS3_residues
NS3_relative["Protein"] = rep(c("NS3_apo"), c((length(NS3_relative$V1))))

NS3Ho_residues <- unique(pdb_holo$atom$resno[NS3Ho_inds$atom]) 
NS3Ho_relative["residue"] = NS3Ho_residues
NS3Ho_relative["Protein"] = rep(c("NS3_holo"), c((length(NS3Ho_relative$V1))))

all_relative_pca_NS3 <- rbind(NS3_relative, NS3Ho_relative)

p_acumulative_NS3 <- ggplot(all_relative_pca_NS3, aes(residue, V1, fill=Protein,colour=Protein)) + theme_bw() +
  ylab("Relative contribution to PC1") +
  xlab("Residue index of protein") +
  theme(axis.text = element_text(family = "Roboto",size = 18),
        axis.title = element_text( family = "Roboto",size = 24),
        legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 5,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5)) + 
  geom_line(size=0.8) + 
  scale_colour_manual(values=color) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(linewidth =7)))

ggsave("~/Desktop/MD/dengue/plots/PCA_NS3.png", plot = p_PCA_NS3, dpi = 300, width = 10, height = 8, units = "in") 

ggsave("~/Desktop/MD/dengue/plots/acumulative_NS3.png", plot = p_acumulative_NS3, dpi = 300, width = 10, height = 8, units = "in") 

