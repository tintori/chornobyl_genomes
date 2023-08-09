# R code for figures in Tintori et al 2023 - Chernobyl O tipulae genomes paper
# Sophia Tintori, Rockman Lab, NYU
# July 2023

#~#~#
# FIGURE 1A - Map of collection sites
#~#~#

library(RColorBrewer) # v.1.1-3
RdYl.brewer <- rev(brewer.pal(11, "RdYlGn")[1:5])
library(ggplot2) # v.3.4.2
library(ggmap) # v.3.0.2
library(ggsn) # v.0.5.0
library(tidyverse) # v.2.0.0
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# Get map background
register_google(key = "AIzaSyD4jn5u4fOevgNsmR4FO2RVPnQQo107Ad4", write = TRUE)
tempMap10_bw <- get_map(location = "51.33748 29.84837", zoom = 10, maptype = "terrain", source = "google", color = "bw")

# read in "temp.finals" object from "../supplementary_data/01_CEZ_QG_table_for_paper.csv"
temp.finals = read.csv2("supplementary_data/01_CEZ_QG_table_for_paper.csv", sep = ",", stringsAsFactors = F)
for(tmp.col in c("GeigerAvg", "ambient.rad", "GPSlatitude", "GPSlatitude", "GPSelevation",
                 "TempCelsius", "PercentH2O", "BqCorrected", "BqPerGram")){
    temp.finals[,tmp.col] = as.numeric(temp.finals[,tmp.col])
}

temp.funnels = temp.finals$funnel
set.seed(3)
p <- ggmap(tempMap10_bw, extent = "device")+
    geom_jitter(data = subset(funnel_table, SampleName %in% temp.finals$funnel),
                aes(x=GPSlongitude, y=GPSlatitude,
                    colour=log10((GeigerAvg+.1)*8.76)), 
                width = .018, height = .01, size=3)+
    geom_point(aes(x=30.096789, y=51.389051), shape=4, size=0.5)+
    scale_color_gradientn(colours = RdYl.brewer, guide = "colourbar", limits = c(-.06,4), name="Ambient\nradiation",
                          breaks=c(0,1,2,3,4),labels=c("1 mSv/y", "10 mSv/y", "0.1 Sv/y", "1 Sv/y", "10 Sv/y"))+
    scalebar(x.min = 29.55476, x.max = 30.20199,
             y.min = 51.18535, y.max = 51.3896,
             dist = 5, dist_unit = "mi",
             box.fill = c("gray", "white"), box.color = "dark gray",
             st.color = "gray",
             height = .06, st.dist = .05, border.size = .075, st.size = 2.5,
             st.bottom = FALSE,
             transform = TRUE)
p

ggsave(paste("plots_out/01A_298_samples_on_map.pdf", sep = ""),  device = "pdf",
       scale = 1, width = 7, height = 7, units = "in", dpi = 600)

#~#~#
# FIGURE 1B - Gene tree of 298 wild isolates 
#~#~#

library(ape) # v.5.7-1
library(tidyverse) # v.2.0.0
library(msa) # v.1.28.0
library(ggtree) # v.3.6.2
source("supplementary_data/00_chornobyl_plot_shortcuts.R")


# Read in input.merged.table from supp file "../supplementary_data/02_table_for_tree.csv"
input.merged.table = read.csv2("supplementary_data/02_table_for_tree.csv", sep = ",", stringsAsFactors = F)
tree.pieces=list()
tree.pieces[["DNAseqs"]] = data.frame(cez.id = character(0), QG.id = character(0), species.guess = character(0), 
                                      clade = character(0), ambient.rad = numeric(0), seq = character(0), direction = character(0))
for(rown in 1:nrow(input.merged.table)){
    if(input.merged.table$fasta.merged[rown]!=""){
        tree.pieces$DNAseqs[nrow(tree.pieces$DNAseqs)+1,] <- unlist(c(input.merged.table[rown, c("LineName", "QG.id", "species.guess", "clade", "ambient.rad", "fasta.merged")], "merge"))
    } else if(input.merged.table$fasta.f[rown]!=""){
        if(input.merged.table$fasta.r[rown]==""){
            tree.pieces$DNAseqs[dim(tree.pieces$DNAseqs)[1]+1,] <- unlist(c(input.merged.table[rown, c("LineName", "QG.id", "species.guess", "clade", "ambient.rad", "fasta.f")], "forward"))
        } else {
            temp.nchar.f <- nchar(input.merged.table$fasta.f[rown])
            temp.nchar.r <- nchar(input.merged.table$fasta.r[rown])
            if(temp.nchar.f < temp.nchar.r) {tree.pieces$DNAseqs[dim(tree.pieces$DNAseqs)[1]+1,] <- unlist(c(input.merged.table[rown, c("LineName", "QG.id", "species.guess", "clade", "ambient.rad", "fasta.r")], "reverse"))}
            else {tree.pieces$DNAseqs[dim(tree.pieces$DNAseqs)[1]+1,] <- unlist(c(input.merged.table[rown, c("LineName", "QG.id", "species.guess", "clade", "ambient.rad", "fasta.f")], "forward"))}
        }
    } else if (input.merged.table$fasta.r[rown]!=""){
        tree.pieces$DNAseqs[dim(tree.pieces$DNAseqs)[1]+1,] <- unlist(c(input.merged.table[rown, c("LineName", "QG.id", "species.guess", "clade", "ambient.rad", "fasta.r")], "reverse"))
    }
}

# Add plectus as outgroup
outgroup.species = "Plectus minimus"
outgroup.genus <- strsplit(outgroup.species, " ")[[1]][1]

tree.pieces$DNAseqs <- rbind(tree.pieces$DNAseqs, c("Plectus", "Plectus", "Pectus minimus", "Plectus", 0,
                                                    "ggaaggcagcaggcgcgcaaattacccactctcggcacgaggaggtagtgacgaaaaataacgaggcggttctctatgaggcccgctatcggaatgggtacaatttaaaccctttaacgaggacctatgagagggcaagtctggtgccagcagccgcggtaattccagctctcaaggtgtatatcgccattgctgcggttaaaaagctcgtagttggatctgcgccttcggactcggtccgcccaacgggtgtgaactgggatccaaggcttatactgctggttttcccttgatgctctttactgggtgtcttgggtggctagcgagtttactttgaaaaaattagagtgcttaacacaggctaacgcctgaatactcgtgcatggaataatagaataagaccacggttctattttattggttttcggaactgtgataatggttaagagggacagacgggggcattcgtatcgctgcgtgagaggtgaaattcttggaccgcagcgagacgccctactgcgaaagcatttgccaagaatgtcttcgttaatcaagaacgaaagtcagagg", 
                                                    "merge"))

tree.pieces$DNAseqs$ambient.rad <- as.numeric(tree.pieces$DNAseqs$ambient.rad)

# set up input sequences for msa
tree.pieces[["fasta"]] <- tree.pieces$DNAseqs$seq
names(tree.pieces$fasta) <- tree.pieces$DNAseqs$cez.id
# turn those sequences into a DNAStringSet object
tree.pieces[["DNAStringSet"]] <- DNAStringSet(tree.pieces$fasta)

# make the alignment before calculating distances
print("Performing MSA (this might take five minutes)")
tree.pieces[["alignment"]] <- msa(tree.pieces$DNAStringSet)

# calculate distance
tree.pieces[["converted.alignment.ape"]] <- msaConvert(tree.pieces$alignment, type="ape::DNAbin")
#tree.pieces[["converted.alignment.seqinr"]] <- msaConvert(tree.pieces$alignment, type="seqinr::alignment")

print("Computing distances")
tree.pieces[["dist.align.ape"]] <- dist.dna(tree.pieces$converted.alignment.ape, 
                                             model = "raw", pairwise.deletion = T) # other models leave NAs in distances, which then make all branch lengths get thrown out

# make black and white tree
dist.tree = njs(tree.pieces$dist.align.ape)
dist.tree = root(dist.tree, outgroup = "Plectus")
dist.tree$clade = dist.tree$tip.label
dist.tree$amb.rad = dist.tree$tip.label
for(i in 1:length(dist.tree$clade)){
    tmp.samp = dist.tree$clade[i]
    if(!tmp.samp %in% input.merged.table$LineName){next}
    tmp.clade = input.merged.table$clade[which(input.merged.table$LineName==tmp.samp)]
    dist.tree$clade[i] = tmp.clade
    tmp.amb.rad = input.merged.table$ambient.rad[which(input.merged.table$LineName==tmp.samp)]
    dist.tree$amb.rad[i] = tmp.amb.rad
}
dist.tree$amb.rad = as.numeric(dist.tree$amb.rad)
dist.tree$edge.length = abs(dist.tree$edge.length)
pdf("plots_out/01B_tree_w_clades.pdf", width = 5, height = 15)
plot(dist.tree, cex=.1, show.tip=F)
tiplabels(dist.tree$clade, cex=.2, frame = "none", adj = 0, offset = .005)
axisPhylo()
dev.off()

# make color tree
require(RColorBrewer)
RdYl.brewer <- rev(brewer.pal(11, "RdYlGn")[1:5])
rad_color_scale = scale_color_gradientn(colours = RdYl.brewer, guide = "colourbar", limits = c(0,4), name="Ambient\nradiation",
                                        breaks=c(0,1,2,3,4),labels=c("1 mSv/y", "10 mSv/y", "0.1 Sv/y", "1 Sv/y", "10 Sv/y"))
rad_fill_scale = scale_fill_gradientn(colours = RdYl.brewer, guide = "colourbar", limits = c(0,4), name="Ambient\nradiation",
                                      breaks=c(0,1,2,3,4),labels=c("1 mSv/y", "10 mSv/y", "0.1 Sv/y", "1 Sv/y", "10 Sv/y"))
dd <-data.frame(taxa=c(dist.tree$tip.label),amb.rad=c(dist.tree$amb.rad), clade = dist.tree$clade)

ggtree(dist.tree) %<+% dd + 
    geom_tippoint(aes(color=log10(amb.rad*8.76)))+
    rad_color_scale+
    geom_treescale()
ggsave("plots_out/01B_tree_w_colors.pdf", device = "pdf", height = 15, width = 5)

ggtree(dist.tree) %<+% dd + 
    geom_tiplab(aes(label=clade, color=log10(amb.rad*8.76)), size=2)+
    rad_color_scale+
    geom_treescale()
ggsave("plots_out/01B_tree_w_color_clades.pdf", device = "pdf", height = 15, width = 5)


#~#~#
# FIGURE 1C
#~#~#

library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# using same table as for Figure 1A (temp.finals, from "../supplementary_data/01_CEZ_QG_table_for_paper.csv")

temp.finals %>%
    filter(clade %in% c("Oscheius", "Panagrolaimus", "Acrobeloides"), SampleType!="") %>%
    ggplot(aes(x=as.numeric(PercentH2O), y=as.numeric(GeigerAvg)*8.76, color=SampleType))+
    xlab(expression(Moisture~("%"~weight~H[2]*O)))+
    ylab("Ambient radiation (mSv/yr)")+
    labs(color="Source\nof worms")+
    geom_point()+
    theme_bw()+
    coord_trans(y = "log10")+
    facet_wrap(~clade, ncol = 3)+
    theme(
        strip.background = element_blank(),
        panel.grid.minor = element_blank()
        #strip.text.x = element_blank()
    )+
    scale_y_continuous(breaks=c(1, 10,100, 1000,10000), limits = c(1,10000))
ggsave("plots_out/01C_environmental_measurements.pdf", 
       width = 8, height = 3)


#~#~#
# Supplementary Figure 1 - just the little orange dots 
#~#~#

library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

temp.finals %>%
    filter(cryo.status==T) %>% select(GeigerAvg) %>%
    mutate(log10_msv_yr = log10(as.numeric(GeigerAvg)*8.76)) %>%
    ggplot(aes(x=log10_msv_yr, y=1, color = log10_msv_yr))+
    rad_color_scale+
    geom_jitter(height = 1, width = .02, size=4.5)+
    theme_sophie
ggsave("plots_out/Supp01_298_dots.pdf", 
       device = "pdf", width = 9, height = 2) 

#~#~#
# Supplementary Figure 2
#~#~#   

library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# using same table as for Figure 1A (temp.finals, from "../supplementary_data/01_CEZ_QG_table_for_paper.csv")

ggplot(data=subset(temp.finals, (!is.na(BqCorrected) & BqCorrected >0)), 
       aes(x=log10(ambient.rad), y = log10(BqCorrected), color = SampleType))+
    xlab("Ambient equivalent dose\nlog10( mSv/hr )") + ylab("Substrate radioactivity\nlog10( Becquerel )") +
    labs(color="Source\nof worms")+
    scale_y_continuous(limits = c(-1,5))+
    scale_x_continuous(limits = c(-1,3))+
    geom_point()+
    theme_bw()+
    stat_smooth(geom='line', alpha=0.5, se=FALSE, method = "lm")

ggsave(filename = "plots_out/Supp02_Bq_vs_Sv.pdf", 
       width = 6, height = 5)

#~#~#
# Supplementary Figure 3
#~#~#

library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# using same table as for Figure 1A (temp.finals, from "../supplementary_data/01_CEZ_QG_table_for_paper.csv")

register_google(key = "AIzaSyD4jn5u4fOevgNsmR4FO2RVPnQQo107Ad4", write = TRUE)
tempMap10_bw <- get_map(location = "51.33748 29.84837", zoom = 10, maptype = "terrain", source = "google", color = "bw")

p = temp.finals %>% 
    filter(clade %in% c("Oscheius", "Panagrolaimus", "Acrobeloides"), SampleType!="") %>%
    ggmap(tempMap10_bw, extent = "device")+
    geom_jitter(data = subset(funnel_table, SampleName %in% temp.finals$funnel),
                aes(x=GPSlongitude, y=GPSlatitude,
                    colour=log10((GeigerAvg+.1)*8.76)), 
                width = .018, height = .01, size=3)+
    geom_point(aes(x=30.096789, y=51.389051), shape=4, size=0.5)+
    facet_wrap(~clade, nrow = 1)+
    scale_color_gradientn(colours = RdYl.brewer, guide = "colourbar", limits = c(-.06,4), name="Ambient\nradiation",
                          breaks=c(0,1,2,3,4),labels=c("1 mSv/y", "10 mSv/y", "0.1 Sv/y", "1 Sv/y", "10 Sv/y"))+
    scalebar(x.min = 29.55476, x.max = 30.20199,
             y.min = 51.18535, y.max = 51.3896,
             dist = 5, dist_unit = "mi",
             box.fill = c("gray", "white"), box.color = "dark gray",
             st.color = "gray",
             height = .06, st.dist = .05, border.size = .075, st.size = 2.5,
             st.bottom = FALSE,
             transform = TRUE)
p

ggsave("plots_out/Supp03_geog_of_3_clades.pdf",
              scale = 1, width = 7, height = 7, units = "in", dpi = 600)

#~#~#
# Figure 2A
#~#~#

library(tidyverse)
library(reshape2)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

diff_counts_anyDP05 = read.csv2("supplementary_data/03_the_distance_table.csv")

# cluster to order axes
actual_clust = diff_counts_anyDP05 %>% 
    dcast(samp1~samp2, value.var = "freq_of_var_any") %>%
    column_to_rownames(var = "samp1") %>% 
    as.dist() %>% hclust()

# make blue plot
diff_counts_anyDP05 %>%
    ggplot(aes(x=factor(samp1, levels = actual_clust$labels[actual_clust$order]),
               y=factor(samp2, levels = actual_clust$labels[actual_clust$order]),
               fill=freq_of_var_any))+
    geom_tile(color="#FFFFFF", linewidth=.005)+
    coord_fixed()+
    theme_sophie+
    labs(x=NULL, y=NULL)+
    theme(legend.title = element_blank())+
    scale_fill_gradient(low="#FFFFFF", high="#08306b")+
    theme(axis.text.x = element_text(face="bold", colour = rad_colors[actual_clust$labels[actual_clust$order]]))+
    theme(axis.text.y = element_text(face="bold", colour = rad_colors[actual_clust$labels[actual_clust$order]]))

ggsave("plots_out/02A_blue_plot.pdf", 
       device = "pdf", width = 8, height = 6)

#~#~#
# Figures 2B + 2C
#~#~#

library(tidyverse)
library(ggrepel)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# On all sites where there is >5 coverage for everyone, and exactly 2 alleles
ju75_lcc_DP05 = read.csv2("supplementary_data/04_vcfs/all_JU75_oriented_lcc_anyDP05_m2.vcf", 
                          comment.char = "#", sep = "\t", header = F, stringsAsFactors = F,  # this takes one minute to load
                          col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                                        "FORMAT", "CEZ.010", "CEZ.047", "CEZ.061", "CEZ.077", "CEZ.083",
                                        "CEZ.107", "CEZ.119", "CEZ.138", "CEZ.144", "CEZ.156", "CEZ.166",
                                        "CEZ.224", "CEZ.253", "CEZ.263", "CEZ.274", "DF5062",  "DF5103",
                                        "DF5110",  "DF5111",  "DF5123"))
ju75_lcc_DP05$INFO = ""
ju75_lcc_DP05 = ju75_lcc_DP05  %>% filter(!grepl(pattern = ",", x = ALT))
dot.samps=intersect(rad_table$dot_names, colnames(ju75_lcc_DP05))
for(tmp.col in dot.samps){  # This will take ~20-25 minutes to process
    print(Sys.time())
    print(tmp.col)
    tmp.vect = ju75_lcc_DP05 %>% 
        select(tmp.col) %>%
        separate(col = tmp.col, into = c("gt", "ad", "rest1"), sep = ":", extra = "merge") %>%
        select(ad) %>%
        separate(col="ad", into = c("ref.d", "alt.d"), sep = ",", extra = "warn") %>%
        mutate(frac_alt = (as.integer(alt.d)/(as.integer(alt.d)+as.integer(ref.d)))) %>%
        mutate(gt = case_when(as.integer(alt.d)+as.integer(ref.d) < 5 ~ NA,
                              frac_alt>.6 ~ 1,
                              frac_alt<.4 ~ 0,
                              .default=0.5 ))
    ju75_lcc_DP05[,tmp.col] = tmp.vect[,"gt"]
}

ju75_lcc_DP05$total = apply(ju75_lcc_DP05[,dot.samps],
                            1, sum)
table(ju75_lcc_DP05$total)
ju75_lcc_DP05 = ju75_lcc_DP05 %>%
    filter(!total %in% c(0,20))
dim(ju75_lcc_DP05)
# down to 21K loci # hmm, this time i'm getting 609K loci. is that because it's any and not all?
ju75_lcc_DP05 %>% na.omit() %>% dim() # goes down to 29K, ok thats better but still a little weird

# now pca
corr_matrix <- ju75_lcc_DP05[,dot.samps] %>% na.omit() %>% cor()
data.pca <- princomp(corr_matrix)
data.pca$toplot = data.pca$scores %>%
    as.data.frame() %>%
    mutate(dot_names = rownames(data.pca$scores)) %>%
    left_join(rad_table)

ggplot(data.pca$toplot, aes(x=Comp.1, y=Comp.2, color=log10_mSv_per_yr))+
    geom_point(size=3) + 
    rad_color_scale+
    theme_sophie+
    geom_text_repel(data=filter(data.pca$toplot, is.na(log10_mSv_per_yr)),
                    aes(label=fancy_names))+
    theme(legend.position = "none")
ggsave("plots_out/02B_PCA_wide.pdf", 
       device = "pdf", width = 5, height = 5)

ggplot(data.pca$toplot, aes(x=Comp.1, y=Comp.2, color=log10_mSv_per_yr))+
    geom_point(size=3) + 
    rad_color_scale+
    theme_sophie+
    geom_text_repel(data=data.pca$toplot,aes(label=fancy_names))+
    xlim(c(0.15,0.45))+ylim(c(-.025,.035))+
    theme(legend.position = "none")
ggsave("plots_out/02C_PCA_close.pdf", 
       device = "pdf", width = 5, height = 5)

#~#~#
# Mantel test (geographic distance vs genetic distance)
#~#~#

library(geosphere)
library(tidyverse)
library(adegenet)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# Load "diff_counts_anyDP05" from "../supplementary_data/03_the_distance_table.csv"

geog_dist = data.frame()
for(i in unique(diff_counts_anyDP05$samp1)){
    for(j in unique(diff_counts_anyDP05$samp2)){
        print(paste(i,j))
        geog_dist[nrow(geog_dist)+1, c("samp1", "samp2", "haversine_dist")] = 
            c(i, j, distHaversine(p1 = rad_table[i,c("long", "lat")], p2 = rad_table[j,c("long", "lat")]))
    }
}

geog_dist.dist = geog_dist %>%
    arrange(samp1, samp2) %>%
    dcast(samp1~samp2, value.var="haversine_dist")  %>%
    column_to_rownames(var="samp1") %>%
    dist()
genet_dist = diff_counts_anyDP05 %>% 
    arrange(samp1, samp2) %>%
    select(samp1, samp2, freq_of_var_any) %>%
    dcast(samp1~samp2, value.var="freq_of_var_any") %>%
    column_to_rownames(var="samp1") %>% 
    dist()
mantel.randtest(m1 = geog_dist.dist, m2 = genet_dist, nrepet = 999)

#~#~#
# Supplementary figure 4  
#~#~#

rm(cew1_M2, tmp.cew1_M2.chrom, tmp.vect, tree.pieces, ju75_lcc_DP05)
rm(list = ls(pattern="^tmp|^temp|^test"))

library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# Load "diff_counts_anyDP05" from "../supplementary_data/03_the_distance_table.csv"

full_join(geog_dist, diff_counts_anyDP05) %>% 
    mutate(comp_type = case_when(samp1 %in% CEZ_samples & samp2 %in% CEZ_samples ~ "cez_to_cez",
                                 !samp1 %in% CEZ_samples & !samp2 %in% CEZ_samples ~ "non_to_non",
                                 .default = "mix")) %>%
    ggplot()+
    geom_point(aes(x=as.numeric(haversine_dist)/1000, y=freq_of_var_any, fill=comp_type), shape=21, color="black")+
    scale_fill_manual(breaks = c("cez_to_cez", "non_to_non", "mix"),
                       values = c("#a80000", "#494949", "#9b6300"),
                       name="")+
    labs(x="geographic distance (km)", y = "genetic distance (var per bp)")+
    ggtitle("all samples")+
    theme_sophie
ggsave("plots_out/Supp04A_mantel_all_samples.pdf",
       device = "pdf", width = 5, height = 4)

full_join(geog_dist, diff_counts_anyDP05) %>%  
    filter(samp1 %in% CEZ_samples, samp2 %in% CEZ_samples) %>%
    mutate(comp_type = case_when(samp1 %in% CEZ_samples & samp2 %in% CEZ_samples ~ "cez_to_cez",
                                 !samp1 %in% CEZ_samples & !samp2 %in% CEZ_samples ~ "non_to_non",
                                 .default = "mix")) %>%
    ggplot()+
    geom_point(aes(x=as.numeric(haversine_dist)/1000, 
                   y=freq_of_var_any), fill="#a80000", shape=21, color="black")+
    labs(x="geographic distance (km)", y = "genetic distance (var per bp)")+
    ggtitle("cez samples")+
    theme_sophie
ggsave("plots_out/Supp04B_mantel_cez_samples.pdf",
       device = "pdf", width = 5, height = 4)

#~#~#
# Supplementary figure 5
#~#~#

library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

load("supplementary_data/04_vcfs/cew1_big_vcf_DP05_lcc_M2.RData")

cew1_M2 = cew1_M2 %>% 
    mutate(total = CEZ.010+CEZ.047+CEZ.061+CEZ.077+CEZ.083+CEZ.107+CEZ.119+CEZ.138+CEZ.144+CEZ.156+CEZ.166+CEZ.224+CEZ.253+CEZ.263+CEZ.274+DF5062+DF5103+DF5110+DF5111+DF5123)

table(cew1_M2$total)
#  0        ...        20
# 48292352   1290652      58427

cew1_M2 = cew1_M2 %>%
    mutate(aligned = 1, 
           var = case_when(!total %in% c(0,20) ~ 1))

# Make 10 kb windows to tally everything up 

# find the subset of sites to calculate
cew1_M2 = cew1_M2 %>% 
    mutate(binpos10kb = paste(CHROM, as.character(10000*floor(POS/10000)), sep = "_"))
subst10kb = match(unique(cew1_M2$binpos10kb), cew1_M2$binpos10kb)
length(subst10kb) # 5898
cew1_M2 = cew1_M2 %>% select(-binpos10kb)

for(tmp.chrom in unique(cew1_M2$CHROM)){   # Might take a half an hour per chrom, might have to increase memory allocation
    print(paste(Sys.time(), tmp.chrom, sep = " ~ "))
    cew1_mini = cew1_M2[subst10kb,] %>%   # a table of just the sites that need calculating
        filter(CHROM==tmp.chrom) %>%
        select(CHROM, POS, aligned, var)
    tmp.cew1_M2.chrom = cew1_M2 %>% filter(CHROM==tmp.chrom) # a table of all sites in the chrom
    print(dim(tmp.cew1_M2.chrom))
    print(dim(cew1_mini))
    
    for(rown in 1:nrow(cew1_mini)){
        if(rown%%25==0){
            print(rown)
            if(rown%%200==0){Sys.sleep(1)}}
        tmp.pos = cew1_mini[rown,"POS"]
        tmp.chrom = cew1_mini[rown,"CHROM"]
        tmp.all.aligned_1MB = sum(tmp.cew1_M2.chrom %>% 
                                      filter(CHROM==tmp.chrom, POS>(tmp.pos-500000), POS<(tmp.pos+500000)) %>%
                                      select(aligned), na.rm = T)
        
        tmp.all.var_1MB = sum(tmp.cew1_M2.chrom %>% 
                                  filter(CHROM==tmp.chrom, POS>(tmp.pos-500000), POS<(tmp.pos+500000)) %>%
                                  select(var), na.rm = T)
        
        cew1_mini[rown, "aligns_and_vars_in_10kb_1mb"] = paste(c(tmp.all.aligned_1MB, tmp.all.var_1MB), collapse = ",")
    }
    print("ok")
    Sys.sleep(2)
    cew1_M2[which(cew1_M2$CHROM==tmp.chrom & cew1_M2$POS %in% cew1_mini$POS),"a_and_v_10kb_1mb"] = 
        cew1_mini$aligns_and_vars_in_10kb_1mb
    rm(cew1_mini, tmp.cew1_M2.chrom)
    gc()
    Sys.sleep(2)
}

cew1_M2_10kb = cew1_M2 %>%
    filter(!is.na(a_and_v_10kb_1mb)) %>%
    separate(col = a_and_v_10kb_1mb, into = c("alignscount_1mb", "varscount_1mb"), sep = ",", convert = T)

cew1_M2_10kb %>%
    filter(CHROM!="mito") %>%
    ggplot(aes(x=POS, y=varscount_1mb/alignscount_1mb))+
    geom_line()+
    facet_grid(.~CHROM, scales = "free_x",space = "free_x")+
    theme_sophie
ggsave("plots_out/Supp05_cew1_distr_of_vars.pdf", 
       device = "pdf", width = 12, height = 4)

#~#~#
# Supplementary figure 6
#~#~#

library(tidyverse)
library(scales)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# intrachromosomal ld
test.ld_cew1 = read.csv2("supplementary_data/05_CEW1_full_anyDP05.ld", 
                         sep="", header = T, stringsAsFactors = F)
test.ld_cew1$CHR_A[which(test.ld_cew1$CHR_A=="23")] = "chrX"
test.ld_cew1$CHR_B[which(test.ld_cew1$CHR_B=="23")] = "chrX"
ref.chrom.levels = c("chrI","chrII","chrIII","chrIV",   
                     "chrV","chrX","mito")
test.ld_cew1$CHR_A = factor(test.ld_cew1$CHR_A, levels = ref.chrom.levels)
test.ld_cew1$CHR_B = factor(test.ld_cew1$CHR_B, levels = ref.chrom.levels)
test.ld_cew1$R2 = as.numeric(test.ld_cew1$R2)
test.ld_cew1 = test.ld_cew1 %>%
    mutate(distance = case_when(CHR_A==CHR_B ~ abs(BP_B - BP_A),
                                CHR_A!=CHR_B ~ BP_B))
avg.ld_1 = test.ld_cew1 %>% 
    filter(distance<2000) %>% 
    group_by(distance) %>%
    mutate(all_chrom_avg_ld = mean(R2)) %>%  
    ungroup() %>%
    group_by(distance, CHR_A, CHR_B, all_chrom_avg_ld) %>%
    summarise(avg_ld = mean(R2))
highest = 200
test.ld_cew1 %>% 
    sample_n(10000) %>% 
    filter(distance < highest) %>% 
    filter(CHR_B != "mito") %>%  
    ggplot(aes(x=distance, y = R2))+
    geom_point(alpha=.20, size=.75) +
    labs(x="distance (bp)", y = bquote(expr = 'LD ('~R^2~')'))+
    theme_sophie+
    scale_x_continuous(breaks = pretty_breaks())+
    geom_line(data = avg.ld_1 %>% filter(distance<highest, CHR_A == "chrI") %>% select(distance, all_chrom_avg_ld) %>% distinct(), 
              aes(x=distance, y=all_chrom_avg_ld), color="#FF0000", linewidth=2)
ggsave("plots_out/Supp06A_ld_200_1_all_genome.pdf", 
       device="pdf", width=5, height = 4)

# interchromosomal ld
interchr.ld = read.csv2("supplementary_data/06_CEW1_full_anyDP05_interchr001.ld", 
                         sep="", header = T, stringsAsFactors = F)
interchr.ld$CHR_A[which(interchr.ld$CHR_A=="23")] = "chrX"
interchr.ld$CHR_B[which(interchr.ld$CHR_B=="23")] = "chrX"

ref.chrom.levels = c("chrI","chrII","chrIII","chrIV",   
                     "chrV","chrX","mito")
interchr.ld$CHR_A = factor(interchr.ld$CHR_A, levels = ref.chrom.levels)
interchr.ld$CHR_B = factor(interchr.ld$CHR_B, levels = ref.chrom.levels)
interchr.ld$R2 = as.numeric(interchr.ld$R2)
interchr.ld = interchr.ld %>%
    mutate(distance = case_when(CHR_A==CHR_B ~ abs(BP_B - BP_A),
                                CHR_A!=CHR_B ~ BP_B))
dim(interchr.ld) # 2 mill rows

interchr.ld = data.frame(CHR_A = c(interchr.ld$CHR_A, interchr.ld$CHR_B),
                         CHR_B = c(interchr.ld$CHR_B, interchr.ld$CHR_A),
                         R2 = c(interchr.ld$R2, interchr.ld$r2),
                         distance = c(interchr.ld$distance, interchr.ld$distance))

interchr.ld.avg = interchr.ld %>%
    filter(CHR_A!=CHR_B) %>%
    group_by(CHR_A) %>%
    mutate(all_chrom_avg_r2 = mean(R2)) %>%
    ungroup() %>%
    group_by(CHR_A, CHR_B, all_chrom_avg_r2) %>%
    summarise(avg_r2 = mean(R2)) 

interchr.ld.avg %>% 
    ungroup %>%
    select(CHR_A, all_chrom_avg_r2) %>% distinct() %>% 
    ggplot(aes(x=CHR_A, y=all_chrom_avg_r2)) +
    geom_point()+
    ylim(c(0,1))+
    labs(x="Chromosome", 
         y=bquote(expr = 'Avg LD with other chromosomes ('~R^2~')'))+
    theme_sophie
ggsave("plots_out/Supp06B_ld_interchrom_avg.pdf", 
       device="pdf", width=5, height = 5)

#~#~#
# Figure 3
#~#~#

# First, made a bunch of pafs (alignment files) using minimap2 v2.17 [minimap2 -x sr -t 4 cew1.fa new_genome.fa -o new_genome.paf]

library(pafr)
library(ggplot2)
library(seqinr)
library(Biostrings)
library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

dir.create("plots_out/03_syntenies")
for(temp.file in list.files("supplementary_data/07_all_pafs_to_CEW1/")){
    temp.samp=strsplit(temp.file, "\\.paf")[[1]][1]
    print(temp.samp)
    if(nchar(temp.samp)>7){next}
    
    # made this paf with 00_scripts/55_align_for_pafs.sh
    paf_in = paste0("supplementary_data/07_all_pafs_to_CEW1/", temp.samp, ".paf")
    ali <- read_paf(paf_in)
    
    # dot plot
    prim_alignment <- filter_secondary_alignments(ali)
    rorder=intersect(c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "mito"), prim_alignment$tname)
    qorder=unique(prim_alignment$qname)
    dotplot(prim_alignment, order_by="provided", ordering = list(qorder, rorder), 
            xlab = temp.samp, ylab = "CEW1", dashes = T)+
        theme_sophie+
        theme(panel.border = element_rect(colour=rad_colors[temp.samp], linewidth=5),
              axis.title = element_blank(),
              axis.text.x=element_blank(), axis.text.y=element_blank())
    ggsave(paste0("plots_out/03_syntenies/", temp.samp, ".pdf"), 
           device = "pdf", width = 2.5, height = 2.5)
}

#~#~#
# Supplementary figure 7
#~#~#

require(seqinr)
require(reshape2)
library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

contiguity.table = data.frame()
depths = read.csv("supplementary_data/08_seq-depth_etc.txt", 
                  comment.char = "#", sep = "\t", stringsAsFactors = F)
for(tmp.samp in unique(rad_table$sample)){
    print(tmp.samp)
    if(tmp.samp %in% c("CEW1", "JU75", "CEZ-164", "JU170")){next}
    tmp.fa = read.fasta(paste0("supplementary_data/09_genomes/", tmp.samp, ".fa"))
    tmp.numb.cont = unlist(lapply(tmp.fa, length))
    tmp.numb.long.cont = length(which(tmp.numb.cont>200000))
    tmp.numb.cont = length(tmp.numb.cont)
    tmp.illum = depths[which(depths$platform=="illumina" & depths$sample==tmp.samp),"xCoverage"]
    tmp.min = depths[which(depths$platform=="minion" & depths$sample==tmp.samp),"xCoverage"]
    contiguity.table[tmp.samp, c("sample", "numb.contigs", "numb.long.contigs", "illum.depth", "min.depth")] = 
        c(tmp.samp, tmp.numb.cont, tmp.numb.long.cont, tmp.illum, tmp.min)
}
contiguity.table$illum.depth = as.numeric(contiguity.table$illum.depth)
contiguity.table$min.depth = as.numeric(contiguity.table$min.depth)
contiguity.table$numb.contigs = as.numeric(contiguity.table$numb.contigs)
contiguity.table$numb.long.contigs = as.numeric(contiguity.table$numb.long.contigs)


contiguity.table %>% 
    melt(id=c("sample", "numb.long.contigs", "numb.contigs"), variable.name = "platform", value.name = "depth") %>%
    ggplot(aes(x=numb.long.contigs, y=depth, shape=platform, color=sample))+
    geom_point() + 
    labs(x="Number of contigs over 200K", y="Avg read depth")+
    theme_sophie +
  #  theme(legend.position = "none")
    theme(legend.direction = "vertical", legend.position = "bottom")

ggsave("plots_out/Supp07_depth_vs_long_contiguity.pdf", # "~/Dropbox/postdoc/023_genome_seq/24_depth_vs_long_contiguity.pdf"
       device = "pdf", height = 8, width = 4)

#~#~#
# Supplementary Figure 8
#~#~#

library(pafr)
library(seqinr)
library(Biostrings)
library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

samps=c("CEZ-263", "CEZ-107")
refs=c("CEZ-107", "CEZ-263")
for(i in 1:length(samps)){
    temp.samp=samps[i] # for second plot: "CEZ-107"
    temp.ref=refs[i] # for second plot: "CEZ-263"
    paf_in = paste0("supplementary_data/10_pafs_denovo_pairs/", temp.ref, "_", temp.samp, ".paf")
    ali <- read_paf(paf_in)
    
    # dot plot
    prim_alignment <- filter_secondary_alignments(ali)
    qorder=unique(unlist(c(as.vector(read.table(paste0("supplementary_data/10_pafs_denovo_pairs/", temp.samp, "_contig_order.txt"))),
                           unique(prim_alignment$qname))))
    rorder=unique(unlist(c(as.vector(read.table(paste0("supplementary_data/10_pafs_denovo_pairs/", temp.ref, "_contig_order.txt"))),
                           unique(prim_alignment$tname))))
    
    dotplot(prim_alignment, # %>% filter(tname %in% rorder, qname %in% qorder), 
            order_by="provided", ordering = list(qorder, rorder), 
            xlab = temp.samp, ylab = temp.ref)+
        theme_sophie+
        theme(panel.border = element_rect(colour=rad_colors[temp.samp], size=10))
    ggsave(paste0("plots_out/Supp08_nonref_alignments_", temp.ref, "_", temp.samp, ".pdf"), 
           device = "pdf", width = 5, height = 5)
}


#~#~#
# Supplementary figure 9
#~#~#

library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

samps=c("CEZ-263", "CEZ-107")
contigs=c("contig_2_pilon", "contig_17_pilon")
starts=c(2033000, 2115000)
ends=c(2048000, 2120000)
for(i in 1:length(samps)){
    temp.samp=samps[i]
    paf_in = paste0("supplementary_data/10_pafs_denovo_pairs/", temp.samp, ".paf")
    ali <- read_paf(paf_in)
    prim_alignment <- filter_secondary_alignments(ali)
    ali = prim_alignment %>% 
        filter(qname==contigs[i], qstart>starts[i], qend<ends[i]) %>% 
        left_join(CEW1_contigs %>% select(contigs, chroms_mid), 
                  by = c("tname" = "contigs"))
    line_size=5
    ggplot()+
        geom_segment(data = ali[ali$strand == "+",], 
                     aes_string(x = "qstart", xend = "qend", 
                                y = "tstart", yend = "tend"), size = line_size #, colour = alignment_colour
        ) + 
        geom_segment(data = ali[ali$strand == "-", ], 
                     aes_string(x = "qstart", xend = "qend", 
                                y = "tstart", yend = "tend"), size = line_size #, colour = alignment_colour
        ) +
        scale_x_continuous(paste("Position on ", temp.samp, " chrV"), labels = Mb_lab) + 
        scale_y_continuous("Position on CEW1 chromosomes", labels = Mb_lab) +
        facet_wrap(~chroms_mid, ncol = 1)+
        theme_sophie
    ggsave(paste0("plots_out/Supp09A_", temp.samp, ".pdf"), 
           device = "pdf", width = 5, height = 5)
}

# Bottom half of figure is a screengrab from genome workbench

#~#~#
# Figure 4B
#~#~#

library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# Read in mut_type_tally_DF5123out from "../supplementary_data/10_DF5123_mut_type_tally.csv"
mut_type_tally_DF5123out = read.csv2("supplementary_data/11_DF5123_mut_type_tally.csv", sep = ",")
for(tmp.col in c("log10_mSv_y_str1", "log10_mSv_y_str2", "Bq", "percH20")){
    mut_type_tally_DF5123out[,tmp.col] = as.numeric(mut_type_tally_DF5123out[,tmp.col])
}
tmp.samp="CEZ-083"
mut_type_tally_DF5123out %>%
	filter(tier2=="SNP", sample1==tmp.samp) %>%
	group_by(sample1, sample2, log10_mSv_y_str1, log10_mSv_y_str2, cez_not_cez_str1, cez_not_cez_str2) %>%
	summarise(sumcount = sum(count),
				sumtotal = sum(total)) %>%
	mutate(rel_MR = (sumcount+50)/(50+sumtotal-sumcount)) %>%
	group_by(sample1) %>%
	mutate(avg_dist = mean(sumtotal),
			avg_rel_MR = mean(rel_MR)) %>%
	ggplot()+
	geom_hline(yintercept = 1)+
	geom_point(aes(x=sumtotal, y=rel_MR, color=log10_mSv_y_str1, group=sample1), size=5)+
	geom_point(aes(x=sumtotal, y=rel_MR,
					color=log10_mSv_y_str2, group=sample1))+
	rad_color_scale+
	coord_trans(y = "log10", x = "log10")+
	theme_sophie+
	ggtitle(tmp.samp)+
	labs(x="Number of mutations between pair",
			y="Ratio of unique mutations from\n each member of pair")+
	theme(legend.position = "none")
ggsave("plots_out/04B_CEZ-083_example.pdf",
		device = "pdf", width = 5, height = 5)

#~#~#
# Figure 4C
#~#~#

library(tidyverse)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

strain_order = unique(c(rad_table$sample[which(is.na(rad_table$log10_mSv_per_yr))],
                        rad_table %>% arrange(log10_mSv_per_yr) %>% select(sample) %>% unlist))
    
mut_type_tally_DF5123out %>%
    filter(tier2=="SNP") %>%
    group_by(sample1, sample2, log10_mSv_y_str1, log10_mSv_y_str2, cez_not_cez_str1, cez_not_cez_str2) %>%
    summarise(sumcount = sum(count),
              sumtotal = sum(total)) %>% 
    mutate(rel_MR = (sumcount+50)/(50+sumtotal-sumcount)) %>%
    group_by(sample1) %>%
    # mutate(avg_dist = mean(sumtotal),
    #        avg_rel_MR = mean(rel_MR)) %>%
    mutate(avg_dist = median(sumtotal),
           avg_rel_MR = median(rel_MR)) %>%
    left_join(rad_table %>% select(sample, twoline), by=c("sample1" = "sample")) %>% 
    mutate(sample1 = factor(sample1, levels = strain_order)) %>%

    ggplot()+
    geom_hline(yintercept = 1)+
    geom_point(aes(x=sample1, y=rel_MR, color=log10_mSv_y_str1, group=sample1), size=2)+
    geom_point(aes(x=sample1, y=rel_MR, color=log10_mSv_y_str2, group=sample1), size=.3)+
    rad_color_scale+
    geom_point(aes(x=sample1, y=avg_rel_MR, fill=log10_mSv_y_str1), shape=21, color="black", size=8, alpha=.05)+
    geom_text(aes(x=sample1, y=avg_rel_MR, label=twoline), 
              lineheight=.75, hjust=.5, vjust=.5, size=2, color="black")+
    rad_fill_scale+
    coord_trans(y = "log10")+  #, x = "log10"
    theme_sophie+
    labs(x="Numerator strain", 
         y="Ratio of unique mutations from\n each member of pair")+
    theme(legend.position = "none")

ggsave(paste0("plots_out/04C_avg_rel_MA_medians.pdf"), 
       device = "pdf", width = 5, height = 5)

#~#~#
# Supplementary Figure 10 - Relative mutation acquisition by mutation type (plus Mantel tests)
#~#~#

library(tidyverse)
library(adegenet)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# Mantel test on all mutations together
mut_acq_table = mut_type_tally_DF5123out %>%
    filter(tier2=="SNP") %>%
    group_by(sample1, sample2) %>%
    summarise(sumcount = sum(count),
              sumtotal = sum(total)) %>% 
    mutate(rel_MR = (sumcount+50)/(50+sumtotal-sumcount)) 
rad_dist = data.frame()
for(i in unique(diff_counts_anyDP05$samp1)){
    for(j in unique(diff_counts_anyDP05$samp2)){
        print(paste(i,j))
        rad_dist[nrow(rad_dist)+1, c("samp1", "samp2", "rad_dist")] = 
            c(i, j, (10^(rad_table[i,"log10_mSv_per_yr"]) / 10^(rad_table[j,"log10_mSv_per_yr"])))
    }
}
sub_samps = unique(mut_acq_table$sample1)
rad_dist.dist = rad_dist %>%
    filter(samp1%in%sub_samps, samp2%in%sub_samps) %>%
    arrange(samp1, samp2) %>%
    dcast(samp1~samp2, value.var="rad_dist")  %>%
    column_to_rownames(var="samp1") %>%
    dist()
mutacq_dist = mut_acq_table %>% 
    filter(sample1%in%sub_samps, sample2%in%sub_samps) %>%
    arrange(sample1, sample2) %>%
    select(sample1, sample2, rel_MR) %>%
    dcast(sample1~sample2, value.var="rel_MR") %>%
    column_to_rownames(var="sample1") %>% 
    dist()
mantel.randtest(m1 = rad_dist.dist, m2 = mutacq_dist, nrepet = 999)
# Quick plot of that relationship
full_join(rad_dist, mut_acq_table %>% rename("samp1" = "sample1", "samp2"="sample2")) %>% 
    filter(samp1 %in% CEZ_samples & samp2 %in% CEZ_samples) %>%
    ggplot()+
    geom_point(aes(x=as.numeric(rad_dist), y=rel_MR), 
               shape=21, color="black", fill="#a80000")+
    labs(x="fold difference in radiation", y = "relative mutation acquisition")+
    ggtitle("cez samples")+
    theme_sophie
ggsave("plots_out/Supp10Z_mantel_for_rad_vs_relMA.pdf", device = "pdf", width = 4, height = 4)

#~# Separate analysis out by mutation type
sub_samps = unique(mut_type_tally_DF5123out$sample1)
rad_dist = data.frame()
for(i in unique(diff_counts_anyDP05$samp1)){
    for(j in unique(diff_counts_anyDP05$samp2)){
        print(paste(i,j))
        rad_dist[nrow(rad_dist)+1, c("samp1", "samp2", "rad_dist")] = 
            c(i, j, (10^(rad_table[i,"log10_mSv_per_yr"]) / 10^(rad_table[j,"log10_mSv_per_yr"])))
    }
}
# Mantel test on mutations split by mut type
rad_dist.dist = rad_dist %>%
    filter(samp1%in%sub_samps, samp2%in%sub_samps) %>%
    arrange(samp1, samp2) %>%
    dcast(samp1~samp2, value.var="rad_dist")  %>%
    column_to_rownames(var="samp1") %>%
    dist()
mantel.type.table = data.frame(type=character(0),pvalue=numeric(0))
for(tmp.type in unique(mut_type_tally_DF5123out$tier2)){
    print(tmp.type)
    mut_acq_table = mut_type_tally_DF5123out %>%
        filter(tier2==tmp.type) %>%
        group_by(sample1, sample2) %>%
        summarise(sumcount = sum(count),
                  sumtotal = sum(total)) %>% 
        mutate(rel_MR = (sumcount+50)/(50+sumtotal-sumcount))
    mutacq_dist = mut_acq_table %>% 
        filter(sample1%in%sub_samps, sample2%in%sub_samps) %>%
        arrange(sample1, sample2) %>%
        select(sample1, sample2, rel_MR) %>%
        dcast(sample1~sample2, value.var="rel_MR") %>%
        column_to_rownames(var="sample1") %>% 
        dist()
    tmp.mantel = mantel.randtest(m1 = rad_dist.dist, m2 = mutacq_dist, nrepet = 999)
    mantel.type.table[tmp.type,c("type", "pvalue")] = c(tmp.type, tmp.mantel$pvalue)
}
for(tmp.type in unique(mut_type_tally_DF5123out$tier1)){
    print(tmp.type)
    mut_acq_table = mut_type_tally_DF5123out %>%
        filter(tier1==tmp.type) %>%
        group_by(sample1, sample2) %>%
        summarise(sumcount = sum(count),
                  sumtotal = sum(total)) %>% 
        mutate(rel_MR = (sumcount+50)/(50+sumtotal-sumcount))
    mutacq_dist = mut_acq_table %>% 
        filter(sample1%in%sub_samps, sample2%in%sub_samps) %>%
        arrange(sample1, sample2) %>%
        select(sample1, sample2, rel_MR) %>%
        dcast(sample1~sample2, value.var="rel_MR") %>%
        column_to_rownames(var="sample1") %>% 
        dist()
    tmp.mantel = mantel.randtest(m1 = rad_dist.dist, m2 = mutacq_dist, nrepet = 999)
    mantel.type.table[tmp.type,c("type", "pvalue")] = c(tmp.type, tmp.mantel$pvalue)
}
for(tmp.type in unique(mut_type_tally_DF5123out$tier0)){
    print(tmp.type)
    mut_acq_table = mut_type_tally_DF5123out %>%
        filter(tier0==tmp.type) %>%
        group_by(sample1, sample2) %>%
        summarise(sumcount = sum(count),
                  sumtotal = sum(total)) %>% 
        mutate(rel_MR = (sumcount+50)/(50+sumtotal-sumcount))
    mutacq_dist = mut_acq_table %>% 
        filter(sample1%in%sub_samps, sample2%in%sub_samps) %>%
        arrange(sample1, sample2) %>%
        select(sample1, sample2, rel_MR) %>%
        dcast(sample1~sample2, value.var="rel_MR") %>%
        column_to_rownames(var="sample1") %>% 
        dist()
    tmp.mantel = mantel.randtest(m1 = rad_dist.dist, m2 = mutacq_dist, nrepet = 999)
    mantel.type.table[tmp.type,c("type", "pvalue")] = c(tmp.type, tmp.mantel$pvalue)
}
write.csv2(mantel.type.table, "plots_out/Supp10Z_mantel_for_mut_types.csv")

# Plot of 10 mutation types
strain_order = unique(c(rad_table$sample[which(is.na(rad_table$log10_mSv_per_yr))],
                        rad_table %>% arrange(log10_mSv_per_yr) %>% select(sample) %>% unlist))
mut_order = unique(c("Ins_1-5", "Ins_6+", "Del_1-5", "Del_6+", unique(mut_type_tally_DF5123out$tier1)))
mut_type_tally_DF5123out %>%
    group_by(sample1, sample2, tier1, log10_mSv_y_str1, log10_mSv_y_str2, cez_not_cez_str1, cez_not_cez_str2) %>%
    summarise(sumcount = sum(count),
              sumtotal = sum(total)) %>% 
    mutate(rel_MR = (sumcount+50)/(50+sumtotal-sumcount)) %>%
    group_by(sample1, tier1) %>%
    mutate(avg_dist = median(sumtotal),
           avg_rel_MR = median(rel_MR)) %>%
    left_join(rad_table %>% select(sample, twoline), by=c("sample1" = "sample")) %>% 
    mutate(sample1 = factor(sample1, levels = strain_order)) %>%
    mutate(tier1 = factor(tier1, levels = mut_order)) %>%
    
    ggplot()+
    geom_hline(yintercept = 1)+
    geom_point(aes(x=sample1, y=rel_MR, color=log10_mSv_y_str1, group=sample1), size=2)+
    geom_point(aes(x=sample1, y=rel_MR, color=log10_mSv_y_str2, group=sample1), size=.3)+
    rad_color_scale+
    geom_point(aes(x=sample1, y=avg_rel_MR, fill=log10_mSv_y_str1), shape=21, color="black", size=8, alpha=.05)+
    geom_text(aes(x=sample1, y=avg_rel_MR, label=twoline), 
              lineheight=.75, hjust=.5, vjust=.5, size=2, color="black")+
    rad_fill_scale+
    coord_trans(y = "log10")+  #, x = "log10"
    facet_wrap(~tier1) +
    theme_sophie+
    labs(x="Numerator strain", 
         y="Ratio of unique mutations from\n each member of pair")+
    theme(legend.position = "none")
ggsave(paste0("plots_out/Supp10_median_by_mut_type.pdf"),
       device = "pdf", width = 16, height = 10)

# glm by mut type
tmp.data.for.stat = mut_type_tally_DF5123out %>%
    group_by(sample1, sample2,tier1, log10_mSv_y_str1, log10_mSv_y_str2, cez_not_cez_str1, cez_not_cez_str2) %>%
    summarise(sumcount = sum(count),
              sumtotal = sum(total)) %>%
    mutate(rel_MR = (sumcount+50)/(50+sumtotal-sumcount)) %>%
    ungroup %>%
    group_by(sample1, tier1, log10_mSv_y_str1, cez_not_cez_str1) %>%
    summarise(avg_dist = mean(sumtotal),
              avg_rel_MR = mean(rel_MR))
glm_for_10_muts = data.frame()
for(tmp.tier in unique(tmp.data.for.stat$tier1)){
    # for(tmp.tier in unique(tmp.data.for.stat$NC_S_NS)){
    #for(tmp.tier in unique(tmp.data.for.stat$tier2)){
    print(tmp.tier)
    temp.stat = glm(formula=(avg_rel_MR~log10_mSv_y_str1),  
                    # temp.stat = glm(formula=(avg_rel_MR~Bq),  
                    family=gaussian, 
                    #  data=tmp.data.for.stat %>% filter(tier2=="SNP", Bq>0))
                    data=tmp.data.for.stat %>% filter(tier1==tmp.tier))
    # data=tmp.data.for.stat %>% filter(NC_S_NS==tmp.tier))
    
    
    #summary(temp.stat)
    tmp.p = summary(temp.stat)$coefficients["log10_mSv_y_str1","Pr(>|t|)"]
    print(tmp.p)
    glm_for_10_muts[nrow(glm_for_10_muts)+1,c("tier1","p-val")] = c(tmp.tier, tmp.p)
}
write.table(glm_for_10_muts, file = "plots_out/Supp10Z_glm_for_mut_types.csv", 
            sep = ",", row.names = F)

# #~#~#
#  Power analysis
# #~#~#

# Simulations of branchlengths under elevated recent mutation rates, conditional on the genealogy of the Oscheius samples.

library(SNPRelate)
library(gdsfmt)
library(ape)
library(phangorn)
require(reshape2)

gds.folder ="supplementary_data/all_CEW1_no3_minis/"
gds.set = list.files(gds.folder, pattern = "*gds$")

### SPECIFY PARAMETERS
n.trees <- length(gds.set); #how many sections of genome are modeled as having their own genealogies. Should be 1400?
# multiplier <- 2 # a strain collected at 1Sv/yr exposure has extra mutations as a function of the mutation rate increasing by this factor since 1986
n.sims <- 1000 # how many simulations should be performed for a given set of parameters
baseline.mutation <- 2.5e-9
years.exposure <- 33 # 2019-1986
gens.per.year <- 91.25 # 365 days/year divided by 4 days/generation
bases.per.tree <- 40000000/n.trees # how many basepairs in each of the 1400 trees? Here approximated as 40M/number of trees

### Get radiation levels
source("supplementary_data/00_chornobyl_plot_shortcuts.R")
radlevels = rad_table
radlevels[c(1:15, 18,21),1:2] -> radlevels
radlevels[16:17,2] <- 0
radlevels[,3] <- 10^radlevels[,2]
radlevels[,3] <- radlevels[,3] -10 
radlevels[which(radlevels[,3] < 0),3] <- 0
radlevels[,4] <- radlevels[,3] / 990 
names(radlevels)[4] <- "radscale"

# Run simulation for a few multipliers
sims_by_multiplier = list()
for(multiplier in c(2,5,10,100)){
    # How many extra mutations do we expect each strain to carry? 
    expected.extra.mutations <- years.exposure * bases.per.tree * gens.per.year * 
        baseline.mutation * radlevels$radscale * (multiplier-1) # (mutations per tree)
    
    ### Prepare some useful functions
    # function for getting tree tip labels
    getedge <- function(x){which.edge(umt, x)}
    # pull out sections of lists of lists (straight from stack overflow)
    listslice <- function(lst, n){ sapply(lst, `[`, n) }
    
    ### Now loop through all the trees and generate n.sims simulations for each 
    # At the end, list.of.result.lists will be a list with one entry for each of the n.trees genealogies
    # Each of these entries will be a list with n.sims entries, each of which is a matrix of simulated pairwise SNP differences
    
    ### Initialize a list to hold results
    list.of.result.lists <- vector("list", n.trees) 
    ### Read in the trees as a gds object 
    set.seed(1)
    for(i in c(1:n.trees)){     # takes 5-20 minutes
        tree <- snpgdsOpen(paste0(gds.folder, gds.set[i]), allow.duplicate = T)
        
        dm <- snpgdsHCluster(snpgdsIBS(tree,autosome.only=FALSE, verbose = F)) 
        # root the tree and convert to ultrametric
        rt <- root(nj(dm$dist), outgroup = "DF5123", resolve.root=T) 
        umt <- nnls.tree(cophenetic(rt), rt, root =T)
        umt$edge.length[which(umt$edge.length < 0)]  <-0 # get rid of negative branch lengths
        simtree <- umt
        # pull out the positions of the tip branches
        strain.ids <- unlist(lapply(umt$tip.label, getedge))
        
        # Here comes the actual simulation part 
        # Each strain receives a Poisson-distributed number of mutations according to its expected number  
        extra.mutations <- matrix(nrow = 17, ncol = n.sims, rpois(17*n.sims, expected.extra.mutations))
        # Scale all the branchlengths to be 1000 SNPs after addition of the extra mutations
        null.edges <- apply(matrix(nrow = n.sims, length(simtree$edge.length), 
                                   data = 1000 - apply(extra.mutations, 2, sum)), 1, function(x){x * simtree$edge.length/sum(simtree$edge.length)})
        alt.edges <- null.edges
        alt.edges[strain.ids,] <- null.edges[strain.ids,] + extra.mutations
        # Now that the branchlengths are adjusted to represent extra mutations, drop 1000 mutations on the trees 
        sim.edge.lengths <- apply(alt.edges, 2, function(x){rmultinom(1,1000,x)})
        # Store the resulting pairwise distances as n.sims matrices
        list.of.result.lists[[i]] <- lapply(split(t(sim.edge.lengths), seq(NCOL(sim.edge.lengths))), function(x){simtree$edge.length <- x; cophenetic.phylo(simtree)})
    }
    
    ### Combine data across the simulated genealogical trees
    # The final result of a simulation is the pairwise distance between two strains summed across all 1400 trees. 
    # To get that we will add the mutations across the 1400 trees and then calculate pairwise branchlengths,
    # with DF5123 as the outgroup in analyses of each pair of strains. 
    
    # Initialize a list to hold results
    sim.branchlengths <- vector("list", n.sims)
    # Step through the n.sims sims. (There's probably a way to do this with map or mapply)
    for(j in 1:n.sims){    # took 1 minute
        full.distance <- Reduce("+", listslice(list.of.result.lists, j))
        sim.branchlengths[[j]] <- matrix(nrow = 16, ncol = 16, 0, dimnames = list(dimnames(full.distance)[[1]][1:16],dimnames(full.distance)[[2]][1:16] ))
        for(m in 1:16){
            for(n in 1:16){
                sim.branchlengths[[j]][m,n] <- full.distance[m,n]/2 + full.distance[m,17]/2 - full.distance[n,17]/2
            }}}
    # This gives a branchlength matrix for each simulation, where entry [m,n] gives the length of the branch leading to strain m 
    # based on the pairwise distances between m, n, and DF5123. 
    sim.pvals = lapply(sim.branchlengths, function(x){ # 5 minutes
        # combine counts and totals into one table
        tmp.counts.df = x %>%
            melt(value.name="count") %>%
            filter(Var1!=Var2)
        tmp.totals.df = (x + t(x)) %>%
            melt(value.name="total") %>%
            filter(Var1!=Var2)
        tmp.df = tmp.counts.df %>%
            left_join(tmp.totals.df) %>%
            left_join( rad_table %>% select(sample, log10_mSv_per_yr) %>%
                           rename("Var1"="sample"))
        # calculate avg_rel_MR for each strain
        tmp.df = tmp.df %>%
            mutate(relMR = (count+50)/(total-count+50)) %>%
            group_by(Var1, log10_mSv_per_yr) %>%
            summarize(avg_rel_MRs = mean(relMR))
        # And glm
        glm.stat = glm(formula=(avg_rel_MRs~log10_mSv_per_yr),
                       family=gaussian,
                       data=tmp.df)
        glm.p = summary(glm.stat)$coefficients["log10_mSv_per_yr","Pr(>|t|)"]
        return(glm.p)
    }) %>% unlist()
    Sys.sleep(5)
    sims_by_multiplier[[paste0("mult.",multiplier)]] = sim.pvals
}

for(tmp.name in names(sims_by_multiplier)){
    print(tmp.name)
    tmp.perc = round(100*length(which(sims_by_multiplier[[tmp.name]]<0.05))/length(sims_by_multiplier[[tmp.name]]), digits = 1)
    print(paste0(as.character(tmp.perc), "%"))
}
# [1] "mult.2"
# [1] "22.8%"
# [1] "mult.5"
# [1] "64.7%"
# [1] "mult.10"
# [1] "98.4%"
# [1] "mult.100"
# [1] "100%"

#~#~#
# Supplementary Figure 11 - JU75 (O. sp 3) as outgroup
#~#~#

load("supplementary_data/12_mut_type_strains_SNS_anc_rad_JU75_check.RData")
mut_type_strains_SNS_anc_rad_JU75_check %>%
    filter(tier2=="SNP") %>%
    group_by(sample1, sample2, log10_mSv_y_str1, log10_mSv_y_str2, cez_not_cez_str1, cez_not_cez_str2) %>%
    summarise(sumcount = sum(count),
              sumtotal = sum(total)) %>% 
    mutate(rel_MR = (sumcount+50)/(50+sumtotal-sumcount)) %>%
    group_by(sample1) %>%
    mutate(avg_rel_MR = mean(rel_MR)) %>%
    left_join(rad_table %>% select(sample, twoline), by=c("sample1" = "sample")) %>%
    mutate(sample1 = factor(sample1, levels = unique(c("DF5062", "DF5103", "DF5110", "DF5111", "DF5123", rad_table %>% arrange(log10_mSv_per_yr) %>% select(sample) %>% unlist())))) %>%
    ggplot()+
    geom_hline(yintercept = 1)+
    geom_point(aes(x=sample1, y=rel_MR, color=log10_mSv_y_str1, group=sample1), size=2)+
    geom_point(aes(x=sample1, y=rel_MR, color=log10_mSv_y_str2, group=sample1), size=.3)+
    geom_point(aes(x=sample1, y=avg_rel_MR), color="black", size=8.4)+
    geom_point(aes(x=sample1, y=avg_rel_MR, color=log10_mSv_y_str1), size=8)+
    # geom_text(aes(x=avg_dist, y=avg_rel_MR, color=log10_mSv_y_str1, label=sample1),
    #           hjust=-4, vjust=.5, size=2, color="black", angle=90)+
    geom_text(aes(x=sample1, y=avg_rel_MR, color=log10_mSv_y_str1, 
                  #label=rad_table$twoline[sample1]
                  label=twoline,
                  # label=sample1
    ), 
    lineheight=.75, hjust=.5, vjust=.5, size=2, color="black")+
    #  geom_smooth(method='lm',formula=y~x, se = F)+
    rad_color_scale+
    coord_trans(y = "log10")+
    theme_sophie+
    labs(x="Numerator strain", 
         y="Ratio of unique mutations from\n each member of pair")+
    ggtitle("all mutations")+
    theme(legend.position = "none")
ggsave("plots_out/Supp11A_nuclear_by_strain.pdf", 
       device = "pdf", width = 5, height = 5)

# and the mitochondria 
load("supplementary_data/12_mito_muts.RData")
mito_muts$sample1 = as.character(mito_muts$sample1)
mito_muts$sample2 = as.character(mito_muts$sample2)
mito_muts$log10_mSv_y_str2 = rad_table[mito_muts$sample2,"log10_mSv_per_yr"]
mito_muts = mito_muts %>% 
    left_join(rad_table %>% select(sample, twoline) %>% rename("sample1"="sample"))

mito_muts %>% 
    filter(tier2=="SNP") %>% 
    group_by(sample1, sample2, log10_mSv_y_str1, log10_mSv_y_str2, cez_not_cez_str1, cez_not_cez_str2) %>%
    summarise(sumcount = sum(count),
              sumtotal = sum(total)) %>% 
    mutate(rel_MR = (sumcount+1)/(1+sumtotal-sumcount)) %>% 
    group_by(sample1) %>%
    mutate(avg_rel_MR = mean(rel_MR)) %>% 
    left_join(rad_table %>% select(sample, twoline), by=c("sample1" = "sample")) %>%
    mutate(sample1 = factor(sample1, levels = unique(c("DF5062", "DF5103", "DF5110", "DF5111", "DF5123", rad_table %>% arrange(log10_mSv_per_yr) %>% select(sample) %>% unlist())))) %>%
    ggplot()+
    geom_hline(yintercept = 1)+
    geom_point(aes(x=sample1, y=rel_MR, color=log10_mSv_y_str1, group=sample1), size=2)+
    geom_point(aes(x=sample1, y=rel_MR, color=log10_mSv_y_str2, group=sample1), size=.3)+
    geom_point(aes(x=sample1, y=avg_rel_MR), color="black", size=8.4)+
    geom_point(aes(x=sample1, y=avg_rel_MR, color=log10_mSv_y_str1), size=8)+
    # geom_text(aes(x=avg_dist, y=avg_rel_MR, color=log10_mSv_y_str1, label=sample1),
    #           hjust=-4, vjust=.5, size=2, color="black", angle=90)+
    geom_text(aes(x=sample1, y=avg_rel_MR, color=log10_mSv_y_str1, 
                  label=twoline), lineheight=.75,
              hjust=.5, vjust=.5, size=2, color="black")+
    #  geom_smooth(method='lm',formula=y~x, se = F)+
    rad_color_scale+
    coord_trans(y = "log10")+
    theme_sophie+    
    labs(x="Numerator strain", 
         y="Ratio of unique mutations from\n each member of pair")+
    ggtitle("mitochondrial mutations")+
    theme(legend.position = "none")
ggsave("plots_out/Supp11B_mito_by_strain.pdf", 
       device = "pdf", width = 5, height = 5)


#~#~#
# Microscatellites analysis
#~#~#

require(tidyverse)
require(reshape2)
source("~/Dropbox/postdoc/012_R/07_chernobyl_plot_shortcuts.R")

samps = c("CEZ-010", "CEZ-047", "CEZ-061", "CEZ-077", "CEZ-083", "CEZ-107", 
          "CEZ-119", "CEZ-138", "CEZ-144", "CEZ-156", "CEZ-166", "CEZ-224", 
          "CEZ-253", "CEZ-263", "CEZ-274", "DF5103")

# load vcfs for tandem repeats and correct the format
hipstr_df5123 = read.table("supplementary_data/13_hpster_DF5123.str.vcf", 
                           col.names = c("CHROM", "POS", "ID", "REF", "ALT", "dot1", 
                                         "dot2", "INFO", "FORMAT", samps ))

dotnames16 = intersect(colnames(hipstr_df5123), rad_table$dot_names)
dashnames16 = unlist(lapply(dotnames16, function(x){rad_table$sample[which(rad_table$dot_names==x)]})) 
for(tmp.samp in dotnames16){
    print(tmp.samp)
    hipstr_df5123 = hipstr_df5123 %>%
        separate(col = as.symbol(tmp.samp), into = c(paste0(tmp.samp, "_GT"), paste0(tmp.samp, "_GB"), paste0(tmp.samp, "_drop1"),
                                                     paste0(tmp.samp, "_drop2"), paste0(tmp.samp, "_DP"), paste0(tmp.samp, "_drop3"),
                                                     paste0(tmp.samp, "_DSTUTTER"), paste0(tmp.samp, "_DFLANKINDEL"),
                                                     paste0(tmp.samp, "_PDP"), paste0(tmp.samp, "_drop4"), paste0(tmp.samp, "_GLDIFF")), 
                 sep = ":", extra = "drop") %>%
        select(-c(paste0(tmp.samp, "_drop1"), paste0(tmp.samp, "_drop2"), paste0(tmp.samp, "_drop3"), paste0(tmp.samp, "_drop4")))
}

# Correct non-convincing hets
for(samp1 in dotnames16){
    print(samp1)
    tmp.mini = data.frame(GT = hipstr_df5123[,paste0(samp1,"_GT")], PDP = hipstr_df5123[,paste0(samp1,"_PDP")])
    tmp.users = which(!is.na(tmp.mini$PDP))
    tmp.mini = tmp.mini %>% 
        slice(tmp.users) %>%
        separate(col = GT, into = c("GT_1", "GT_2"), sep = "\\|") %>%
        separate(col = PDP, into = c("PDP_1", "PDP_2"), sep = "\\|") 
    for(tmp.col in c("GT_1", "GT_2", "PDP_1", "PDP_2")){
        tmp.mini[,tmp.col] = as.numeric(tmp.mini[,tmp.col])}
    # get rid of hets we don't believe
    tmp.mini = tmp.mini %>%
        mutate(GT_1 = case_when( ((PDP_1/(PDP_1+PDP_2))<.4) ~ GT_2, .default = GT_1)) %>%
        mutate(GT_2 = case_when( ((PDP_2/(PDP_1+PDP_2))<.4) ~ GT_1, .default = GT_2)) 
    # now put them back
    hipstr_df5123[tmp.users,paste0(samp1,"_GT")] = paste(tmp.mini$GT_1, tmp.mini$GT_2, sep = "|")
}

# Now get some distances
microsat.table_df5123 = data.frame(samp1 = character(0), samp2 = character(0), numb_sites = integer(0), numb_diffs = integer(0))
for(samp1 in dotnames16){
    print(samp1)
    for(samp2 in dotnames16){
        if(samp1==samp2){next}
        mini.vcf = hipstr_df5123 %>% select(c("CHROM", "POS", "ID", "REF", "ALT", 
                                              paste0(samp1, "_GT"), paste0(samp1, "_DP"), paste0(samp1, "_PDP"), 
                                              paste0(samp2, "_GT"), paste0(samp2, "_DP"), paste0(samp2, "_PDP")))
        colnames(mini.vcf) = c("CHROM", "POS", "ID", "REF", "ALT", "samp1_GT", "samp1_DP", "samp1_PDP", "samp2_GT", "samp2_DP", "samp2_PDP")
        
        mini.vcf = mini.vcf %>%
            filter(as.numeric(samp1_DP)>4, as.numeric(samp2_DP)>4)
        
        # For each pair, how many sites to compare, and how many differences
        tmp.diffs = length(which(mini.vcf$samp1_GT!=mini.vcf$samp2_GT))
        
        newrown = nrow(microsat.table_df5123)+1
        microsat.table_df5123[newrown, c("samp1", "samp2")] = c(samp1, samp2) 
        microsat.table_df5123[newrown, c("numb_sites", "numb_diffs")] = c(nrow(mini.vcf),tmp.diffs)
    }
}
microsat.table_df5123$frac_var = microsat.table_df5123$numb_diffs / microsat.table_df5123$numb_sites
microsat.table_df5123$samp1 = unlist(lapply(microsat.table_df5123$samp1, function(x){rad_table$sample[which(rad_table$dot_names==x)]}))
microsat.table_df5123$samp2 = unlist(lapply(microsat.table_df5123$samp2, function(x){rad_table$sample[which(rad_table$dot_names==x)]}))

# How many loci have reads for all 16 samples?
gt_table_df5123 = hipstr_df5123[,c(10,17,24,31,38,45,52,59,66,73,80,87,94,101,108,115)] 
gt_table_df5123[gt_table_df5123 == "."] <- NA
gt_table_df5123.levels = unique(unlist(gt_table_df5123))
for(tmp.col in colnames(gt_table_df5123)){
    gt_table_df5123[,tmp.col] = factor(gt_table_df5123[,tmp.col], levels = gt_table_df5123.levels)}
# Reduce each of these to only lines where there are reads for all 16 samples
gt_table_df5123 = na.omit(gt_table_df5123)   # 469 -> 147
# That is an extremely  small number

# MAKE BLUE MATRIX
actual_clust = microsat.table_df5123 %>% 
    dcast(samp1~samp2, value.var = "frac_var") %>%
    column_to_rownames(var = "samp1") %>% 
    as.dist() %>% hclust()
snp.order = c("CEZ-061", "CEZ-107", "CEZ-224", "CEZ-138", "CEZ-144", "CEZ-083", "DF5103", "CEZ-077", 
              "CEZ-166", "CEZ-047", "CEZ-156", "CEZ-263", "CEZ-010", "CEZ-274", "CEZ-119", "CEZ-253")
microsat.table_df5123 %>%
    ggplot(aes(x=factor(samp1, levels =   snp.order),# actual_clust$labels[actual_clust$order]),
               y=factor(samp2, levels =   snp.order), # actual_clust$labels[actual_clust$order]),
               fill=frac_var))+
    geom_tile(color="#FFFFFF", linewidth=.005)+
    coord_fixed()+
    theme_sophie+
    labs(x=NULL, y=NULL)+
    theme(legend.title = element_blank())+
    scale_fill_gradient(low="#FFFFFF", high="#08306b") +
    theme(axis.text.x = element_text(face="bold", colour = rad_colors[snp.order]))+
    theme(axis.text.y = element_text(face="bold", colour = rad_colors[snp.order]))
ggsave("plots_out/SuppXA_microsat_distances.pdf", 
       device = "pdf", width = 8, height = 6)
# all looks the same broadly as SNP data

# PLOT mut distance BY RAD DISTANCE?
microsat.table_df5123$log10raddiff = unlist(apply(microsat.table_df5123, 1, function(x){
    rad1 = rad_table[x["samp1"], "log10_mSv_per_yr"]
    rad2 = rad_table[x["samp2"], "log10_mSv_per_yr"]
    return(log10(rad1/rad2))
}))
microsat.table_df5123 %>%
    filter(log10raddiff>0) %>%
    ggplot(aes(x=log10raddiff, y=frac_var)) +
    geom_point()+
    labs(x="log10(rad1/rad2)", y="variant/total microsatellites")+
    theme_sophie
ggsave("plots_out/SuppXB_rad_vs_genet_df5123.pdf",
       device = "pdf", width = 6, height = 6)

# OR THE ORANGE DOTS (relative mutation acquisition), BY STRAIN
for(samp1 in dotnames16){
    print(samp1)
    for(samp2 in dotnames16){
        if(samp1==samp2){next}
        mini.vcf = hipstr_df5123 %>% select(c("CHROM", "POS", "ID", "REF", "ALT",
                                              paste0(samp1, "_GT"), paste0(samp1, "_DP"), paste0(samp1, "_PDP"),
                                              paste0(samp2, "_GT"), paste0(samp2, "_DP"), paste0(samp2, "_PDP")))
        colnames(mini.vcf) = c("CHROM", "POS", "ID", "REF", "ALT", "samp1_GT", "samp1_DP", "samp1_PDP", "samp2_GT", "samp2_DP", "samp2_PDP")
        
        mini.vcf = mini.vcf %>%
            filter(as.numeric(samp1_DP)>4, as.numeric(samp2_DP)>4) %>%
            filter(samp1_GT=="0|0" | samp2_GT == "0|0") %>%
            filter(samp1_GT!=samp2_GT)
        
        dashsamp1 = rad_table$sample[which(rad_table$dot_names==samp1)]
        dashsamp2 = rad_table$sample[which(rad_table$dot_names==samp2)]
        tmp.rown = which(microsat.table_df5123$samp1==dashsamp1 & microsat.table_df5123$samp2==dashsamp2)
        microsat.table_df5123[tmp.rown,c("der_micros", "anc_micros")] = c(1+length(which(mini.vcf$samp2_GT=="0|0")), 
                                                                          1+length(which(mini.vcf$samp1_GT=="0|0")))
    }
}
microsat.table_df5123$rel_MS_MR = microsat.table_df5123$der_micros / microsat.table_df5123$anc_micros
microsat.table_df5123 = microsat.table_df5123 %>%
    left_join(rad_table %>%
                  select(sample, log10_mSv_per_yr) %>%
                  rename("samp1"="sample", "log10_mSv_per_yr_1" = "log10_mSv_per_yr")) %>%
    left_join(rad_table %>% 
                  select(sample, log10_mSv_per_yr) %>%
                  rename("samp2"="sample", "log10_mSv_per_yr_2" = "log10_mSv_per_yr"))
microsat.table_df5123 %>%
    ggplot()+
    geom_point(aes(x=factor(samp1, levels = c("DF5103",rad_order)),
                   y=rel_MS_MR,
                   color=log10_mSv_per_yr_1))+
    geom_point(aes(x=factor(samp1, levels = c("DF5103",rad_order)),
                   y=rel_MS_MR,
                   color=log10_mSv_per_yr_2), size=.5, shape=18)+
    rad_color_scale+
    coord_trans(y="log10")+
    labs(x="strain", y="relative mutation aquisition - microsatellites")+
    theme_sophie
ggsave("plots_out/SuppXC_rel_MA.pdf", 
       device = "pdf", width = 7, height = 5)

#~#~#
# Supplementary figure 12 - cosmic signatures
#~#~#

require(tidyverse)
require(sigminer)
require(NMF)
require(reshape2)
require(adegenet)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# Read in mut_type_tally_DF5123out from "../supplementary_data/10_DF5123_mut_type_tally.csv"
# ? mut_type_tally_DF5123out = read.csv2("supplementary_data/11_DF5123_mut_type_tally.csv", sep = ",")

# Split into matrix
muttypes = unique(mut_type_tally_DF5123out$tier0[which(mut_type_tally_DF5123out$tier2=="SNP")])
muttypes = data.frame(row.names = muttypes, oldnames = muttypes, newnames=muttypes)
muttypes$newnames = unlist(lapply(muttypes$oldnames, function(x){
    tmp.parts = strsplit(x,"")[[1]]
    return(paste0(tmp.parts[1], "[", paste(tmp.parts[c(2,4,5)], collapse = ""), "]", tmp.parts[3]))}))
mut_matrix = data.frame()
for(tmp.samp1 in unique(mut_type_tally_DF5123out$sample1)){  # this takes 5 minutes
    for(tmp.samp2 in unique(mut_type_tally_DF5123out$sample2)){
        if(tmp.samp1==tmp.samp2){ next }
        tmp.rowname = paste0(tmp.samp1, "_to_", tmp.samp2)
        print(tmp.rowname)
        for(tmp.muttype in muttypes$oldnames){
            tmp.colname = muttypes[tmp.muttype, "newnames"]
            tmp.count = mut_type_tally_DF5123out %>% 
                filter(sample1==tmp.samp1, sample2==tmp.samp2, tier0==tmp.muttype) %>% 
                select(count) %>% sum()
            mut_matrix[tmp.rowname, tmp.colname] = tmp.count
}}}

# Top of figure (just three closely related pairs):
mini_matrix = mut_matrix[c("CEZ-047_to_CEZ-156", "CEZ-156_to_CEZ-047", "CEZ-119_to_CEZ-253", "CEZ-253_to_CEZ-119"),]
norm_mini = apply(mini_matrix, 1, function(x){x/sum(x)})
mini_sim = get_sig_similarity(Signature = norm_mini, set_order = T)
mini_sims_plottable = mini_sim$similarity %>%
    melt(varnames = c("pair", "signature"), value.name = "similarity") %>%
    separate(col = pair, into = c("samp1", "samp2"), sep = "_to_", remove = F) %>% 
    mutate(set = case_when(samp1 %in% c("CEZ-047", "CEZ-156") ~ "CEZ-047/CEZ-156",
                           samp1 %in% c("CEZ-119", "CEZ-253") ~ "CEZ-119/CEZ-253")) %>%
    mutate(hilo = case_when(samp1 %in% c("CEZ-047", "CEZ-119") ~ "high mut",
                           samp1 %in% c("CEZ-156", "CEZ-253") ~ "low mut"))
# make order for signatures
mini_sims_order = mini_sims_plottable %>% 
    group_by(signature) %>%
    summarise(avg_sim = mean(similarity)) %>% 
    arrange(desc(avg_sim)) %>% select(signature)
mini_sims_order = as.character(mini_sims_order$signature)
mini_sims_plottable$signature = factor(mini_sims_plottable$signature, levels = mini_sims_order)
mini_sims_plottable$set = factor(mini_sims_plottable$set, levels = c("CEZ-047/CEZ-156", "CEZ-119/CEZ-253"))
mini_sims_plottable %>%
    ggplot(aes(x=as.integer(signature)+.15*as.integer(set)-.225, y=similarity, color=hilo, shape=set)) +
    geom_vline(xintercept = c(1:72), color="#AAAAAA", alpha=.1)+
    geom_point() +
    scale_color_manual(values=c("high mut" = "#e08214", "low mut"  = "#8073ac"))+
    theme_sophie +
    scale_x_continuous(name="COSMIC signatures", 
                       breaks = c(1:72), labels=mini_sims_order)+
    theme(legend.position = "none")
ggsave("plots_out/Supp12A_sigs_for_three_pairs.pdf", 
       device = "pdf", width = 12, height = 4) 

# Bottom of figure (all pairwise comparisons)
norm_mut_matrix = apply(mut_matrix, 1, function(x){x/sum(x)})
sim = get_sig_similarity(Signature = norm_mut_matrix, set_order = T)
sims_plottable = sim$similarity %>%
    melt(varnames = c("pair", "signature"), value.name = "similarity") %>%
    separate(col = pair, into = c("samp1", "samp2"), sep = "_to_", remove = F) %>%
    left_join(rad_table%>% rename("samp1"="sample"))
# make order for signatures
sims_order = sims_plottable %>% 
    group_by(signature) %>%
    summarise(avg_sim = mean(similarity)) %>% 
    arrange(desc(avg_sim)) %>% select(signature)
sims_order = as.character(sims_order$signature)
sims_plottable$signature = factor(sims_plottable$signature, levels = sims_order)
sims_plottable %>%
    group_by(samp1,signature, log10_mSv_per_yr) %>%
    summarise(avg_similarity = mean(similarity)) %>%
    ggplot(aes(x=as.integer(signature), y=avg_similarity, color=log10_mSv_per_yr)) +
    geom_vline(xintercept = c(1:72),  alpha=.1)+
    rad_color_scale+
    geom_jitter(width = .15, size=.5) +
    theme_sophie +
    labs(x="COSMIC signatues")+
    theme(legend.position = "none")+
    scale_x_continuous(name="COSMIC signatures", breaks = c(1:72), labels=sims_order)
ggsave("plots_out/Supp12B_sigs_for_all_pairs.pdf", 
       device = "pdf", width = 12, height = 4) 

# Mantel test on each SBS to see if any correspond with average rel. mut aqcuisition
sims_plottable = sims_plottable %>% 
    left_join(rad_table %>% 
                  select(sample, log10_mSv_per_yr) %>% 
                  rename("samp1"="sample", "log10_mSv_per_yr_str1" = "log10_mSv_per_yr")) %>%
    left_join(rad_table %>% 
                  select(sample, log10_mSv_per_yr) %>% 
                  rename("samp2"="sample", "log10_mSv_per_yr_str2" = "log10_mSv_per_yr"))%>% 
    mutate(fold_change_rad = log10_mSv_per_yr_str1/log10_mSv_per_yr_str2)
# mantel on each SBS
sbs_mantel = data.frame()
# make the rad distance matrix
rad_fc_dist = sims_plottable %>%
    select(samp1, samp2, fold_change_rad) %>%
    distinct() %>%
    dcast(samp1~samp2) %>% 
    column_to_rownames(var="samp1") %>%
    dist()
for(tmp.sbs in unique(sims_plottable$signature)){
    # make the sbs distance matrix
    sbs_dist = sims_plottable %>%
        filter(signature==tmp.sbs) %>%
        select(samp1, samp2, similarity) %>%
        distinct() %>%
        dcast(samp1~samp2) %>%
        column_to_rownames(var="samp1") %>%
        dist()
    # test mantel and push to dataframe
    tmp.mantel = mantel.randtest(m1=rad_fc_dist, m2=sbs_dist, nrepet = 999)
    sbs_mantel[nrow(sbs_mantel)+1,c("signatue", "pval")] = 
        c(tmp.sbs, tmp.mantel$pvalue)
}
range(sbs_mantel$pval) # 1, all 1

#~#~#
# Figure 5 + Supplementary figure 13 + Supplementary figure 15
#~#~#

require(tidyverse)
require(RColorBrewer)
library(lubridate)
library(reshape2)
library(scales)
library(broom)  
library(cowplot)
library(forcats)
library(ggbeeswarm)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# useful function
tm <- function(tmp.time){ as.POSIXct(strptime(tmp.time, format="%y-%m-%d_%H-%M")) }

# Supp figures 12 and 14, and all of figure 5 are generated with this loop:
dir.create("plots_out/05_mutagen_tolerance")
dir.create("plots_out/Supp13_population_growth_traces")
dir.create("plots_out/Supp15_dose_response_elegans_tipulae")

finals.table = data.frame(file=c("supplementary_data/14_data_for_scanner_experiments/Cisp_20220118_DR.RData", 
                                "supplementary_data/14_data_for_scanner_experiments/Cisp_20220816_big21Ot.RData", 
                                 "supplementary_data/14_data_for_scanner_experiments/EMS_20220107_DR.RData", 
                                 "supplementary_data/14_data_for_scanner_experiments/EMS_20220628_big22Ot.RData",
                                 "supplementary_data/14_data_for_scanner_experiments/HU_20220316_DR.RData", 
                                 "supplementary_data/14_data_for_scanner_experiments/HU_20221019_big21Ot.RData"),
                          object=c("cisp_DR", "cisp_big21", "EMS_DR", "EMS_big21", "HU_DR", "HU_big21"),
                          title=c("Cisp_1DR", "Cisp_2big", "EMS_1DR", "EMS_2big", "HU_1DR", "HU_2big"),
                          doses=c(paste(c("0uM", "5uM", "20uM", "50uM", "100uM"), collapse = ","),
                                  paste(c("0uM",  "75uM"), collapse = ","),
                                  paste(c("0mM", "0.2mM","1mM","5mM","25mM","50mM","100mM","200mM"), collapse = ","),
                                  paste(c("0mM" , "10mM"), collapse = ","),
                                  paste(c("0mM","0.2mM","1mM","5mM","10mM","25mM","50mM","100mM"), collapse = ","),
                                  paste(c("0mM","10mM"), collapse = ",")))
pinkgreen8 = c("#4d9221", "#7fbc41", "#b8e186", "#e6f5d0", "#fde0ef", "#f1b6da", "#de77ae", "#c51b7d")
for(k in 1:nrow(finals.table)){
    print(finals.table$title[k])
    # Load all datasets
    load(finals.table$file[k])     
	temp.output <- get(finals.table$object[k])
    temp.doses = strsplit(finals.table$doses[k], ",")[[1]]
    # set how many/which colors
    custom_color = NULL
    color_steps=c(1,8,6,3,7,2,5,4)
    if(length(temp.doses)<9){ custom_color = scale_color_manual(values = pinkgreen8[sort(color_steps[1:length(temp.doses)])], limits=temp.doses)}
    if(length(temp.doses)==9){ custom_color = scale_color_manual(values = c(pinkgreen8, "#821252"), limits=temp.doses)}
    if(length(temp.doses)==13){ custom_color = scale_color_manual(values = c(pinkgreen8, "#960EBE", "#6600FF", "#5555FF", "#44AAFF", "#33FFFF"), limits=temp.doses)}
    
    # make plots of times to starvation (left part of figure 5, and supp 16 for dose responses)
    tmp.angle=0
    if(length(temp.doses)!=2){tmp.angle<-60}
    p = temp.output$summary.sd.trimmed.fit.stats %>%
        filter(starve==T, contamination==F) %>%
        ggplot(aes(x=factor(dose, levels = temp.doses),
                   y=sd_hours_to_starve, color=factor(dose, levels = temp.doses)))+
        geom_point(size=3)+
        custom_color +
        geom_point(shape=1, color="black", size=3)+
        labs(x="Dose", y="Hours to starvation", color="Dose")+
        scale_x_discrete(limits=temp.doses)+
        theme_sophie+
        geom_text(aes(label=scales::scientific(time_p_val, digits = 2), x=dose,
                      y = mean(c(rep(min(sd_hours_to_starve, na.rm = T),8),
                                 rep(max(sd_hours_to_starve, na.rm = T),1)))),
                  size=2, color="dark gray", angle=tmp.angle)
    if(str_sub(finals.table$title[k], -3, -1)=="big"){ p = p + facet_wrap(~strain, dir = "v", nrow = 3)
    } else {p = p + facet_wrap(~strain, dir = "v")}
    print(p)
    if(length(temp.doses)==2){
        ggsave(paste0("plots_out/05_mutagen_tolerance/left_", finals.table$title[k], "_Bdots.pdf"), device = "pdf", width = 9, height = 5)
    } else {
        ggsave(paste0("plots_out/Supp15_dose_response_elegans_tipulae/", finals.table$title[k], "_Bdots.pdf"), device = "pdf", width = 8, height = 6) }
    
    # if there are only two conditions, 
    #   make ratio and p value plot with rad color (right half of figure 5), and
    #   make growth curve plots (Supp 14)
    if(length(temp.doses)==2){
        
        # plot growth curves (supp 14)
        temp.output$summary.sd = temp.output$summary.sd %>% mutate(time_elapsed = difftime(tm(timepoint), tm(time_fed), units = "days"))
        temp.output$summary.sd.trimmed = temp.output$summary.sd.trimmed %>% mutate(time_elapsed = difftime(tm(timepoint), tm(time_fed), units = "days"))
        p = temp.output$summary.sd %>%
            ggplot(aes(x=time_elapsed, y=sd,
                       color=factor(dose, levels = temp.doses), group=sample)) +
            geom_line(alpha=.1) +
            geom_line(data = temp.output$summary.sd.trimmed %>% filter(contamination==F, starve==T),
                      aes(x = time_elapsed, y=sd,
                          color=factor(dose, levels = temp.doses), group=sample))+
            geom_point(data=temp.output$summary.sd.traits %>% filter(contamination==F, starve==T),
                       aes(x=sd_hours_to_starve/24,
                           y=max(temp.output$summary.sd$sd)-as.numeric(factor(dose, levels = temp.doses)),
                           color=factor(dose, levels = temp.doses)))+
            custom_color +
            facet_wrap(~strain, dir = "v", nrow = 3)+
            labs(x="Time (Days)", y="Worm density (St. dev. of pixel intensity)", color="Treatment") +
            theme_sophie
        print(p)
        ggsave(paste0("plots_out/Supp13_population_growth_traces/", finals.table$title[k], "_Atraces.pdf"), device = "pdf", width = 9, height = 5)

        # make ratio table (right half of figure 5)
        ratio_table = temp.output$summary.sd.trimmed.fit.stats %>%
            filter(dose %in% temp.doses) %>%
            group_by(strain, dose) %>%
            summarise(avg_time_to_starve = mean(sd_hours_to_starve)) %>%
            dcast(strain ~ dose)  %>% 
            full_join(temp.output$summary.sd.trimmed.fit.stats %>% 
                          filter(dose %in% temp.doses) %>%
                          filter(!is.na(time_p_val)) %>%
                          select(strain, time_p_val)%>%
                          distinct()) %>% 
            select(-dose) %>% 
           rename(ref = temp.doses[1], exp = temp.doses[2])
        require(ggrepel)
        p = ratio_table %>%
            left_join(rad_table, by = c("strain" = "sample")) %>%
            filter(!is.na(time_p_val)) %>%
            ggplot(aes(y = -log10(time_p_val), x = exp/ref, label = strain, color=log10_mSv_per_yr))+
            geom_point()+
            rad_color_scale+
            geom_text_repel(hjust=0.5, vjust=0, force = .5, nudge_y = .1, nudge_x = -.005, size=2)+
            expand_limits(x=1,y=0)+
            geom_vline(xintercept = 1)+
            ylab("-log10( P value of difference )")+
            xlab("Treatment / Ctrl (hours to starvation)")+
            theme_sophie+
            theme(legend.position = "none")
        print(p)
        ggsave(filename = paste0("plots_out/05_mutagen_tolerance/right_", finals.table$title[k], "_Dcolorratio.pdf"), device = "pdf", width = 4, height = 4)
    }
}

#~#~#
# Supplementary figure 14
#~#~#

library(tidyverse)
library(reshape2)
source("supplementary_data/00_chornobyl_plot_shortcuts.R")

# Load "clean_endpoints" from "../supplementary_data/14_all_endpoints.csv"
in_endpoints = read.csv2("supplementary_data/15_all_endpoints.csv", sep = ",")
clean_endpoints = in_endpoints %>% 
    filter(contamination==F, starve==T) %>% 
    group_by(strain, dosed, mutagen, time_t_stat, time_p_val) %>% 
    summarise(avg_time_to_starve = mean(as.numeric(sd_hours_to_starve))) %>% 
    dcast(formula = strain+mutagen~dosed, value.var = "avg_time_to_starve") %>% 
    rename("undosed"="FALSE", "dosed"="TRUE") %>% 
    group_by(strain, mutagen) %>%
    summarise(starve_ratio= dosed/undosed) %>%
    left_join(in_endpoints %>% select(strain, time_t_stat, time_p_val, mutagen) %>% 
                  filter(!is.na(time_t_stat)) %>% distinct()) %>%
    left_join(rad_table %>% select(sample, log10_mSv_per_yr) %>% rename("strain"="sample"))

strain_order = rad_table %>% arrange(log10_mSv_per_yr) %>% select(sample) %>% unlist()
clean_endpoints$strain = factor(clean_endpoints$strain, levels = strain_order)
clean_endpoints %>%
	ggplot(aes_string(x="strain", y="starve_ratio", color="log10_mSv_per_yr")) +   # OR time_p_val OR time_t_stat
	geom_point()+
	rad_color_scale+
	facet_wrap(~mutagen, ncol = 1, scales="free_y")+
	theme_sophie + 
	theme(legend.position = "none") +
	xlab("Strain")+ylab("Avg starve time (with mutagen / without mutagen)")
ggsave("plots_out/Supp14_strain_by_mutagen_tolerance.pdf", device = "pdf", width = 5, height = 5)

# glm p-values for each facet
stats_table_normal = data.frame()
clean_endpoints = clean_endpoints %>% 
    mutate(cez_non_cez = case_when(is.na(log10_mSv_per_yr) ~ "nonCEZ",
                                   !is.na(log10_mSv_per_yr) ~ "CEZ"))
for(tmp.mut in unique(clean_endpoints$mutagen)){
    print(tmp.mut)
    # cez_non_CEZ
    tmp.table = clean_endpoints %>% filter(mutagen==tmp.mut) %>% 
        select(strain, starve_ratio, cez_non_cez)
    tmp.fit = glm(formula = starve_ratio ~ cez_non_cez, data = tmp.table)
    stats_table_normal[nrow(stats_table_normal)+1, c("mutagen", "metric", "var", "pval")] = c(tmp.mut, "starve_ratio", "cez_non_cez", summary(tmp.fit)$coefficients[2,"Pr(>|t|)"])
    # rad
    tmp.table = clean_endpoints %>% filter(mutagen==tmp.mut, cez_non_cez=="CEZ") %>% select(strain, starve_ratio, log10_mSv_per_yr)
    tmp.fit = glm(formula = starve_ratio ~ log10_mSv_per_yr, data = tmp.table)
    stats_table_normal[nrow(stats_table_normal)+1, c("mutagen", "metric", "var", "pval")] = c(tmp.mut, "starve_ratio", "rad", summary(tmp.fit)$coefficients[2,"Pr(>|t|)"])
}
write.table(x = stats_table_normal, file = "plots_out/Supp14_glms_on_strain_by_mutagen_tolerance.csv", 
            sep = ",")
