
library(ggplot2) # data visualization 
library(viridis) # data visualization 
library(ggh4x) # data visualization 
library(dplyr) # data manipulation
library(effsize) # stats
library(reshape2)
library(readr)
library(gridExtra)
library(cowplot)
library(desiderata)

########### tag mut control ###########
filepath<- "../../symbulation_experiments/2025.05.02_control_fixed_mutualists/"
a <- read_csv(paste0(filepath,"org_dump_file.dat"),col_names=T)
a$sym_tag[a$sym_tag==""]<-NA
a$cond<-"mut"# <- org_dump_df

filepath<- "../../symbulation_experiments/2025.05.02_control_fixed_parasites/"
b <- read_csv(paste0(filepath,"org_dump_file.dat"),col_names=T)
b$sym_tag[b$sym_tag==""]<-NA
b$cond<-"para"# <- org_dump_df

small <- rbind(subset(b, !is.na(sym_tag)), subset(a, !is.na(sym_tag)))
small <- subset(small, tag_mut > 0)
identifying_cols <- 1:3
for (x in identifying_cols){
  small[[x]]<- as.factor(small[[x]])
}

if("_none" %in% small$tag_dist){
  small$tag_dist <- factor(small$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625"))
} else if("-1" %in% small$tag_dist){
  small$tag_dist <- factor(small$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"-1"))
  
}



small$sym_bin <- ifelse(small$sym_int < -0.6, "-1 to -0.6 (Parasitic)",
                        ifelse(small$sym_int < -0.2 & small$sym_int >= -0.6, "-0.6 to -0.2 (Detrimental)",
                               ifelse(small$sym_int < 0.2 & small$sym_int >= -0.2, "-0.2 to 0.2 (Nearly neutral)",
                                      ifelse(small$sym_int < 0.6 & small$sym_int >= 0.2, "0.2 to 0.6 (Positive)",
                                             ifelse(small$sym_int <= 1 & small$sym_int >= 0.6, "0.6 to 1.0 (Mutualistic)", "error"
                                             )))))


facet_nested_theme <- theme(
  strip.background = element_rect(fill = "white", colour = "grey", linetype="dotted", size = 0.2, linewidth=30), 
  panel.background = element_rect(fill='white', color='grey'),
  panel.border = element_blank(),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  legend.key = element_rect(color = "transparent"),
  panel.spacing.x = unit(0,"line"),
  panel.spacing.y = unit(0,"line"),
  axis.text.x = element_text(angle = 60, vjust = 0.2, hjust=0.2)) 

small <- small %>% group_by(sym_bin,tag_dist,tag_mut) %>% mutate(td_mean =  mean(tag_distance))


plot <- ggplot(small, aes(x=tag_distance)) + 
  geom_histogram(data=subset(small, sym_bin=="-1 to -0.6 (Parasitic)"), aes(fill="Fixed parasites",y=..count../10000),alpha=0.7,  bins=16)+
  geom_histogram(data=subset(small, sym_bin=="0.6 to 1.0 (Mutualistic)"), aes(fill="Fixed mutualists",y=..count../10000), alpha=0.7, bins=16) +
  geom_vline(data=subset(small, sym_bin=="0.6 to 1.0 (Mutualistic)"), aes(xintercept = td_mean), color="#7A0403FF",linetype="dotted") +
  geom_vline(data=subset(small, sym_bin=="-1 to -0.6 (Parasitic)"), aes(xintercept = td_mean), color="#3E9BFEFF",linetype="dotted") +
  xlab("Tag proportional mismatch") + 
  scale_fill_manual(name="Symbiont behavior",
                    values=c("Fixed parasites"="#3E9BFEFF", "Fixed mutualists"="#7A0403FF"),
                    breaks=c("Fixed parasites", "Fixed mutualists")) + 
  ylab("Pooled count across replicates / 10,000")+
  facet_nested("Tag permissiveness" + tag_dist ~ "Tag mutation rate" + tag_mut, drop=TRUE) + 
  rotate_y_facet_text(angle = 0) +
  facet_nested_theme + theme(legend.position = 'bottom') #+ theme(axis.text.y = element_text(angle = 60))
plot

unique(subset(small, tag_dist != 0.125 & tag_dist != 0.1875 & tag_dist != 0.25)$tag_dist)
plot <- ggplot(subset(small, tag_dist != 0.125 & tag_dist != 0.1875 & tag_dist != 0.25), aes(x=tag_distance)) + 
  geom_histogram(data=subset(small, tag_dist != 0.125 & tag_dist != 0.1875 & tag_dist != 0.25 & sym_bin=="-1 to -0.6 (Parasitic)"), aes(fill="Fixed parasites",y=..count../10000),alpha=0.7,  bins=16)+
  geom_histogram(data=subset(small, tag_dist != 0.125 & tag_dist != 0.1875 & tag_dist != 0.25 & sym_bin=="0.6 to 1.0 (Mutualistic)"), aes(fill="Fixed mutualists",y=..count../10000), alpha=0.7, bins=16) +
  geom_vline(data=subset(small, tag_dist != 0.125 & tag_dist != 0.1875 & tag_dist != 0.25 & sym_bin=="0.6 to 1.0 (Mutualistic)"), aes(xintercept = td_mean), color="#7A0403FF",linetype="dotted") +
  geom_vline(data=subset(small, tag_dist != 0.125 & tag_dist != 0.1875 & tag_dist != 0.25 & sym_bin=="-1 to -0.6 (Parasitic)"), aes(xintercept = td_mean), color="#3E9BFEFF",linetype="dotted") +
  xlab("Tag proportional mismatch") + 
  scale_fill_manual(name="Symbiont behavior",
                    values=c("Fixed parasites"="#3E9BFEFF", "Fixed mutualists"="#7A0403FF"),
                    breaks=c("Fixed parasites", "Fixed mutualists")) + 
  ylab("Pooled count across replicates / 10,000")+
  facet_nested("ab" + tag_dist ~ "Tag mutation rate" + tag_mut, drop=TRUE) + 
  rotate_y_facet_text(angle = 0) +
  facet_nested_theme + theme(legend.position = 'bottom') +
  scale_y_continuous(breaks = c(0,5))
plot
 parastab_control_small <- small

########### vt control ###########
filepath<- "../../symbulation_experiments/2025.05.04_control_vtsweep_fixed_mutualists/"
a <- read_csv(paste0(filepath,"org_dump_file.dat"),col_names=T)
a$sym_tag[a$sym_tag==""]<-NA
a$cond<-"mut"# <- org_dump_df

filepath<- "../../symbulation_experiments/2025.05.04_control_vtsweep_fixed_parasites/"
b <- read_csv(paste0(filepath,"org_dump_file.dat"),col_names=T)
b$sym_tag[b$sym_tag==""]<-NA
b$cond<-"para"# <- org_dump_df

small <- rbind(subset(b, !is.na(sym_tag)), subset(a, !is.na(sym_tag)))

identifying_cols <- 1:3
for (x in identifying_cols){
  small[[x]]<- as.factor(small[[x]])
}

if("_none" %in% small$tag_dist){
  small$tag_dist <- factor(small$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625"))
} else if("-1" %in% small$tag_dist){
  small$tag_dist <- factor(small$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"-1"))
  
}



small$sym_bin <- ifelse(small$sym_int < -0.6, "-1 to -0.6 (Parasitic)",
                        ifelse(small$sym_int < -0.2 & small$sym_int >= -0.6, "-0.6 to -0.2 (Detrimental)",
                               ifelse(small$sym_int < 0.2 & small$sym_int >= -0.2, "-0.2 to 0.2 (Nearly neutral)",
                                      ifelse(small$sym_int < 0.6 & small$sym_int >= 0.2, "0.2 to 0.6 (Positive)",
                                             ifelse(small$sym_int <= 1 & small$sym_int >= 0.6, "0.6 to 1.0 (Mutualistic)", "error"
                                             )))))


facet_nested_theme <- theme(
  strip.background = element_rect(fill = "white", colour = "grey", linetype="dotted", size = 0.2, linewidth=30), 
  panel.background = element_rect(fill='white', color='grey'),
  panel.border = element_blank(),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  legend.key = element_rect(color = "transparent"),
  panel.spacing.x = unit(0,"line"),
  panel.spacing.y = unit(0,"line"),
  axis.text.x = element_text(angle = 60, vjust = 0.2, hjust=0.2))

small <- small %>% group_by(sym_bin,tag_dist,vt) %>% mutate(td_mean =  mean(as.numeric(tag_distance)))

plot <- ggplot(subset(small,  tag_dist != 0.125 & tag_dist != 0.1875 & tag_dist != 0.25), aes(x=as.numeric(tag_distance))) + 
  geom_histogram(data=subset(small,  tag_dist != 0.125 & tag_dist != 0.1875 & tag_dist != 0.25 & sym_bin=="-1 to -0.6 (Parasitic)"), aes(fill="Fixed parasites",y=..count../10000),alpha=0.7,  bins=16)+
  geom_histogram(data=subset(small,  tag_dist != 0.125 & tag_dist != 0.1875 & tag_dist != 0.25 & sym_bin=="0.6 to 1.0 (Mutualistic)"), aes(fill="Fixed mutualists",y=..count../10000), alpha=0.7, bins=16) +
  geom_vline(data=subset(small,  tag_dist != 0.125 & tag_dist != 0.1875 & tag_dist != 0.25 & sym_bin=="0.6 to 1.0 (Mutualistic)"), aes(xintercept = td_mean), color="#7A0403FF",linetype="dotted") +
  geom_vline(data=subset(small,  tag_dist != 0.125 & tag_dist != 0.1875 & tag_dist != 0.25 & sym_bin=="-1 to -0.6 (Parasitic)"), aes(xintercept = td_mean), color="#3E9BFEFF",linetype="dotted") +
  xlab("Tag proportional mismatch") + 
  scale_fill_manual(name="Symbiont behavior",
                    values=c("Fixed parasites"="#3E9BFEFF", "Fixed mutualists"="#7A0403FF"),
                    breaks=c("Fixed parasites", "Fixed mutualists")) + 
  ylab("Pooled count across replicates / 10,000")+
  facet_nested("ab" + tag_dist ~ "Vertical transmission rate" + vt, drop=TRUE) + 
  rotate_y_facet_text(angle = 0) + 
  facet_nested_theme + theme(legend.position = 'bottom') + 
  scale_y_continuous(breaks = c(0,4,8)) + 
  scale_x_continuous(breaks = c(0,0.4,0.8))

plot
#vtsweep_control_small <- small

########### tag mut evo ###########
filepath<- "../../symbulation_experiments/2025.04.17_parastab_tagmat_fixed_phylointeractions/"
small <- read_csv(paste0(filepath,"inc_syms_only_dump_file.dat"),col_names=T)
small$sym_tag[small$sym_tag==""]<-NA
small$tag_distance <- as.numeric(small$tag_distance)
small <- subset(small, !is.na(sym_tag))
identifying_cols <- 1:3
for (x in identifying_cols){
  small[[x]]<- as.factor(small[[x]])
}

if("_none" %in% small$tag_dist){
  small$tag_dist <- factor(small$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625"))
} else if("-1" %in% small$tag_dist){
  small$tag_dist <- factor(small$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"-1"))
  
}
small$sym_bin <- ifelse(small$sym_int < -0.6, "-1 to -0.6 (Parasitic)",
                        ifelse(small$sym_int < -0.2 & small$sym_int >= -0.6, "-0.6 to -0.2 (Detrimental)",
                       ifelse(small$sym_int < 0.2 & small$sym_int >= -0.2, "-0.2 to 0.2 (Nearly neutral)",
                      ifelse(small$sym_int < 0.6 & small$sym_int >= 0.2, "0.2 to 0.6 (Positive)",
                     ifelse(small$sym_int <= 1 & small$sym_int >= 0.6, "0.6 to 1.0 (Mutualistic)", "error"
                     )))))

plot <- ggplot(small, aes(x=tag_distance)) + 
  geom_histogram(data=subset(small, sym_bin=="-1 to -0.6 (Parasitic)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16)+
  geom_histogram(data=subset(small, sym_bin=="-0.6 to -0.2 (Detrimental)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
  geom_histogram(data=subset(small, sym_bin=="-0.2 to 0.2 (Nearly neutral)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
   geom_histogram(data=subset(small, sym_bin=="0.2 to 0.6 (Positive)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
  geom_histogram(data=subset(small, sym_bin=="0.6 to 1.0 (Mutualistic)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
  xlab("Tag proportional mismatch") + 
  scale_fill_manual(name="Symbiont\nbehavior",
                    values=c("-1 to -0.6 (Parasitic)"="#3E9BFEFF", '-0.6 to -0.2 (Detrimental)'="#46F884FF", 
                             '-0.2 to 0.2 (Nearly neutral)'="#E1DD37FF",
                             "0.2 to 0.6 (Positive)"="#F05B12FF", "0.6 to 1.0 (Mutualistic)"="#7A0403FF"),
                    breaks=c("-1 to -0.6 (Parasitic)",'-0.6 to -0.2 (Detrimental)',
                             '-0.2 to 0.2 (Nearly neutral)',
                             "0.2 to 0.6 (Positive)", "0.6 to 1.0 (Mutualistic)")) + 
  ylab("Pooled count across replicates / 10,000")+
  facet_nested("a" + tag_dist ~ "Tag mutation rate" + tag_mut, drop=TRUE) + 
  rotate_y_facet_text(angle = 0) + 
  scale_y_continuous(breaks = c(0,5)) + 
  facet_nested_theme+ theme(legend.position = 'bottom')+guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.key.size = unit(0.4, 'cm'), #change legend key size
        #legend.key.height = unit(0.1, 'cm'), #change legend key height
        #legend.key.width = unit(0.1, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8)) #change legend text font size
plot

tagmut_evo_small <- small





########### vt evo ########### 
filepath<- "../../symbulation_experiments/2025.04.17_vtsweep_tagmat_fixed_phylointeractions/"
small <- read_csv(paste0(filepath,"inc_syms_only_dump_file.dat"),col_names=T)
colnames(small)
small <- subset(small, vt <= 0.5)
small$sym_tag[small$sym_tag==""]<-NA
small$tag_distance <- as.numeric(small$tag_distance)
small <- subset(small, !is.na(sym_tag))
identifying_cols <- 1:3
for (x in identifying_cols){
  small[[x]]<- as.factor(small[[x]])
}

if("_none" %in% small$tag_dist){
  small$tag_dist <- factor(small$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625"))
} else if("-1" %in% small$tag_dist){
  small$tag_dist <- factor(small$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"-1"))
  
}
small$sym_bin <- ifelse(small$sym_int < -0.6, "-1 to -0.6 (Parasitic)",
                        ifelse(small$sym_int < -0.2 & small$sym_int >= -0.6, "-0.6 to -0.2 (Detrimental)",
                               ifelse(small$sym_int < 0.2 & small$sym_int >= -0.2, "-0.2 to 0.2 (Nearly neutral)",
                                      ifelse(small$sym_int < 0.6 & small$sym_int >= 0.2, "0.2 to 0.6 (Positive)",
                                             ifelse(small$sym_int <= 1 & small$sym_int >= 0.6, "0.6 to 1.0 (Mutualistic)", "error"
                                             )))))


plot <- ggplot(small, aes(x=tag_distance)) + 
  geom_histogram(data=subset(small, sym_bin=="-1 to -0.6 (Parasitic)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16)+
  geom_histogram(data=subset(small, sym_bin=="-0.6 to -0.2 (Detrimental)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
  geom_histogram(data=subset(small, sym_bin=="-0.2 to 0.2 (Nearly neutral)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
  geom_histogram(data=subset(small, sym_bin=="0.2 to 0.6 (Positive)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
  geom_histogram(data=subset(small, sym_bin=="0.6 to 1.0 (Mutualistic)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
  xlab("Tag proportional mismatch") + 
  scale_fill_manual(name="Symbiont behavior",
                    values=c("-1 to -0.6 (Parasitic)"="#3E9BFEFF", '-0.6 to -0.2 (Detrimental)'="#46F884FF", 
                             '-0.2 to 0.2 (Nearly neutral)'="#E1DD37FF",
                             "0.2 to 0.6 (Positive)"="#F05B12FF", "0.6 to 1.0 (Mutualistic)"="#7A0403FF"),
                    breaks=c("-1 to -0.6 (Parasitic)",'-0.6 to -0.2 (Detrimental)',
                             '-0.2 to 0.2 (Nearly neutral)',
                             "0.2 to 0.6 (Positive)", "0.6 to 1.0 (Mutualistic)")) + 
  ylab("Pooled count across replicates / 10,000")+
  facet_nested("a" + tag_dist ~ "Vertical transmission rate" + vt, drop=TRUE) + 
  facet_nested_theme+ theme(legend.position = 'bottom')+guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  rotate_y_facet_text(angle = 0) + 
  scale_y_continuous(breaks = c(0,5)) + 
  facet_nested_theme+ theme(legend.position = 'bottom')+guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.key.size = unit(0.4, 'cm'), #change legend key size
        #legend.key.height = unit(0.1, 'cm'), #change legend key height
        #legend.key.width = unit(0.1, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8)) #change legend text font size
plot



########### vt strat barplot ########### 
filepath <- "../../symbulation_experiments/2025.04.17_vtsweep_tagmat_fixed_phylointeractions/"
syms_df <- read_csv(paste0(filepath,"sym_vals.dat"))
syms_df <- subset(syms_df, update==max(update))
hosts_df <- read_csv(paste0(filepath,"host_vals.dat"))
hosts_df <- subset(hosts_df, update==max(update))

hosts_df$tag_dist <- ifelse(hosts_df$tag_dist == '_none', "no tag",hosts_df$tag_dist)
syms_df$tag_dist <- ifelse(syms_df$tag_dist == '_none', "no tag",syms_df$tag_dist)

identifying_cols <- 1:3
for (x in identifying_cols){
  syms_df[[x]]<- as.factor(syms_df[[x]])
  hosts_df[[x]]<- as.factor(hosts_df[[x]])
}

str(hosts_df)
hosts_df$tag_dist <- factor(hosts_df$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"no tag"))
syms_df$tag_dist <- factor(syms_df$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"no tag"))

rm(identifying_cols)

end_syms_df <- subset(syms_df, update==max(update))

end_syms_df <- end_syms_df %>%
  group_by(tag_dist,vt) %>% 
  mutate(com_rep = match(seed, unique(seed)))

start_i <- 9
colnames(end_syms_df)[start_i:(start_i+19)]
hist_df <- end_syms_df %>% select(c(com_rep, tag_dist, vt), colnames(end_syms_df)[start_i:(start_i+19)])

# Plotting
# load our color palette 
tenhelix <- c("#891901", "#B50142", "#D506AD", "#AB08FF", "#5731FD", "#4755FF", "#42bcf5", "#42e0f5", "#9af1fc", "#c5f6fc")

# now let's collapse all those columns...
hist_df <- melt(hist_df, id=c("com_rep","vt","tag_dist"))

# Group bins into larger categories and give them descriptive names
hist_df$variable <- as.character(hist_df$variable)

hist_df$variable[hist_df$variable == 'Hist_-1'] <- "-1 to -0.8 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.9'] <- "-1 to -0.8 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.8'] <- "-0.8 to -0.6 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.7'] <- "-0.8 to -0.6 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.6'] <- "-0.6 to -0.4 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.5'] <- "-0.6 to -0.4 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.4'] <- "-0.4 to -0.2 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.3'] <- "-0.4 to -0.2 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.2'] <- "-0.2 to 0 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_-0.1'] <- "-0.2 to 0 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.0'] <- "0 to 0.2 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.1'] <- "0 to 0.2 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.2'] <- "0.2 to 0.4 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.3'] <- "0.2 to 0.4 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.4'] <- "0.4 to 0.6 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.5'] <- "0.4 to 0.6 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.6'] <- "0.6 to 0.8 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.7'] <- "0.6 to 0.8 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.8'] <- "0.8 to 1.0 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.9'] <- "0.8 to 1.0 (Mutualistic)"

# force set the order of the variable column
hist_df$variable <- factor(hist_df$variable, levels=rev(unique(hist_df$variable)))

# compress bin counts
hist_df <- hist_df %>%
  group_by(com_rep, tag_dist, vt, variable) %>% 
  summarise(value = sum(value))

hist_df$tag_dist <- factor(hist_df$tag_dist, levels=rev(unique(hist_df$tag_dist)))

# order by proportion of ubermutualists
hist_df$max_mut_count <- rep(hist_df$value[hist_df$variable=="0.8 to 1.0 (Mutualistic)"],1,each=10)
hist_df <- hist_df %>% group_by(tag_dist,vt,com_rep) %>% mutate(max_mut_prop = max_mut_count/sum(value))
hist_df<-hist_df %>% group_by(tag_dist, vt) %>% arrange(desc(max_mut_prop)) 
hist_df<-hist_df %>% group_by(tag_dist, vt) %>%mutate(group_id = rep(1:30,1,each=10))

a <- ggplot(hist_df, aes(fill=variable, x=group_id,y=value), na.rm=T) + 
  geom_col(position="fill", na.rm=T, width = 1) +
  scale_fill_manual(name="Behavior",values=tenhelix,guide="none")+
  xlab("Replicate") + ylab("Proportion of symbiont population") +
  facet_nested("a" + tag_dist ~ "Vertical transmission rate" + vt) + 
  facet_nested_theme + rotate_y_facet_text(angle = 0) 

syms_hist <- hist_df
end_syms_df <- subset(hosts_df, update==max(update))

end_syms_df <- end_syms_df %>%
  group_by(tag_dist,vt) %>% 
  mutate(com_rep = match(seed, unique(seed)))

start_i <- 10
colnames(end_syms_df)[start_i:(start_i+19)]
hist_df <- end_syms_df %>% select(c(com_rep, tag_dist, vt), colnames(end_syms_df)[start_i:(start_i+19)])

# Plotting
# load our color palette 
tenhelix <- c("#891901", "#B50142", "#D506AD", "#AB08FF", "#5731FD", "#4755FF", "#42bcf5", "#42e0f5", "#9af1fc", "#c5f6fc")

# now let's collapse all those columns...
hist_df <- melt(hist_df, id=c("com_rep","vt","tag_dist"))

# Group bins into larger categories and give them descriptive names
hist_df$variable <- as.character(hist_df$variable)

hist_df$variable[hist_df$variable == 'Hist_-1'] <- "-1 to -0.8 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.9'] <- "-1 to -0.8 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.8'] <- "-0.8 to -0.6 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.7'] <- "-0.8 to -0.6 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.6'] <- "-0.6 to -0.4 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.5'] <- "-0.6 to -0.4 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.4'] <- "-0.4 to -0.2 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.3'] <- "-0.4 to -0.2 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.2'] <- "-0.2 to 0 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_-0.1'] <- "-0.2 to 0 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.0'] <- "0 to 0.2 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.1'] <- "0 to 0.2 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.2'] <- "0.2 to 0.4 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.3'] <- "0.2 to 0.4 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.4'] <- "0.4 to 0.6 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.5'] <- "0.4 to 0.6 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.6'] <- "0.6 to 0.8 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.7'] <- "0.6 to 0.8 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.8'] <- "0.8 to 1.0 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.9'] <- "0.8 to 1.0 (Mutualistic)"

# force set the order of the variable column
hist_df$variable <- factor(hist_df$variable, levels=rev(unique(hist_df$variable)))

# compress bin counts
hist_df <- hist_df %>%
  group_by(com_rep, tag_dist, vt, variable) %>% 
  summarise(value = sum(value))

hist_df$tag_dist <- factor(hist_df$tag_dist, levels=rev(unique(hist_df$tag_dist)))

hist_df <- hist_df %>%left_join(select(syms_hist, c(tag_dist, vt, group_id, variable, com_rep)), 
                                by= join_by(vt,variable,tag_dist,com_rep))

b <- ggplot(hist_df, aes(fill=variable, x=group_id,y=value), na.rm=T) + 
  geom_col(position="fill", na.rm=T, width = 1) +
  scale_fill_manual(name="Behavior",values=tenhelix)+
  xlab("Replicate") + ylab("Proportion of host population") +
  facet_nested("a" + tag_dist ~ "Vertical transmission rate" + vt) + 
  facet_nested_theme + rotate_y_facet_text(angle = 0) 

legend <- cowplot::get_legend(b+ theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))
plots <- grid.arrange(a+theme(axis.text.y=element_blank(), axis.text.x=element_blank(), axis.ticks = element_blank()),b+theme(legend.position = "none",axis.text.y=element_blank(), axis.text.x=element_blank(),axis.ticks = element_blank()), ncol=2,widths=c(1,1))
grid.arrange(plots,legend, ncol=1, nrow=2,heights=c(1,0.2))

plots


########### tag mut strat barplot ########### 
options(scipen = 999)
filepath <- "../../symbulation_experiments/2025.04.17_parastab_tagmat_fixed_phylointeractions/"
syms_df <- read_csv(paste0(filepath,"sym_vals.dat"))
syms_df <- subset(syms_df, update==max(update))
hosts_df <- read_csv(paste0(filepath,"host_vals.dat"))
hosts_df <- subset(hosts_df, update==max(update))

a <- syms_df 
b <- hosts_df

print((hosts_df$tag_mut))
hosts_df$tag_dist <- as.character(hosts_df$tag_dist)
syms_df$tag_dist <- as.character(syms_df$tag_dist)
hosts_df$tag_mut <- as.character(hosts_df$tag_mut)
syms_df$tag_mut <- as.character(syms_df$tag_mut)

hosts_df$tag_dist <- ifelse(hosts_df$tag_dist == '-1', "no tag",hosts_df$tag_dist)
syms_df$tag_dist <- ifelse(syms_df$tag_dist == '-1', "no tag",syms_df$tag_dist)

hosts_df$tag_mut <- ifelse(hosts_df$tag_mut == '-1', "no tag",hosts_df$tag_mut)
syms_df$tag_mut <- ifelse(syms_df$tag_mut == '-1', "no tag",syms_df$tag_mut)

identifying_cols <- 1:3
for (x in identifying_cols){
  syms_df[[x]]<- as.factor(syms_df[[x]])
  hosts_df[[x]]<- as.factor(hosts_df[[x]])
}

str(hosts_df$tag_dist)

hosts_df$tag_dist <- factor(hosts_df$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"no tag"))
syms_df$tag_dist <- factor(syms_df$tag_dist, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"no tag"))


hosts_df$tag_mut <- factor(hosts_df$tag_mut, levels=c("no tag", "0", "0.0001", "0.0005", "0.001", "0.005", "0.01", "0.05", "0.1"))
syms_df$tag_mut <- factor(syms_df$tag_mut, levels=c("no tag", "0", "0.0001", "0.0005", "0.001", "0.005", "0.01", "0.05", "0.1"))



rm(identifying_cols)

end_syms_df <- subset(syms_df, update==max(update))

end_syms_df <- end_syms_df %>%
  group_by(tag_dist,tag_mut) %>% 
  mutate(com_rep = match(seed, unique(seed)))

start_i <- 9
colnames(end_syms_df)[start_i:(start_i+19)]
hist_df <- end_syms_df %>% select(c(com_rep, tag_dist, tag_mut), colnames(end_syms_df)[start_i:(start_i+19)])

# Plotting
# load our color palette 
tenhelix <- c("#891901", "#B50142", "#D506AD", "#AB08FF", "#5731FD", "#4755FF", "#42bcf5", "#42e0f5", "#9af1fc", "#c5f6fc")

# now let's collapse all those columns...
hist_df <- melt(hist_df, id=c("com_rep","tag_mut","tag_dist"))

# Group bins into larger categories and give them descriptive names
hist_df$variable <- as.character(hist_df$variable)

hist_df$variable[hist_df$variable == 'Hist_-1'] <- "-1 to -0.8 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.9'] <- "-1 to -0.8 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.8'] <- "-0.8 to -0.6 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.7'] <- "-0.8 to -0.6 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.6'] <- "-0.6 to -0.4 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.5'] <- "-0.6 to -0.4 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.4'] <- "-0.4 to -0.2 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.3'] <- "-0.4 to -0.2 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.2'] <- "-0.2 to 0 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_-0.1'] <- "-0.2 to 0 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.0'] <- "0 to 0.2 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.1'] <- "0 to 0.2 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.2'] <- "0.2 to 0.4 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.3'] <- "0.2 to 0.4 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.4'] <- "0.4 to 0.6 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.5'] <- "0.4 to 0.6 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.6'] <- "0.6 to 0.8 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.7'] <- "0.6 to 0.8 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.8'] <- "0.8 to 1.0 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.9'] <- "0.8 to 1.0 (Mutualistic)"

# force set the order of the variable column
hist_df$variable <- factor(hist_df$variable, levels=rev(unique(hist_df$variable)))

# compress bin counts
hist_df <- hist_df %>%
  group_by(com_rep, tag_dist, tag_mut, variable) %>% 
  summarise(value = sum(value))

hist_df$tag_dist <- factor(hist_df$tag_dist, levels=rev(unique(hist_df$tag_dist)))

# order by proportion of ubermutualists
hist_df$max_mut_count <- rep(hist_df$value[hist_df$variable=="0.8 to 1.0 (Mutualistic)"],1,each=10)
hist_df <- hist_df %>% group_by(tag_dist,tag_mut,com_rep) %>% mutate(max_mut_prop = max_mut_count/sum(value))
hist_df<-hist_df %>% group_by(tag_dist, tag_mut) %>% arrange(desc(max_mut_prop)) 
hist_df<-hist_df %>% group_by(tag_dist, tag_mut) %>%mutate(group_id = rep(1:30,1,each=10))

a <- ggplot(hist_df, aes(fill=variable, x=group_id,y=value), na.rm=T) + 
  geom_col(position="fill", na.rm=T, width = 1) +
  scale_fill_manual(name="Behavior",values=tenhelix,guide="none")+
  xlab("Replicate") + ylab("Proportion of symbiont population") +
  facet_nested("a" + tag_dist ~ "Tag mutation rate" + tag_mut) + 
  facet_nested_theme+ rotate_y_facet_text(angle = 0) 

syms_hist <- hist_df
end_syms_df <- subset(hosts_df, update==max(update))

end_syms_df <- end_syms_df %>%
  group_by(tag_dist,tag_mut) %>% 
  mutate(com_rep = match(seed, unique(seed)))

start_i <- 10
colnames(end_syms_df)[start_i:(start_i+19)]
hist_df <- end_syms_df %>% select(c(com_rep, tag_dist, tag_mut), colnames(end_syms_df)[start_i:(start_i+19)])

# Plotting
# load our color palette 
tenhelix <- c("#891901", "#B50142", "#D506AD", "#AB08FF", "#5731FD", "#4755FF", "#42bcf5", "#42e0f5", "#9af1fc", "#c5f6fc")

# now let's collapse all those columns...
hist_df <- melt(hist_df, id=c("com_rep","tag_mut","tag_dist"))

# Group bins into larger categories and give them descriptive names
hist_df$variable <- as.character(hist_df$variable)

hist_df$variable[hist_df$variable == 'Hist_-1'] <- "-1 to -0.8 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.9'] <- "-1 to -0.8 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.8'] <- "-0.8 to -0.6 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.7'] <- "-0.8 to -0.6 (Parasitic)"
hist_df$variable[hist_df$variable == 'Hist_-0.6'] <- "-0.6 to -0.4 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.5'] <- "-0.6 to -0.4 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.4'] <- "-0.4 to -0.2 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.3'] <- "-0.4 to -0.2 (Detrimental)"
hist_df$variable[hist_df$variable == 'Hist_-0.2'] <- "-0.2 to 0 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_-0.1'] <- "-0.2 to 0 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.0'] <- "0 to 0.2 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.1'] <- "0 to 0.2 (Nearly Neutral)"
hist_df$variable[hist_df$variable == 'Hist_0.2'] <- "0.2 to 0.4 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.3'] <- "0.2 to 0.4 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.4'] <- "0.4 to 0.6 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.5'] <- "0.4 to 0.6 (Positive)"
hist_df$variable[hist_df$variable == 'Hist_0.6'] <- "0.6 to 0.8 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.7'] <- "0.6 to 0.8 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.8'] <- "0.8 to 1.0 (Mutualistic)"
hist_df$variable[hist_df$variable == 'Hist_0.9'] <- "0.8 to 1.0 (Mutualistic)"

# force set the order of the variable column
hist_df$variable <- factor(hist_df$variable, levels=rev(unique(hist_df$variable)))

# compress bin counts
hist_df <- hist_df %>%
  group_by(com_rep, tag_dist, tag_mut, variable) %>% 
  summarise(value = sum(value))

hist_df$tag_dist <- factor(hist_df$tag_dist, levels=rev(unique(hist_df$tag_dist)))

hist_df <- hist_df %>%left_join(select(syms_hist, c(tag_dist, tag_mut, group_id, variable, com_rep)), 
                                by= join_by(tag_mut,variable,tag_dist,com_rep))

b <- ggplot(hist_df, aes(fill=variable, x=group_id,y=value), na.rm=T) + 
  geom_col(position="fill", na.rm=T, width = 1) +
  scale_fill_manual(name="Behavior",values=tenhelix)+
  xlab("Replicate") + ylab("Proportion of host population") +
  facet_nested("a" + tag_dist ~ "Tag mutation rate" + tag_mut) + 
  facet_nested_theme + rotate_y_facet_text(angle = 0) 

legend <- cowplot::get_legend(b+ theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))
plots <- grid.arrange(a+theme(axis.text.y=element_blank(), axis.text.x=element_blank(), axis.ticks = element_blank()),b+theme(legend.position = "none",axis.text.y=element_blank(), axis.text.x=element_blank(),axis.ticks = element_blank()), ncol=2,widths=c(1,1))
grid.arrange(plots,legend, ncol=1, nrow=2,heights=c(1,0.2))

