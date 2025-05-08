setwd("~/ALIFE2025-TagPartnerChoice/stats")

filepath <- "../data/tagmut_evolving/"

syms_df <- read.csv(paste0(filepath,"sym_vals.dat"),h=T)

syms_df <- subset(syms_df, update==max(update))

start_i <- 9
colnames(syms_df[(start_i+16):(start_i+19)])
syms_df$mut_present <- rowSums(syms_df[(start_i+16):(start_i+19)]) > 0
syms_df <- syms_df %>% group_by(tag_dist,tag_mut) %>% summarize(reps_muts_present = sum(mut_present))
syms_df$reps_muts_absent <- 30 - syms_df$reps_muts_present
no_tags <- subset(syms_df, tag_dist=="-1")
syms_df <- subset(syms_df, tag_dist!="-1")

nm <- c()
md <- c()
low <- c()
high <- c()
n_test <- 0
for (td in unique(syms_df$tag_dist)){
  for(v in unique(syms_df$tag_mut)){
    n_test <- n_test + 1
    print(paste0("td",td,"_tm",v))
    data = matrix(c(as.numeric(subset(syms_df, tag_dist == td & tag_mut == v)[3:4]),as.numeric(no_tags[3:4])), nrow = 2)
    test <- fisher.test(data)
    nm <- append(nm, paste0("td",td,"_tm",v))
    md <- append(md, test$p.value)
    low <- append(low, test$conf.int[1])
    high <- append(high, test$conf.int[2])
  }
}

models <- data.frame(names = nm, p_vals=md, conf_low = low, conf_high = high)
write.csv(models, "tagmut_evolve_mutpres", row.names=F)






######################### VT ######################### 

filepath <- "../data/vt_evolving/"

syms_df <- read.csv(paste0(filepath,"sym_vals.dat"),h=T)

syms_df <- subset(syms_df, update==max(update))

start_i <- 9
colnames(syms_df[(start_i+16):(start_i+19)])
syms_df$mut_present <- rowSums(syms_df[(start_i+16):(start_i+19)]) > 0
syms_df <- syms_df %>% group_by(tag_dist,vt) %>% summarize(reps_muts_present = sum(mut_present))
syms_df$reps_muts_absent <- 30 - syms_df$reps_muts_present
no_tags <- subset(syms_df, tag_dist=="-1")
syms_df <- subset(syms_df, tag_dist!="-1")

nm <- c()
md <- c()
low <- c()
high <- c()
n_test <- 0
for (td in unique(syms_df$tag_dist)){
  for(v in unique(syms_df$vt)){
    n_test <- n_test + 1
    print(paste0("td",td,"_tm",v))
    data = matrix(c(as.numeric(subset(syms_df, tag_dist == td & vt == v)[3:4]),as.numeric(no_tags[3:4])), nrow = 2)
    test <- fisher.test(data)
    nm <- append(nm, paste0("td",td,"_vt",v))
    md <- append(md, test$p.value)
    low <- append(low, test$conf.int[1])
    high <- append(high, test$conf.int[2])
  }
}

models <- data.frame(names = nm, p_vals=md, conf_low = low, conf_high = high)
write.csv(models, "tagmut_evolve_mutpres", row.names=F)

