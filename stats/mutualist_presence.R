setwd("~/ALIFE2025-TagPartnerChoice/stats")

filepath <- "../data/tagmut_evolving/"

syms_df <- read.csv(paste0(filepath,"sym_vals.dat"),h=T)

syms_df <- subset(syms_df, update==max(update))
syms_df$mut_present <- rowSums(syms_df[(start_i+16):(start_i+19)]) > 0
syms_df <- syms_df %>% group_by(tag_dist,tag_mut) %>% summarize(reps_muts_present = sum(mut_present))
syms_df$reps_muts_absent <- 30 - syms_df$reps_muts_present
no_tags <- subset(syms_df, tag_dist=="-1")


nm <- c()
md <- c()
effect <- c()
n_test <- 0
for (td in unique(small$tag_dist)){
  for(v in unique(small$tag_mut)){
    a <- 
    data = matrix(c(subset(small, tag_dist == td & tag_mut == v)[3:4],no_tags[3:4]), nrow = 2)
    if("para" %in% a$cond){
      n_test <- n_test + 1
      nm <- append(nm, paste0("td",td,"_tm",v))
      md <- append(md, wilcox.test(tag_distance ~ cond, data=a)$p.value)
      
      print("\n*******************")
      print(nm[n_test])
      print(md[n_test])
    }
  }
}
print(n_test)


models <- data.frame(names = nm, p_vals=md)
write.csv(models, "parastab_evolve_paravsmut_pvals.csv", row.names=F)