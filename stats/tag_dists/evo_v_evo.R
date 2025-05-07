
library(dplyr)
library(ggplot2)
library(viridis)
library(effsize) # stats

print("about to read mut")

filepath<- ""
small <- read.csv(paste0(filepath,"inc_syms_only_dump_file.dat"),h=T,
                                    #seed,tag_dist,tag_mut
                        colClasses=c("factor","factor","factor",
                                     #host_int,sym_int
                                    "double","double",
                                    # repro/to/from counts
                                    "double","double","double","double","double","double",
                                    #host_tag,sym_tag,tag_distance
                                    "character","character","double"))

small <- subset(small, sym_int >= 0.6 | sym_int < -0.6)
small$cond <- "mut"
small$cond <-ifelse(small$sym_int < -0.6, "para", "mut")

print("done binding and deleting; about to transform")
identifying_cols <- 1:3
for (x in identifying_cols){
  small[[x]]<- as.factor(small[[x]])
}
small <- small %>% group_by(tag_dist,tag_mut) %>% mutate(para_pres = "para"%in%cond)
small <- small %>% group_by(tag_dist,tag_mut) %>% mutate(mut_pres = "mut"%in%cond)
small <- subset(small, para_pres == T & mut_pres == T & tag_mut != '0')
small$cond <- as.factor(small$cond)

print("done transforming; about to calculate stats")

nm <- c()
md <- c()
effect <- c()
n_test <- 0
for (td in unique(small$tag_dist)){
  for(v in unique(small$tag_mut)){
    a <- subset(small, tag_dist == td & tag_mut == v)
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

print("done calculating p vals!")

nm <- c()
md <- c()
effect <- c()
n_test <- 0
for (td in unique(small$tag_dist)){
  for(v in unique(small$tag_mut)){
    a <- subset(small, tag_dist == td & tag_mut == v)
    if("para" %in% a$cond){
      n_test <- n_test + 1
      nm <- append(nm, paste0("td",td,"_tm",v))
      md <- append(md, wilcox.test(tag_distance ~ cond, data=a)$p.value)
      effect <- append(effect, cliff.delta(tag_distance ~ cond, data=a)$estimate)

      print("\n*******************")
      print(nm[n_test])
      print(md[n_test])
      print(effect[n_test])
    }
  }
}
print(n_test)


models <- data.frame(names = nm, p_vals=md, deltas=effect)
write.csv(models, "parastab_evolve_paravsmut_deltas.csv", row.names=F)

print(str(small$cond))
