
library(dplyr)
library(ggplot2)
library(viridis)
library(effsize) # stats

# inc_sym_only_dump_file are org_dump_files with rows sans syms removed 

filepath<- ""
a <- read.csv(paste0(filepath,"inc_syms_only_dump_file.dat"),h=T,
                                    #seed,tag_dist,tag_mut
                        colClasses=c("factor","factor","factor",
                                     #host_int,sym_int
                                    "double","double",
                                    # repro/to/from counts
                                    "double","double","double","double","double","double",
                                    #host_tag,sym_tag,tag_distance
                                    "character","character","double"))
a$cond<-"evolved"# <- org_dump_df
a <- subset(a, sym_int >= 0.6)

print("done with evolved; about to read fixed")
b <- read.csv(paste0(filepath,"inc_syms_only_dump_file.dat"),h=T,
                                    #seed,tag_dist,tag_mut
                        colClasses=c("factor","factor","factor",
                                     #host_int,sym_int
                                    "double","double",
                                    # repro/to/from counts
                                    "double","double","double","double","double","double",
                                    #host_tag,sym_tag,tag_distance
                                    "character","character","double"))
#b$sym_tag[b$sym_tag==""]<-NA
b$cond<-"fixed"# <- org_dump_df
print("done with fixed; about to rbind and del")
small <- rbind(b, a)
rm(a)
rm(b)

print("done binding and deleting; about to transform")
identifying_cols <- 1:3
for (x in identifying_cols){
  small[[x]]<- as.factor(small[[x]])
}
small <- small %>% group_by(tag_dist,vt) %>% mutate(fixed_pres = "fixed"%in%cond)
small <- small %>% group_by(tag_dist,vt) %>% mutate(evolved_pres = "evolved"%in%cond)
small <- subset(small, evolved_pres == T & fixed_pres == T)
small$cond <- as.factor(small$cond)

print("done transforming; about to calculate stats")

nm <- c()
md <- c()
effect <- c()
n_test <- 0
for (td in unique(small$tag_dist)){
  for(v in unique(small$vt)){
    a <- subset(small, tag_dist == td & vt == v)
    if("fixed" %in% a$cond & "evolved"%in% a$cond){
      n_test <- n_test + 1
      nm <- append(nm, paste0("td",td,"_vt",v))
      md <- append(md, wilcox.test(tag_distance ~ cond, data=a)$p.value)

      print("\n*******************")
      print(nm[n_test])
      print(md[n_test])
    }
  }
}
print(n_test)


models <- data.frame(names = nm, p_vals=md)
write.csv(models, "vtsweep_mutualists_pvals.csv", row.names=F)

print("done calculating p vals!")

nm <- c()
md <- c()
effect <- c()
n_test <- 0
for (td in unique(small$tag_dist)){
  for(v in unique(small$vt)){
    a <- subset(small, tag_dist == td & vt == v)
    if("fixed" %in% a$cond & "evolved"%in% a$cond){
      n_test <- n_test + 1
      nm <- append(nm, paste0("td",td,"_vt",v))
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

print(str(small$cond))
models <- data.frame(names = nm, p_vals=md, deltas=effect)

write.csv(models, "vtsweep_mutualists_deltas.csv", row.names=F)
