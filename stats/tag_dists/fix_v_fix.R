
library(dplyr)
library(ggplot2)
library(viridis)
library(effsize) # stats

print("about to read mut")

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
a$sym_tag[a$sym_tag==""]<-NA
a$cond<-"mut"# <- org_dump_df

print("done with mut; about to read para")
b <- read.csv(paste0(filepath,"inc_syms_only_dump_file.dat"),h=T,
                                    #seed,tag_dist,tag_mut
                        colClasses=c("factor","factor","factor",
                                     #host_int,sym_int
                                    "double","double",
                                    # repro/to/from counts
                                    "double","double","double","double","double","double",
                                    #host_tag,sym_tag,tag_distance
                                    "character","character","double"))
b$sym_tag[b$sym_tag==""]<-NA
b$cond<-"para"# <- org_dump_df
print("done with para; about to rbind and del")
small <- rbind(subset(b, !is.na(sym_tag)), subset(a, !is.na(sym_tag)))
rm(a)
rm(b)

print("done binding and deleting; about to transform")
identifying_cols <- 1:3
for (x in identifying_cols){
  small[[x]]<- as.factor(small[[x]])
}
small <- small %>% group_by(tag_dist,vt) %>% mutate(para_pres = "para"%in%cond)
small <- small %>% group_by(tag_dist,vt) %>% mutate(mut_pres = "mut"%in%cond)
small <- subset(small, para_pres == T & mut_pres == T)
small$cond <- as.factor(small$cond)

print("done transforming; about to calculate stats")

nm <- c()
md <- c()
effect <- c()
n_test <- 0
for (td in unique(small$tag_dist)){
  for(v in unique(small$vt)){
    a <- subset(small, tag_dist == td & vt == v)
    if("para" %in% a$cond){
      n_test <- n_test + 1
      print(paste0("td",td,"_vt",v))
      nm <- append(nm, paste0("td",td,"_vt",v))
      md <- append(md, wilcox.test(tag_distance ~ cond, data=a)$p.value)
      effect <- append(effect, cliff.delta(tag_distance ~ cond, data=a)$estimate)

      print("\n*******************")
      print(nm)
      print(md)
      print(effect)
    }
  }
}
print(n_test)


models <- data.frame(names = nm, p_vals=md, deltas=effect)

write.csv(models, "vtsweep_control.csv", row.names=F)
