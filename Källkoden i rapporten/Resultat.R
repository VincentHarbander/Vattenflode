source("BM utan trender.R")
source("Block maxima med trend i my.R")
source("Block maxima med trend i phi.R")
source("PoT utan trend.R")
source("PoT med trend i phi.R")

BH_function <- function(data, column, k){
  BH <- list((data["name"]), p.adjust(t(data[column]), method = "BH"))
  names(BH) <- c("name", k)
  data = merge(data, data.frame(BH), by="name")
  return(data)
} 

tester <- function(data, column, replacer, k){
  c <-c()
  for (i in 1:length(data[,column])){
    if (is.nan(data[column][i,1])){
      data[column][i,1] <- 1 
    }
    if (data[column][i,1] < 0.05){
      c <-c(c, data[replacer][i,1])
    }
    else{
      c <- c(c, 0)
    }
  }
  ans <- list(data["name"], c)
  names(ans) <- c("name", k)
  return(merge(data, data.frame(ans), by = "name"))
}

list <- list.files("data")

list2 <- read.csv("Data reglerat.csv", sep = ";")
list <- c()
p <- 0
for (i in 0:53){
  p <- p+1
  print(t(list2["S"])[p])
  if (t(list2["S"])[p] == "o"){
    list <- c(list, t(list2["name"])[p])
  }
}

df = NULL
# Loopar igenom alla filer
for (i in list) {
  df = rbind(df, data.frame(BM_function(i)))
  
}

df = BH_function(df, "ad_bm0", "adbh_bm0")
df = BH_function(df, "xi_p_bm0", "xi_bh_bm0")

df2 = NULL

for (i in list) {
  df2 = rbind(df2, data.frame(BM_my_function(i)))
}


#df2 = BH_function(df2, "ad_bm1", "adbh_bm1")

df2 = BH_function(df2, "loc1_p_bm1", "loc1_bh_bm1")

df2 = BH_function(df2, "xi_p_bm1", "xi_bh_bm1")

df3 = NULL

for (i in list) {
  df3 = rbind(df3, data.frame(BM_phi_function(i)))
}


#df3 = BH_function(df3, "ad_bm2", "adbh_bm2")

df3 = BH_function(df3, "phi1_p_bm2", "phi1_bh_bm2")

df3 = BH_function(df3, "xi_p_bm2", "xi_bh_bm2")

#rejct

df2 = tester(df2, "loc1_bh_bm1","loc1_bm1","loc1_tested_bm1")
df3 = tester(df3, "phi1_bh_bm2","phi1_bm2","phi1_tested_bm2")

df4 = NULL

for (i in list) {
  df4 = rbind(df4, data.frame(PoT_function(i)))
}


df4 = BH_function(df4, "ad_pot0", "adbh_pot0")
df4 = BH_function(df4, "xi_p_pot0", "xi_bh_pot0")
df4 = BH_function(df4, "lambda1_p", "lambda1_bh")

df5 = NULL

for (i in list) {
  df5 = rbind(df5, data.frame(PoT_phi_function(i)))
}

#df5 = rbind(df5, data.frame(PoT_phi_function(list[4])))

#df5 = BH_function(df5, "ad_pot2", "adbh_pot2")
df5 = BH_function(df5, "phi1_p_pot2", "phi1_bh_pot2")
df5 = BH_function(df5, "xi_p_pot2", "xi_bh_pot2")


df_merge <- merge(df,df2,by="name")
df_merge <- merge(df_merge,df3,by="name")
df_merge <- merge(df_merge,df4,by="name")
df_merge <- merge(df_merge,df5,by="name")

# Skriver till en csv
write.csv2(df_merge, "All data o.csv", row.names=FALSE)
