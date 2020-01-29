## Full credit to Yang Liu who posted this code online at:
## https://liuyanguu.github.io/post/2019/04/17/ggplot-heatmap-us-50-states-map-and-china-province-map/

# China -------------------------------------------------------------------
get_map_shape <- function(local_file_dir="~/Documents/case_to_infection/data/map_data/"){
  china_map <- rgdal::readOGR(paste0(local_file_dir, "bou2_4p.shp"))
  # extract province information from shap file
  china_province = setDT(china_map@data)
  setnames(china_province, "NAME", "province")
  
  # transform to UTF-8 coding format
  china_province[, province:=iconv(province, from = "GBK", to = "UTF-8")] 
  # create id to join province back to lat and long, id = 0 ~ 924
  china_province[, id:= .I-1] 
  # there are more shapes for one province due to small islands
  china_province[, table(province)]
  
  china_province[, province:= as.factor(province)]
  
  dt_china = setDT(fortify(china_map))
  dt_china[, id:= as.numeric(id)]
  
  setkey(china_province, id); setkey(dt_china, id)
  dt_china <- china_province[dt_china]
  
  # make the province EN, CH label file
  province_CH <- china_province[, levels(province)] # the CH are in UTF-8 code
  province_EN <- c("Shanghai", "Yunnan", "Inner Mongolia", "Beijing", "Taiwan",
                   "Jilin", "Sichuan", "Tianjin City", "Ningxia", "Anhui",
                   "Shandong", "Shanxi", "Guangdong", "Guangxi ", "Xinjiang",
                   "Jiangsu", "Jiangxi", "Hebei", "Henan", "Zhejiang",
                   "Hainan", "Hubei", "Hunan", "Gansu", "Fujian",
                   "Tibet", "Guizhou", "Liaoning", "Chongqing", "Shaanxi",
                   "Qinghai", "Hong Kong", "Heilongjiang"
  )
  
  # some population data (from years ago too)
  value <- c(8893483, 12695396,  8470472,  7355291, 23193638,  9162183, 26383458,  3963604,  1945064, 19322432, 30794664, 10654162, 32222752, 13467663,  6902850, 25635291, 11847841, 20813492, 26404973, 20060115, 2451819, 17253385, 19029894,  7113833, 11971873,   689521, 10745630, 15334912, 10272559, 11084516, 1586635,  7026400, 13192935)
  input_data <- data.table(province_CH, province_EN, value)
  setkey(input_data, province_CH)
  setkey(dt_china, province)
  
  # remove small islands on the South China Sea
  china_map_pop <- input_data[dt_china[AREA>0.1,]]
  
  # create label file of province `label_dt`
  label_dt <- china_map_pop[, .(x = mean(range(long)), y = mean(range(lat)), province_EN, province_CH), by = id]
  label_dt <- unique(label_dt)
  setkey(label_dt, province_EN)
  # I have fine-tuned the label position of some provinces
  label_dt['Inner Mongolia', `:=` (x = 110, y = 42)]
  label_dt['Gansu', `:=` (x = 96.3, y = 40)]
  label_dt['Hebei', `:=` (x = 115.5, y = 38.5)]
  label_dt['Liaoning', `:=` (x = 123, y = 41.5)]
  
  
  p <- ggplot(china_map_pop, aes(x = long, y = lat, group = group, fill=value)) +
    labs(fill = "Population (outdated)")+
    geom_polygon()+
    geom_path()+
    scale_fill_gradientn(colours=rev(heat.colors(10)),na.value="grey90",
                         guide = guide_colourbar(barwidth = 0.8, barheight = 10)) + 
    geom_text(data = label_dt, aes(x=x, y=y, label = province_EN),inherit.aes = F)
  return(list(plot=p,data=china_map_pop,labels=label_dt))
}




###########################
## PLOT STUFF ON A MAP
###########################
## Want to work with cumulative infections
dat_cumu <- province_data  %>% group_by(date, var, province) %>% mutate(cumu_n=cumsum(n))
dat_cumu_quantiles <- dat_cumu %>% group_by(date, var, province) %>% 
  do(data.frame(t(quantile(.$cumu_n, probs = c(0.025,0.5,0.975),na.rm=TRUE)))) %>% ungroup()

dat_cumu_mean <- province_data  %>% group_by(date, var, province) %>% 
  mutate(cumu_n=cumsum(n)) %>% summarise(mean_n = mean(cumu_n))

dat_cumu_quantiles$var <- variable_key2[dat_cumu_quantiles$var]
dat_cumu_mean$var <- variable_key2[dat_cumu_mean$var]
colnames(dat_cumu_quantiles) <- c("date","Variable","province","lower","median","upper")
colnames(dat_cumu_mean) <- c("date","Variable","province","mean")

total_confirmed_prov <- confirm_dat_province %>% group_by(province) %>% summarise(n=sum(n))
total_confirmed_prov <- total_confirmed_prov[order(-total_confirmed_prov$n),]
factor_order <- as.character(total_confirmed_prov$province)

confirm_dat_province$province <- factor(as.character(confirm_dat_province$province), 
                                        levels=factor_order)
dat_cumu_quantiles$province <- factor(as.character(dat_cumu_quantiles$province), 
                                      levels=factor_order)
dat_cumu_mean$province <- factor(as.character(dat_cumu_mean$province), 
                                 levels=factor_order)





map_stuff <- get_map_shape()
china_map_dat <- map_stuff$data
label_dt <- map_stuff$labels

## Do the provinces match up?
in_map_prov <- unique(china_map_dat$province_EN)
in_dat_prov <- as.character(unique(sim_data_quantiles_province$province))
setdiff(in_map_prov, in_dat_prov)
## Missing non-mainland territories and two are misnamed
china_map_dat[china_map_dat$province_EN == "Tianjin City", "province_EN"] <- "Tianjin"
china_map_dat[china_map_dat$province_EN == "Guangxi ", "province_EN"] <- "Guangxi"

missing_provinces <- c("Taiwan", "Hong Kong","Tibet")
colnames(china_map_dat)[2] <- "province"
china_map_dat <- as_tibble(china_map_dat)
china_map_dat <- china_map_dat %>%  mutate(province=ifelse(id==177,"Tibet",province))
china_map_dat <- china_map_dat %>% mutate(province=ifelse(id== 538,"Taiwan",province))
china_map_dat <- china_map_dat %>% mutate(province=ifelse(id== 538,"Hong Kong",province))
## Hong Kong and Taiwan missing so will add as points

in_map_prov <- unique(china_map_dat$province)
in_dat_prov <- as.character(unique(dat_cumu_quantiles$province))
setdiff(in_map_prov, in_dat_prov)

china_map_dat$province <- factor(china_map_dat$province, levels=c(levels(dat_cumu_quantiles$province),"Tibet"))

#dat_cumu_quantiles$province <- as.character(dat_cumu_quantiles$province)
#dat_cumu_quantiles <- dat_cumu_quantiles %>% ungroup() %>% mutate(province=ifelse(is.na(province),"Tibet",province))
#dat_cumu_quantiles$province <- factor(dat_cumu_quantiles$province, levels=c(levels(china_map_dat$province)))


dat_cumu_mean$province <- as.character(dat_cumu_mean$province)
dat_cumu_mean <- dat_cumu_mean %>% ungroup() %>% mutate(ifelse(is.na(province),"Tibet",province))
dat_cumu_mean$province <- factor(dat_cumu_mean$province, levels=c(levels(china_map_dat$province)))

#dat_cumu_quantiles_fill <- dat_cumu_quantiles %>% ungroup() %>% 
#  select(date, Variable, province, median) %>%
#  complete(date, Variable, province, fill=list(median=0,lower=0,upper=0))

dat_cumu_mean_fill <- dat_cumu_mean %>% ungroup() %>% 
  select(date, Variable, province, mean) %>%
  complete(date, Variable, province, fill=list(mean=0))

all_map_dat <- full_join(china_map_dat, dat_cumu_mean_fill, by="province")
all_map_dat <- all_map_dat %>% drop_na()

label_dt[label_dt$id == "898","province_EN"] <- "Hong Kong"

dates <- unique(all_map_dat$date)
dates <- dates[28:length(dates)]
for(i in seq_along(dates)){
  print(i)
  p <- ggplot(all_map_dat[all_map_dat$date == dates[i],], aes(x = long, y = lat, group = group,fill=log10(mean))) +
    geom_polygon()+
    geom_path()+
    scale_fill_gradient2(low="blue",high="red",na.value="blue",mid="orange",limits=c(0,5),midpoint=3,#midpoint=120,limits=c(0,170),
                         guide = guide_colourbar(barwidth = 0.8, barheight = 10)) + 
    geom_text(data = label_dt, aes(x=x, y=y, label = province_EN),inherit.aes = F, col="white") +
    ggtitle(paste0("Date: ", dates[i])) +
    theme_classic() + theme(panel.background = element_rect(fill="black"), legend.position=c(0.95,0.3)) + facet_wrap(~Variable,ncol=1)
  #png(paste0("~/tmp/", i,"_infections.png"),height=10,width=15,res=300,units="in")
  png(paste0("plots/",i,"_infections.png"),height=20,width=15,units="in",res=90)
  plot(p)
  dev.off()
}

