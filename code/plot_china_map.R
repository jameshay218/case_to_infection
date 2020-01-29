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

