#
options(timeout = 400000)
source("~/Documents/GitHub/Fire-Adapted-Forests/code/ref_region.R")

carbon <- '/Users/eyackulic/workspace/Tuolumne_Layers' |> list.files(full.names = T)
ca_rast <- carbon[29:33]

projection = get_projection(ca_rast)
#building reference regions for an example shapefile 
# tuolumne <- '/Users/eyackulic/Downloads/tl_2024_us_county/tl_2024_us_county.shp' |> sf::read_sf() |>
#   dplyr::filter(NAME == 'Tuolumne') |>
#   sf::write_sf('/Users/eyackulic/Downloads/tl_2024_us_county/tuolumne.gpkg')

tcsi_eco_regions <- prepRR(
  aoi = '/Users/eyackulic/Downloads/phoenix_lake.gpkg',
  projection = projection,
  l4eco =  '/Users/eyackulic/Downloads/us_eco_l4_state_boundaries/us_eco_l4.shp'
) 

evt_path = '/Users/eyackulic/Downloads/LF2016_EVT_200_CONUS/Tif/LC16_EVT_200.tif'
evt <- conform_raster(l4_ref_region = tcsi_eco_regions,raster_path = evt_path,carbon_path = ca_rast)

rs_carbon <- terra::rast(ca_rast) * evt
rs_vals <- terra::values(rs_carbon/2, na.rm = T) |> data.frame() |> magrittr::set_colnames('Value')
colnames(rs_vals) <- rep('Value',5)
rs_vals <- rs_vals |> tidyr::pivot_longer(dplyr::starts_with('Value'),values_to = 'Value')|> dplyr::select(Value)

rs_vals$Source <- 'CMS'
# 
# sample_data <- seq(0,150,by = 1) |> data.frame() 
# colnames(sample_data) <- 'Value'
# biomass <- sample(sample_data$Value, size = 10000, replace = T) |> data.frame()
# biomass$source <- 'CMS'
# colnames(biomass) <- c('Value','Source')

#take in a project location, expand it to l4 ecoregion, read in rFIA data, mask to l4 ecoregion
#crop and mask biomass layer to l4 ecoregion. calculate biomass at every FIA plot (per ha) ,
#mask the biomass layer to forested areas (EVT == unique forest types in project area). 
#run a rope on biomass layer distribution by comparing each FIA plot. 
#crop to l4 ecoregion
#CA <- rFIA::getFIA(state = 'CA', nCores = 10, dir = tempdir())
clip <- rFIA::clipFIA(CA)

sf_obj <- 
  clip$PLOT |>
#  dplyr::filter(is.na(RDDISTCD) | RDDISTCD == 6) |>
  sf::st_as_sf(coords = c('LON','LAT'), crs = 'EPSG:4326') |>
  sf::st_transform(sf::st_crs(projection)) |> 
  sf::st_intersection(y = tcsi_eco_regions) |>
  data.frame()

#filter plots exceed max(bins) * .10 *mean(bins)
carb_vals <- 
  rFIA::carbon(rFIA::clipFIA(CA), byPlot = TRUE,byComponent = T) |>
  dplyr::filter(COMPONENT %in% 'AG_OVER_LIVE') |> #, PROP_FOREST >= .8) |>
  dplyr::filter(POOL %in% 'AG_LIVE') |> #, PROP_FOREST >= .8) |>
  dplyr::group_by(PLT_CN) |> 
  dplyr::reframe(
    TotalCarbon = 2.47105 * sum(CARB_ACRE,na.rm = T),
    YEAR = YEAR) |>
  dplyr::distinct() |>
  dplyr::filter(YEAR >= 2018) |>
  dplyr::left_join(clip$COND[c('PLT_CN','FORTYPCD')], by = 'PLT_CN') |>
  dplyr::filter(FORTYPCD %in% c(221,371),PLT_CN %in% sf_obj$PLT_CN)

ggplot(carb_vals, aes(x = factor(YEAR), y = TotalCarbon)) + geom_boxplot() + geom_jitter()
carb_vals$YEAR |> table()
carb_vals$YEAR |> table() |> sum()

summary(carb_vals$TotalCarbon)
summary(cms_samp$Value)
summary(rs_vals$Value)

carb_data <- dplyr::bind_cols(carb_vals$TotalCarbon, 'FIA')
colnames(carb_data) <- c('Value','Source')
cms_samp <- dplyr::sample_n(rs_vals, 50000)

data <- dplyr::bind_rows(cms_samp, carb_data)
hist(rs_vals$Value)
hist(carb_data$Value)
a <- 
  ggplot() +
  geom_histogram(data = data |> dplyr::filter(Source %in% 'FIA'), aes(x = Value), fill = 'blue') +
  xlim(0,400) + 
  geom_text(aes(x = 300, y = 12, label = 'FIA Data for Tuolumne / Central Sierra (5h)'))
b <- 
  ggplot() +
  geom_histogram(data = data |> dplyr::filter(Source %in% 'CMS'), aes(x = Value), fill = 'red') +
  xlim(0,400) +
  geom_text(aes(x = 300, y = 30000, label = 'CMS Data for Tuolumne / Central Sierra (5h)'))
gridExtra::grid.arrange(a,a2,b,b2)

library(easystats)
library(rstanarm)
library(ggplot2)

b_mod <- stan_glm(Value ~ Source, data = data) #area is either 'project' or 'ref_region'
#b_mod |> saveRDS('/Users/eyackulic/workspace/model_runs/comp_td_model_all.rda') #optional local saving 
describe_posterior(b_mod)
ps <- get_parameters(b_mod)

ggplot(ps) + 
  geom_density(aes(x = `SourceFIA`), color = 'blue')
#density plot of parameters

#rope_value <- 0.1 * sd(med$max_rdnbr)
rope_range <- rope_range(b_mod) # rope range is +/- 10% of the standard deviation ^^

#calculate how many values within the 89th percentile  fall into the rope range
# 0 = highly significant / 100 = not significant at all
rope(ps$SourceFIA, range = rope_range, ci = 0.89) 
ks.test(rs_vals$Value, carb_data$Value)
t.test(cms_samp$Value, carb_data$Value)
equivalence::tost(cms_samp$Value, carb_data$Value,conf.level = 0.90)

#3-4 more test sites - Mogollon Rim, Rocky Mountains (CO), Dryland E. Washington
