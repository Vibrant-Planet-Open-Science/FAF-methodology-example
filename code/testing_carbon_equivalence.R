#
options(timeout = 400000)

#take in a project location, expand it to l4 ecoregion, read in rFIA data, mask to l4 ecoregion
#crop and mask biomass layer to l4 ecoregion. calculate biomass at every FIA plot (per ha) ,
#mask the biomass layer to forested areas (EVT == unique forest types in project area). 
#run a rope on biomass layer distribution by comparing each FIA plot. 
#crop to l4 ecoregion
CA <- rFIA::getFIA(state = 'VT', nCores = 4, dir = tempdir())

#filter plots exceed max(bins) * .10 *mean(bins)
carb <- rFIA::carbon(rFIA::clipFIA(CA), byPlot = TRUE,byComponent = T) |> dplyr::filter(POOL %in% 'AG_LIVE') |>
  dplyr::group_by(PLT_CN) |> dplyr::reframe(TotalCarbon = 2.47105 * sum(CARB_ACRE,na.rm = T), YEAR = YEAR) |> dplyr::distinct()

sample_data <- seq(0,150,by = 1) |> data.frame() 
colnames(sample_data) <- 'Value'
biomass <- sample(sample_data$Value, size = 10000, replace = T) |> data.frame()
biomass$source <- 'CMS'
colnames(biomass) <- c('Value','Source')

summary(carb$TotalCarbon)
summary(sample_data$Value)

carb_data <- dplyr::bind_cols(carb$TotalCarbon, 'FIA')
colnames(carb_data) <- c('Value','Source')

data <- dplyr::bind_rows(biomass, carb_data)

ggplot(data = sample_data, aes(x = Value)) +geom_histogram()
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
