
df <- '/Users/eyackulic/Desktop/final_cms_data_2025-03-31_111.rds' |> 
  readRDS() |> 
  dplyr::mutate(
    designation = dplyr::case_when(
      dataset %in% 'real' & area %in% 'project' ~ 'treatment',
      dataset %in% 'fake' & area %in% 'project' ~ 'placebo',
      dataset %in% 'real' & area %in% 'ref_region' ~ 'ref region',
      dataset %in% 'fake' & area %in% 'ref_region' ~ 'placebo ref region'
    ),
    designation = factor(designation, levels = c( 'treatment','ref region','placebo ref region', 'placebo')),
    composite_2015_median = composite_2015_median/2, #change biomass to carbon for all years (half of biomass is carbon)
    composite_2016_median = composite_2016_median/2,
    composite_2017_median = composite_2017_median/2,
    composite_2018_median = composite_2018_median/2,
    composite_2019_median = composite_2019_median/2,
    composite_2020_median = composite_2020_median/2,
    composite_2021_median = composite_2021_median/2,
    composite_2022_median = composite_2022_median/2,
    composite_2023_median = composite_2023_median/2,
    max_rdnbr_year = max_rdnbr_year - 1,
    fire_change = NA,
    fire_change = dplyr::case_when(
      max_rdnbr_year %in% 2015 ~ NA,
      max_rdnbr_year %in% 2016 ~ (composite_2017_median - composite_2015_median) / composite_2015_median,
      max_rdnbr_year %in% 2017 ~ (composite_2018_median - composite_2016_median) / composite_2016_median,
      max_rdnbr_year %in% 2018 ~ (composite_2019_median - composite_2017_median) / composite_2017_median,
      max_rdnbr_year %in% 2019 ~ (composite_2020_median - composite_2018_median) / composite_2018_median,
      max_rdnbr_year %in% 2020 ~ (composite_2021_median - composite_2019_median) / composite_2019_median,
      max_rdnbr_year %in% 2021 ~ (composite_2022_median - composite_2020_median) / composite_2020_median,
      max_rdnbr_year %in% 2022 ~ (composite_2023_median - composite_2021_median) / composite_2021_median
    ),    
    pre_fire = NA,
    pre_fire = dplyr::case_when(
      max_rdnbr_year %in% 2015 ~ NA,
      max_rdnbr_year %in% 2016 ~ composite_2015_median,
      max_rdnbr_year %in% 2017 ~ composite_2016_median,
      max_rdnbr_year %in% 2018 ~ composite_2017_median,
      max_rdnbr_year %in% 2019 ~  composite_2018_median,
      max_rdnbr_year %in% 2020 ~  composite_2019_median,
      max_rdnbr_year %in% 2021 ~ composite_2020_median,
      max_rdnbr_year %in% 2022 ~  composite_2021_median
    )
  )


df |> 
  dplyr::filter(dataset %in% 'real') |>
  dplyr::group_by(
    distid, area
  ) |>
  dplyr::reframe(
    m = mean((composite_2023_median - composite_2015_median) / composite_2015_median),
    max = max(max_rdnbr, na.rm = T)
  ) |> 
  dplyr::filter(max > 0) |> #, area %in% 'project') |>
  dplyr::arrange(
    -m
  ) |>
  dplyr::filter(distid %in% c(47509))


test <- 
  df |> 
  dplyr::filter(distid %in% 47509, dataset %in% 'real')

test |>
  dplyr::group_by(area) |> 
  dplyr::reframe(
    n = dplyr::n(),
    m2015 = mean(composite_2015_median),
    m2016 = mean(composite_2016_median),
    m2017 = mean(composite_2017_median),
    m2018 = mean(composite_2018_median),
    m2019 = mean(composite_2019_median),
    m2020 = mean(composite_2020_median),
    m2021 = mean(composite_2021_median),
    m2022 = mean(composite_2022_median),
    m2023 = mean(composite_2023_median)
  )


means <-
  test |>
  dplyr::group_by(area) |>
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean, na.rm = TRUE))

library(ggplot2)
means |>
  tidyr::pivot_longer(dplyr::starts_with('composite')) |>
  dplyr::mutate(year = stringr::str_sub(name,11,14)) |>
ggplot(aes(y = value, x = year, color = area, group = area)) + geom_point() + geom_line()

means |>
  tidyr::pivot_longer(dplyr::starts_with('composite')) |> 
  dplyr::mutate(year = stringr::str_sub(name,11,14)) |>
  write.csv('/Users/eyackulic/Downloads/test_data.csv')

colnames(test)
test_long <-
  test |> 
  dplyr::select(-c(dataset, designation, pre_fire, fire_change, lat, lon, tolerance, disttype, trt_yr, ecoregion)) |>
  tidyr::pivot_wider(names_from = max_rdnbr_year, values_from = max_rdnbr, names_prefix = 'fireyear_') |>
  tidyr::pivot_longer(dplyr::starts_with(c('composite','fireyear'))) |>
  dplyr::group_by(area, name) |>
  dplyr::mutate(
#    vals = mean(unlist(value),na.rm = T),
#    n = dplyr::n(),
    variable = stringr::str_sub(name, 1,4),
    variable_name = 
      dplyr::case_when(
        variable %in% 'fire' ~ 'cbi',
        variable %in% 'comp' ~ 'carbon'
      ),
    obs_year = dplyr::case_when(
      variable %in% 'comp' ~ stringr::str_sub(name, 11,14),
      variable %in% 'fire' ~ stringr::str_sub(name, 10,13)
    ),
    obs_year = dplyr::if_else(obs_year == -1, NA, as.numeric(obs_year))
    ) |>
  dplyr::select(-variable) |>
  dplyr::rename(unique_id = Location) 

carbon <- test_long |> dplyr::filter(variable_name %in% 'carbon')
fire <- test_long |> dplyr::filter(variable_name %in% 'cbi', !is.na(obs_year))


full_dat <- 
  carbon |> 
  dplyr::left_join(fire, by = c('unique_id','area','distid','obs_year')) |>
  dplyr::select(unique_id, area, distid, obs_year, value.x, value.y) |>
  dplyr::rename(
    carbon = value.x,
    cbi = value.y
  )

ggplot(full_dat, aes(x = factor(obs_year), fill = area, y = carbon)) + geom_boxplot()
ggplot(full_dat, aes(x = factor(obs_year), fill = area, y = cbi)) + geom_boxplot()

rand_proj_locations = 
  full_dat[which(full_dat$area %in% 'project'),]$unique_id |>
  unique() |> 
  data.frame() |> 
  dplyr::sample_n(20)|>
  dplyr::rename(id = 1)

rand_ref_locations <- 
  full_dat[which(full_dat$area %in% 'ref_region'),]$unique_id |> 
  unique() |> 
  data.frame() |> 
  dplyr::sample_n(20) |>
  dplyr::rename(id = 1)
#filter for unique id and area match and recombine, send to kat

ref_samp <-
  full_dat |>
  dplyr::filter(unique_id %in% rand_ref_locations$id)

proj_samp <-
  full_dat |>
  dplyr::filter(unique_id %in% rand_proj_locations$id)

table(ref_samp$obs_year)
table(proj_samp$obs_year)

samples <- 
  dplyr::bind_rows(
  ref_samp, 
  proj_samp
)

samples[!is.na(samples$cbi),]
write.csv(samples, '/Users/eyackulic/Downloads/proj47509_sample.csv')
