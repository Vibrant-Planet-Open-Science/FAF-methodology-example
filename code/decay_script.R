
samp <- 
  samples |>
  dplyr::group_by(unique_id, area) |>
  dplyr::arrange(unique_id, obs_year, area) |>
  dplyr::mutate(
    n = seq(1,dplyr::n(),1),
    c_change = carbon - dplyr::lag(carbon),
    emission_factor = NA,
    emission_factor = dplyr::case_when(
      cbi > 0.01 & cbi <= 1 ~ 0.3,
      cbi > 1 & cbi <= 2 ~ 0.4,
      cbi > 2  ~ 0.5
    ),
    fire_emissions = dplyr::lead(c_change) * emission_factor,
    fire_emissions = dplyr::if_else(is.na(fire_emissions), 0, fire_emissions),
    dw_contrib = 0,
    dw_contrib = dplyr::if_else( c_change < 0 & is.na(dplyr::lag(cbi)), c_change, dw_contrib),
    dw_contrib = dplyr::if_else( c_change < 0 & !is.na(dplyr::lag(cbi)), c_change - dplyr::lag(fire_emissions), dw_contrib),
    dw_contrib = dplyr::if_else(obs_year == 2015, 0, dw_contrib),
    decay1 = dplyr::lag(dw_contrib) * 0.915,
    decay2 = dplyr::lag(decay1) * 0.915,
    decay3 = dplyr::lag(decay2) * 0.915,
    decay4 = dplyr::lag(decay3) * 0.915,
    decay5 = dplyr::lag(decay4) * 0.915,
    decay6 = dplyr::lag(decay5) * 0.915,
    decay7 = dplyr::lag(decay6) * 0.915,
    decay8 = dplyr::lag(decay7) * 0.915,
    decay9 = dplyr::lag(decay8) * 0.915,
    decay1 = dplyr::if_else(is.na(decay1),0,decay1),
    decay2 = dplyr::if_else(is.na(decay2),0,decay2),
    decay3 = dplyr::if_else(is.na(decay3),0,decay3),
    decay4 = dplyr::if_else(is.na(decay4),0,decay4),
    decay5 = dplyr::if_else(is.na(decay5),0,decay5),
    decay6 = dplyr::if_else(is.na(decay6),0,decay6),
    decay7 = dplyr::if_else(is.na(decay7),0,decay7),
    decay8 = dplyr::if_else(is.na(decay8),0,decay8),
    decay9 = dplyr::if_else(is.na(decay9),0,decay9),
    decay_pool = decay1 + decay2 + decay3 + decay4 + decay5 + decay6 + decay7 + decay8,
    dw_pool = dw_contrib + decay_pool,
    cum_dw_pool = cumsum(dw_contrib),
    emissions_pool = (cum_dw_pool - dw_pool)# + fire_emissions
    # n2 = which(cum_deadwood_pool < 0)[1],
    # decay_pool =  cum_deadwood_pool - deadwood_pool,#decay rate is 0.085, so post decay is 0.915
    # dw_decay_pool = 0.085 * decay_pool,
    # remaining_dw_pool = 0.915 * decay_pool,
    # remaining_dw_pool =  dplyr::if_else(deadwood_pool < 0 & decay_pool == 0, deadwood_pool, remaining_dw_pool),
  ) |> 
  dplyr::select(-c(decay1,decay2,decay3,decay4,decay5,decay6,decay7,decay8,decay9, cum_dw_pool)) |>
  #dplyr::filter(unique_id %in% 434817) |> 
  data.frame()

#Do our summed, pixel levels losses equal our contributions to the dw pool,  emissions, and fire?
samp[samp$c_change < 0,]$c_change |> sum(na.rm = T) == samp$dw_pool[9] + samp$emissions_pool[9] - 12.8
samp[samp$c_change < 0,]$c_change
samp$fire_emissions = c(0,0,0,0,rep(-12.8,5))
samp$fire_emissions = c(0,0,0,0,-12.8,0,0,0,0)
samp$neg_change <- dplyr::if_else(samp$c_change < 0, samp$c_change, 0)
samp$neg_change <- dplyr::if_else(is.na(samp$neg_change),0,samp$neg_change)
samp$cum_loss = cumsum(samp$neg_change)

samp |> tidyr::pivot_longer(cols = c('emissions_pool','dw_pool','fire_emissions')) |>
  ggplot() +
  geom_bar(aes(x = obs_year, y = value, fill = name), stat = 'identity')+
  geom_line(aes(x = obs_year, y = cum_loss), size = 2) +
  theme_bw() + 
  tidyquant::scale_fill_tq(theme = 'dark') +
  ylab('Lost Carbon')

#live wood graph / deadwood graph, both as time series
#fire vs deadwood pool decay wildfire emissions pulse
lm(carbon ~ obs_year,)
colnames(samp)
ggplot() + geom_point(data = samp, aes (x = obs_year, y = carbon, shape = area, color = area), alpha= .4) +
  geom_point(data = samp, aes (x = obs_year, y = -dw_pool, shape = area, color = area), alpha = .4) + theme_bw() +
  # facet_grid(~area) + 
  geom_smooth(data = samp, aes (x = obs_year, y = carbon, linetype = area), color = 'darkgreen', method = 'loess') +
  geom_smooth(data = samp, aes (x = obs_year, y = -dw_pool, linetype = area), color = 'darkred',method = 'loess')

#deadwood pool vs livewood pool
ggplot() + geom_point(data = samp, aes (x = obs_year, y = carbon), color = 'green') +
  geom_point(data = samp, aes (x = obs_year, y = -dw_pool), color = 'red')


samp |> tidyr::pivot_longer(cols = c('emissions_pool','fire_emissions')) |>
  ggplot() +
  geom_bar(aes(x = obs_year, y = value, fill = name), stat = 'identity')+
  #  geom_line(aes(x = obs_year, y = cum_loss), size = 2) +
  theme_bw() + 
  tidyquant::scale_fill_tq(theme = 'dark') +
  ylab('Lost Carbon') + 
  geom_path(data = samp, aes (x = obs_year, y = carbon), color = 'darkgreen') +
  geom_path(data = samp, aes (x = obs_year, y = -dw_pool), color = 'red')


