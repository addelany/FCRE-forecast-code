library(tidyverse)

lake_directory <- here::here()
config_set_name <- "glm_aed_flare_v3"
configure_run_file <- "configure_run.yml"
config <- FLAREr::set_up_simulation(configure_run_file,lake_directory, config_set_name = config_set_name)

print('read VERA targets...')

targets_vera <- readr::read_csv("https://amnh1.osn.mghpcc.org/bio230121-bucket01/vera4cast/targets/project_id=vera4cast/duration=P1D/daily-inflow-targets.csv.gz",
                                show_col_types = FALSE)

inflow_hist_dates <- tibble(datetime = seq(min(targets_vera$datetime), max(targets_vera$datetime), by = "1 day"))

variables <- c("datetime", "FLOW", "TEMP", "SALT",
               'OXY_oxy',
               'CAR_dic',
               'CAR_ch4',
               'SIL_rsi',
               'NIT_amm',
               'NIT_nit',
               'PHS_frp',
               'OGM_doc',
               'OGM_docr',
               'OGM_poc',
               'OGM_don',
               'OGM_donr',
               'OGM_pon',
               'OGM_dop',
               'OGM_dopr',
               'OGM_pop',
               'PHY_cyano',
               'PHY_green',
               'PHY_diatom')

df <- targets_vera |>
  filter(!variable %in% c("DN_mgL_sample", "DC_mgL_sample")) |>
  select(datetime, variable, observation) |>
  pivot_wider(names_from = variable, values_from = observation) |>
  right_join(inflow_hist_dates, by = "datetime") |>
  mutate(across(Flow_cms_mean:DIC_mgL_sample, imputeTS::na_interpolation)) |>
  tidyr::fill(Flow_cms_mean:DIC_mgL_sample, .direction = "up") |>
  tidyr::fill(Flow_cms_mean:DIC_mgL_sample, .direction = "down") |>
  dplyr::rename(TEMP = Temp_C_mean,
                FLOW = Flow_cms_mean) |>
  dplyr::mutate(NIT_amm = NH4_ugL_sample*1000*0.001*(1/18.04),
                NIT_nit = NO3NO2_ugL_sample*1000*0.001*(1/62.00), #as all NO2 is converted to NO3
                PHS_frp = SRP_ugL_sample*1000*0.001*(1/94.9714),
                OGM_doc = DOC_mgL_sample*1000*(1/12.01)* 0.10,  #assuming 10% of total DOC is in labile DOC pool (Wetzel page 753)
                OGM_docr = 1.5*DOC_mgL_sample*1000*(1/12.01)* 0.90, #assuming 90% of total DOC is in recalcitrant DOC pool
                TN_ugL = TN_ugL_sample*1000*0.001*(1/14),
                TP_ugL = TP_ugL_sample*1000*0.001*(1/30.97),
                OGM_poc = 0.1*(OGM_doc+OGM_docr), #assuming that 10% of DOC is POC (Wetzel page 755
                OGM_don = (5/6)*(TN_ugL_sample-(NIT_amm+NIT_nit))*0.10, #DON is ~5x greater than PON (Wetzel page 220)
                OGM_donr = (5/6)*(TN_ugL_sample-(NIT_amm+NIT_nit))*0.90, #to keep mass balance with DOC, DONr is 90% of total DON
                OGM_pon = (1/6)*(TN_ugL_sample-(NIT_amm+NIT_nit)), #detemined by subtraction
                OGM_dop = 0.3*(TP_ugL_sample-PHS_frp)*0.10, #Wetzel page 241, 70% of total organic P = particulate organic; 30% = dissolved organic P
                OGM_dopr = 0.3*(TP_ugL_sample-PHS_frp)*0.90,#to keep mass balance with DOC & DON, DOPr is 90% of total DOP
                OGM_pop = TP_ugL_sample-(OGM_dop+OGM_dopr+PHS_frp), # #In lieu of having the adsorbed P pool activated in the model, need to have higher complexed P
                CAR_dic = DIC_mgL_sample*1000*(1/52.515),
                OXY_oxy = rMR::Eq.Ox.conc(TEMP, elevation.m = 506, #creating OXY_oxy column using RMR package, assuming that oxygen is at 100% saturation in this very well-mixed stream
                                          bar.press = NULL, bar.units = NULL,
                                          out.DO.meas = "mg/L",
                                          salinity = 0, salinity.units = "pp.thou"),
                OXY_oxy = OXY_oxy *1000*(1/32),
                CAR_ch4 = CH4_umolL_sample,
                PHY_cyano = 0,
                PHY_green = 0,
                PHY_diatom = 0,
                SIL_rsi = DRSI_mgL_sample*1000*(1/60.08),
                SALT = 0) |>
  dplyr::mutate_if(is.numeric, round, 4) |>
  dplyr::select(dplyr::any_of(variables)) |>
  tidyr::pivot_longer(-c("datetime"), names_to = "variable", values_to = "observation") |>
  dplyr::select(datetime, variable, observation) |>
  dplyr::rename(prediction = observation) |>
  dplyr::mutate(parameter = 1,
                flow_number = 1)

print('finished inflow conversions...')

  if(max(df$datetime) != lubridate::as_date(config$run_config$forecast_start_datetime) - lubridate::days(1)){

    variables <- unique(df$variable)
    full_time <- seq(min(df$datetime), lubridate::as_datetime(config$run_config$forecast_start_datetime), by = "1 day")

    full_data <- list()

    for(i in 1:length(variables)){

      new_data <- tibble(datetime = full_time,
                         variable = rep(variables[i], length(full_time)),
                         parameter = 1,
                         flow_number = 1)

      full_data <- bind_rows(full_data, new_data)

    }

    df <- df |>
      dplyr::right_join(full_data, by = join_by(datetime, variable, parameter, flow_number)) |>
      dplyr::arrange(variable, parameter, flow_number, datetime) |>
      dplyr::group_by(parameter, flow_number, variable) |>
      tidyr::fill(prediction, .direction = "down") |>
      dplyr::ungroup()

  }

print('saving inflow forecasts for FLARE...')

arrow::write_dataset(df, path = file.path(lake_directory, "drivers/inflow/historical/model_id=historical_interp_inflow/site_id=fcre"))

  df |>
    dplyr::filter(variable %in% c("TEMP", "FLOW")) |>
    arrow::write_dataset(path = file.path(lake_directory, "drivers/inflow/historical/model_id=historical_interp_outflow/site_id=fcre"))

inflow_forecast_dir <- "inflow"
print(inflow_forecast_dir)

convert_vera4cast_inflow(reference_date = lubridate::as_date(config$run_config$forecast_start_datetime),
                         model_id = "inflow_gefsClimAED",
                         save_path = file.path(lake_directory, "drivers", inflow_forecast_dir, "future"))
print('inflows converted...')
