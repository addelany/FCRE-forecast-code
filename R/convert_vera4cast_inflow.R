convert_vera4cast_inflow <- function(reference_date, model_id, save_path){

variables <- c("TP_ugL_sample", "NH4_ugL_sample","NO3NO2_ugL_sample",
               "SRP_ugL_sample","DOC_mgL_sample","DRSI_mgL_sample",
               "TN_ugL_sample", "CH4_umolL_sample", "DIC_mgL_sample",
               "Flow_cms_mean", "Temp_C_mean")

forecast_df <- NULL

for(i in 1:length(variables)){

  s3_archive <- arrow::s3_bucket(bucket = glue::glue("bio230121-bucket01/vera4cast/forecasts/archive-parquet/project_id=vera4cast/duration=P1D/variable={variables[i]}/model_id={model_id}"),
                                 endpoint_override = "https://amnh1.osn.mghpcc.org",
                                 anonymous = TRUE)
  avail_dates <- gsub("reference_date=", "", s3_archive$ls())

  if(as.character(reference_date) %in% avail_dates) {
    s3 <- arrow::s3_bucket(bucket = glue::glue("bio230121-bucket01/vera4cast/forecasts/archive-parquet/project_id=vera4cast/duration=P1D/variable={variables[i]}/model_id={model_id}/reference_date={reference_date}"),
                           endpoint_override = "https://amnh1.osn.mghpcc.org",
                           anonymous = TRUE)
  }else{
    # fall back to the non-archive bucket if the date isn't in the archive yet
    s3 <- arrow::s3_bucket(bucket = glue::glue("bio230121-bucket01/vera4cast/forecasts/parquet/project_id=vera4cast/duration=P1D/variable={variables[i]}/model_id={model_id}/reference_date={reference_date}"),
                           endpoint_override = "https://amnh1.osn.mghpcc.org",
                           anonymous = TRUE)
  }

  ## test to see if inflow forecast exists ##
  tryCatch({
    df <- arrow::open_dataset(s3) |>
      dplyr::filter(site_id == "tubr") |>
      dplyr::collect() |>
      dplyr::mutate(variable = variables[i])
  }, error = function(e) {
    stop(paste('\nInflow forecasts were not found at the following path: ',
                glue::glue("\"bio230121-bucket01/vera4cast/forecasts/parquet/project_id=vera4cast/duration=P1D/variable={variables[i]}/model_id={model_id}/reference_date={reference_date}\"\n"),
         'Check that inflow model is submitting all needed variables for reference_date value',
         'Stopping workflow...',
         sep='\n'))

  })

  # df <- arrow::open_dataset(s3) |>
  #   dplyr::filter(site_id == "tubr") |>
  #   dplyr::collect() |>
  #   dplyr::mutate(variable = variables[i])

  forecast_df <- dplyr::bind_rows(forecast_df, df)

}

VARS <- c("parameter", "datetime", "FLOW", "TEMP", "SALT",
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

glm_df_inflow <- forecast_df |>
  dplyr::select(datetime, variable, prediction, parameter) |>
  tidyr::pivot_wider(names_from = variable, values_from = prediction) |>
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
  dplyr::select(dplyr::any_of(VARS)) |>
  tidyr::pivot_longer(-c("datetime","parameter"), names_to = "variable", values_to = "prediction") |>
  dplyr::mutate(datetime = lubridate::as_date(datetime),
                model_id = paste0("inflow-aed"),
                site_id = "fcre",
                family = "ensemble",
                flow_number = 1,
                reference_datetime = lubridate::as_datetime(reference_date),
                reference_date = as.character(reference_date)) |>
  dplyr::select(model_id, site_id, reference_datetime, datetime, family, parameter, variable, prediction, flow_number, reference_date)

arrow::write_dataset(glm_df_inflow, path = save_path, partitioning = c("model_id", "reference_date", "site_id"))


glm_df_outflow <- glm_df_inflow |>
  dplyr::select(datetime, parameter, variable, prediction) |>
  dplyr::filter(variable %in% c("FLOW","TEMP")) |>
  dplyr::mutate(datetime = lubridate::as_date(datetime),
                model_id = paste0("outflow-aed"),
                site_id = "fcre",
                family = "ensemble",
                flow_number = 1,
                reference_datetime = lubridate::as_datetime(reference_date),
                reference_date = as.character(reference_date)) |>
  dplyr::select(model_id, site_id, reference_datetime, datetime, family, parameter, variable, prediction, flow_number, reference_date)


arrow::write_dataset(glm_df_outflow, path = save_path, partitioning = c("model_id", "reference_date", "site_id"))
}
