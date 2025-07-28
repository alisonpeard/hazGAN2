fit_tails <- function(
    footprint, train, distn,
    two_tailed = FALSE, grid_i = NULL
) {
  
  # Fit upper tail with error handling
  upper_result <- tryCatch({
    log_info(paste0("Fitting upper tail for grid cell ", grid_i))
    fit <- distn$threshold_selector(train$variable)
    log_info(paste0("Upper tail parameters: ", fit$params))
    
    # append "_upper" suffix to params
    params <- setNames(fit$params, paste0(names(fit$params), "_upper"))
    log_info(paste(names(params), collapse = ", "))
    
    list(
      params = params,
      p_value = fit$p.value,
      pk = fit$pk,
      success = TRUE
    )
  }, error = function(e) {
    log_error(paste0("Upper tail fitting failed for grid cell ", grid_i, ": ", e$message))
    list(
      params = list(thresh_upper = Inf, scale_upper = NA, shape_upper = NA),
      p_value = 0,
      pk = 0,
      success = FALSE
    )
  })
  
  # Add upper tail results to footprint
  footprint <- footprint |> mutate(!!!upper_result$params)
  footprint$p_upper <- upper_result$p_value
  footprint$pk_upper <- upper_result$pk
  
  # Fit lower tail if requested
  if (two_tailed) {
    lower_result <- tryCatch({
      log_info(paste0("Fitting lower tail for grid cell ", grid_i))
      fit_lower <- distn$threshold_selector(-train$variable)
      log_info(paste0("Lower tail parameters: ", fit_lower$params))
      
      # append "_lower" suffix to params
      params_lower <- setNames(fit_lower$params, paste0(names(fit_lower$params), "_lower"))
      
      list(
        params = params_lower,
        p_value = fit_lower$p.value,
        pk = fit_lower$pk,
        success = TRUE
      )
    }, error = function(e) {
      log_error(paste0("Lower tail fitting failed for grid cell ", grid_i, ": ", e$message))
      list(
        params = list(thresh_lower = Inf, scale_lower = NA, shape_lower = NA),
        p_value = 0,
        pk = 0,
        success = FALSE
      )
    })
    
    # Add lower tail results to footprint
    footprint <- footprint |> mutate(!!!lower_result$params)
    footprint$p_lower <- lower_result$p_value
    footprint$pk_lower <- lower_result$pk
  } else {
    # Set lower tail params to NA if not two-tailed
    footprint$thresh_lower <- NA
    footprint$scale_lower <- NA  
    footprint$shape_lower <- NA
    footprint$p_lower <- NA
    footprint$pk_lower <- NA
  }
  
  # Transform variable to uniform (only if upper tail fitting succeeded)
  if (upper_result$success) {
    footprint$scdf <- scdf(
      train$variable, upper_result$params, cdf = distn$cdf, two_tailed = two_tailed
    )(footprint$variable)
  } else {
    footprint$scdf <- ecdf(train$variable)(footprint$variable)
  }
  
  footprint$ecdf <- ecdf(train$variable)(footprint$variable)
  
  return(footprint)
}
