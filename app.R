# dna shearing protocol builder v2
# built by bella pfeiffer (bpfeiffer@covaris.com)
# available here: https://bellabluebird.shinyapps.io/covaris_protocol_builder/
# last updated: 7/21/25

# notes from 7/17 meeting
# needs for future development (from Julie + Katelyn)
#   - lack of protocol availability
#     - large fragment shearing (400 - 1500 bp)
#     - AFA tube TPX 520291
#       - talk with Vanessa about truCOVER validation data 
#     - pulsing protocols 
#     - high priority instruments: R230, LE220, ML230
#   - ability to map protocols between different devices + consumables, ie:
#   - "we shear to xx on the ME220, how can we do that on the R230"
#   - more consideration for sample volume
#       - 50 in 130 !!
#   - fragment analyzer comparison
#       - planning an experiment to compare between different analyzers 
#   - adding additional variables to protocols  (ie water level, z-height, etc)

# bella - to deploy this again
# 

# load libraries
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}
pacman::p_load(shiny, rsconnect)

# define column mappings 
COLS <- list(
  device = 2,
  vessel = 3,
  volume = 4,  
  bp_mode = 10,
  cycles = 11,
  iterations = 12,
  pip = 13,
  duty_factor = 14,
  duration = 16,
  energy = 17
)

# ui definition
ui <- fluidPage(
  tags$head(
    tags$title("🧬 Covaris Protocol Builder") # tab title for web
  ),
  
  titlePanel("Covaris Protocol Builder"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("datafile", "Choose CSV File", accept = c(".csv")),
      hr(),
      
      # dynamic dropdowns
      selectInput("device", "Select Device:", 
                  choices = list("Upload data first" = "")),
      selectInput("vessel", "Select Vessel:", 
                  choices = list("Select device first" = "")),
      selectInput("volume", "Select Sample Volume:",
                  choices = list("Select vessel first" = "")),
      radioButtons("protocol_type", "Protocol Type:",
                   choices = c("Continuous" = "Continuous",
                               "Pulsing (Majority Not Availab)" = "Pulsing"),
                   selected = "Continuous"),
      
      # target selection
      sliderInput("target_bp", "Target Base Pair (bp):", 
                  min = 100, max = 1500, value = 500, step = 25),
      actionButton("calculate", "Calculate Energy", class = "btn-primary"),
      
      hr(),
      # simplified status display
      h4("Model Status"),
      verbatimTextOutput("model_summary")
    ),
    
    mainPanel(
      tabsetPanel(
        # main analysis tab
        tabPanel("Analysis",
                 plotOutput("energy_plot", height = "600px"),
                 
                 wellPanel(
                   h4("Prediction"),
                   verbatimTextOutput("prediction"),
                   hr(),
                   h4("Protocol"),
                   verbatimTextOutput("protocol"),
                   br(),
                   downloadButton("download_protocol", "Download Protocol", 
                                  class = "btn-success", style = "width: 100%;")
                 )
          ),
        
        # equations overview tab
        tabPanel("All Models",
                 h4("Model Summary"),
                 verbatimTextOutput("summary_stats"),
                 hr(),
                 tableOutput("equations_table"),
                 downloadButton("download_equations", "Download Equations"),
                 hr(),
                 h4("Example Model Curves"),
                 plotOutput("example_curves", height = "600px"))
      )
    )
  )
)

# server logic
server <- function(input, output, session) {
  
  # reactive values to store state
  values <- reactiveValues(
    data = NULL,
    current_model = NULL,
    all_models = NULL
  )
  
  # helper functions
  
  # validate uploaded csv
  validate_csv <- function(df) {
    if(ncol(df) < 17) {
      return(list(valid = FALSE, message = "CSV must have at least 17 columns"))
    }
    
    # check numeric columns
    numeric_cols <- c(COLS$volume, COLS$bp_mode, COLS$pip, COLS$duty_factor, COLS$duration, COLS$energy)
    for(col in numeric_cols) {
      if(!all(is.na(df[, col]) | !is.na(as.numeric(df[, col])))) {
        return(list(valid = FALSE, message = paste("Column", col, "must be numeric")))
      }
    }
    
    return(list(valid = TRUE))
  }
  
  # add protocol type to data
  add_protocol_type <- function(df) {
    df$protocol_type <- ifelse(!is.na(df[, COLS$cycles]) & !is.na(df[, COLS$iterations]), 
                               "Pulsing", "Continuous")
    df
  }
  
  # filter data by selections 
  filter_data <- function(df, device, vessel, volume, protocol = NULL) {
    filtered <- df[df[, COLS$device] == device & 
                     df[, COLS$vessel] == vessel & 
                     df[, COLS$volume] == volume, ]
    
    if(!is.null(protocol)) {
      filtered <- filtered[filtered$protocol_type == protocol, ]
    }
    
    # extract and clean data
    bp <- as.numeric(filtered[, COLS$bp_mode])
    energy <- as.numeric(filtered[, COLS$energy])
    
    # remove invalid values
    valid <- !is.na(bp) & !is.na(energy) & energy > 0 & bp > 0
    
    list(
      raw = filtered,
      bp = bp[valid],
      energy = energy[valid],
      n_points = sum(valid),
      n_pulsing = sum(filtered$protocol_type == "Pulsing" & valid),
      n_continuous = sum(filtered$protocol_type == "Continuous" & valid)
    )
  }
  
  # fit best model to data
  fit_best_model <- function(bp, energy) {
    if(length(bp) < 3) {
      return(list(success = FALSE, message = "Need at least 3 data points"))
    }
    
    df_model <- data.frame(bp_mode = bp, energy_j = energy)
    
    # fit candidate models
    models <- list()
    tryCatch({
      models$power <- lm(log(energy_j) ~ log(bp_mode), data = df_model)
    }, error = function(e) {
      return(list(success = FALSE, message = e$message))
    })
    
    # select best by aic
    aics <- sapply(models, function(m) tryCatch(AIC(m), error = function(e) Inf))
    best_name <- names(which.min(aics))
    best_model <- models[[best_name]]
    
    # add metadata
    attr(best_model, "model_type") <- best_name
    attr(best_model, "data_range") <- range(bp)
    attr(best_model, "bp_data") <- bp
    attr(best_model, "energy_data") <- energy
    
    list(success = TRUE, model = best_model, name = best_name)
  }
  
  # predict with intervals
  predict_energy <- function(model, bp_target) {
    model_type <- attr(model, "model_type")
    
    # create appropriate prediction dataframe
    if(model_type == "power") {
      pred_df <- data.frame(bp_mode = bp_target)
    } else {
      pred_df <- data.frame(bp_mode = bp_target)
    }
    
    # get predictions with both confidence and prediction intervals
    log_conf <- predict(model, newdata = pred_df, interval = "confidence", level = 0.95)
    log_pred <- predict(model, newdata = pred_df, interval = "prediction", level = 0.95)
    
    # back-transform
    list(
      fit = exp(log_conf[, "fit"]),
      conf_lower = exp(log_conf[, "lwr"]),
      conf_upper = exp(log_conf[, "upr"]),
      pred_lower = exp(log_pred[, "lwr"]),
      pred_upper = exp(log_pred[, "upr"]),
      extrapolating = bp_target < min(attr(model, "data_range")) | 
        bp_target > max(attr(model, "data_range"))
    )
  }
  
  # format equation string
  format_equation <- function(model) {
    coefs <- coef(model)
    model_type <- attr(model, "model_type")
    
    switch(model_type,
           "power" = sprintf("E = exp(%.3f) * bp^%.3f", coefs[1], coefs[2]),
           "unknown"
    )
  }
  
  # calculate protocol parameters
  calculate_protocol <- function(target_energy, filtered_data) {
    raw_data <- filtered_data$raw
    
    # get parameter columns
    pip_values <- as.numeric(raw_data[, COLS$pip])
    duration_values <- as.numeric(raw_data[, COLS$duration])
    energy_values <- as.numeric(raw_data[, COLS$energy])
    
    # find valid protocols
    valid <- !is.na(pip_values) & !is.na(duration_values) & !is.na(energy_values) & 
      pip_values > 0 & duration_values > 0
    
    if(!any(valid)) {
      return(list(success = FALSE, message = "No valid protocols found"))
    }
    
    # find closest energy match
    energy_diff <- abs(energy_values[valid] - target_energy)
    closest_idx <- which.min(energy_diff)
    
    # get closest protocol parameters
    pip <- pip_values[valid][closest_idx]
    duration <- duration_values[valid][closest_idx]
    
    # calculate adjusted duty factor
    duty_factor <- (target_energy * 100) / (pip * duration)
    
    # validate duty factor
    warnings <- character()
    if(duty_factor > 100) {
      warnings <- c(warnings, "Duty factor exceeds 100% - consider increasing PIP or duration")
      duty_factor <- 100
    } else if(duty_factor < 1) {
      warnings <- c(warnings, "Duty factor very low (<1%) - consider decreasing PIP or duration")
    }
    
    list(
      success = TRUE,
      pip = pip,
      duration = duration,
      duty_factor = duty_factor,
      calculated_energy = pip * duration * (duty_factor / 100),
      warnings = warnings
    )
  }
  
  # file upload handler
  observeEvent(input$datafile, {
    req(input$datafile)
    
    df <- read.csv(input$datafile$datapath, stringsAsFactors = FALSE)
    
    # validate
    validation <- validate_csv(df)
    if(!validation$valid) {
      showNotification(validation$message, type = "error")
      return()
    }
    
    # add protocol type and store
    values$data <- add_protocol_type(df)
    
    # update device dropdown
    devices <- unique(df[, COLS$device])
    updateSelectInput(session, "device", choices = devices, selected = devices[1])
  })
  
  # device change handler
  observeEvent(input$device, {
    req(values$data, input$device)
    
    # get vessels for this device
    device_data <- values$data[values$data[, COLS$device] == input$device, ]
    vessels <- unique(device_data[, COLS$vessel])
    updateSelectInput(session, "vessel", choices = vessels, selected = vessels[1])
  })
  
  # vessel change handler
  observeEvent(input$vessel, {
    req(values$data, input$device, input$vessel)
    
    # get volumes for this device/vessel combination
    device_vessel_data <- values$data[values$data[, COLS$device] == input$device & 
                                        values$data[, COLS$vessel] == input$vessel, ]
    volumes <- sort(unique(as.numeric(device_vessel_data[, COLS$volume])))
    updateSelectInput(session, "volume", choices = volumes, selected = volumes[1])
  })
  
  # calculate button handler (updated)
  observeEvent(input$calculate, {
    req(values$data, input$device, input$vessel, input$volume)
    
    # filter data
    filtered <- filter_data(values$data, input$device, input$vessel, 
                            as.numeric(input$volume), input$protocol_type)
    
    if(filtered$n_points < 3) {
      showNotification("Insufficient data for modeling", type = "error")
      return()
    }
    
    # fit model
    fit_result <- fit_best_model(filtered$bp, filtered$energy)
    
    if(!fit_result$success) {
      showNotification(fit_result$message, type = "error")
      return()
    }
    
    values$current_model <- fit_result$model
    
    # calculate all models in background
    values$all_models <- calculate_all_models()
  })
  
  # calculate all models for export (updated)
  calculate_all_models <- reactive({
    req(values$data)
    
    combos <- unique(values$data[, c(COLS$device, COLS$vessel, COLS$volume)])
    results <- list()
    
    # show progress
    withProgress(message = 'Calculating equations', value = 0, {
      total_steps <- nrow(combos) * 3  # 3 protocol types per combo
      current_step <- 0
      
      for(i in 1:nrow(combos)) {
        device <- combos[i, 1]
        vessel <- combos[i, 2]
        volume <- as.numeric(combos[i, 3])
        
        for(protocol in c("All", "Continuous", "Pulsing")) {
          current_step <- current_step + 1
          incProgress(1/total_steps, detail = paste(device, "-", vessel, "-", volume, "μL -", protocol))
          
          # filter data
          if(protocol == "All") {
            filtered <- filter_data(values$data, device, vessel, volume)
          } else {
            filtered <- filter_data(values$data, device, vessel, volume, protocol)
          }
          
          # Initialize result entry with ALL possible columns
          result_row <- data.frame(
            device = device,
            vessel = vessel,
            volume = volume,
            protocol = protocol,
            n_points = filtered$n_points,
            n_pulsing = filtered$n_pulsing,
            n_continuous = filtered$n_continuous,
            model_type = NA_character_,
            aic = NA_real_,
            r_squared = NA_real_,
            adj_r_squared = NA_real_,
            rse = NA_real_,
            bp_min = NA_real_,
            bp_max = NA_real_,
            equation = NA_character_,
            coefficients = NA_character_,
            status = NA_character_,
            stringsAsFactors = FALSE
          )
          
          # fit model if enough data
          if(filtered$n_points >= 3) {
            fit_result <- fit_best_model(filtered$bp, filtered$energy)
            
            if(fit_result$success) {
              model <- fit_result$model
              coefs <- coef(model)
              
              result_row$model_type <- fit_result$name
              result_row$aic <- round(AIC(model), 2)
              result_row$r_squared <- round(summary(model)$r.squared, 3)
              result_row$adj_r_squared <- round(summary(model)$adj.r.squared, 3)
              result_row$rse <- round(summary(model)$sigma, 2)
              result_row$bp_min <- min(filtered$bp)
              result_row$bp_max <- max(filtered$bp)
              result_row$equation <- format_equation(model)
              result_row$coefficients <- paste(names(coefs), "=", round(coefs, 6), collapse = "; ")
              result_row$status <- "success"
            } else {
              result_row$status <- paste("error:", fit_result$message)
            }
          } else {
            result_row$status <- "insufficient data"
          }
          
          results[[length(results) + 1]] <- result_row
        }
      }
    })
    
    # Now all data frames have the same structure, so rbind will work
    do.call(rbind, results)
  })
  
  # model summary
  output$model_summary <- renderPrint({
    req(values$current_model)
    model <- values$current_model
    
    cat("Model:", attr(model, "model_type"), "\n")
    cat("R²:", round(summary(model)$r.squared, 3), "\n")
    cat("AIC:", round(AIC(model), 1), "\n")
    cat("Data points:", length(attr(model, "bp_data")), "\n")
    cat("BP range:", paste(round(attr(model, "data_range")), collapse = "-"), "\n")
  })
  
  # main plot (updated)
  # main plot (updated with error handling)
  output$energy_plot <- renderPlot({
    req(values$data, input$device, input$vessel, input$volume)
    
    # get all data for device/vessel/volume (not filtered by protocol)
    all_data <- values$data[values$data[, COLS$device] == input$device & 
                              values$data[, COLS$vessel] == input$vessel &
                              values$data[, COLS$volume] == as.numeric(input$volume), ]
    
    # extract bp and energy
    bp_all <- as.numeric(all_data[, COLS$bp_mode])
    energy_all <- as.numeric(all_data[, COLS$energy])
    protocol_all <- all_data$protocol_type
    
    # remove invalid
    valid <- !is.na(bp_all) & !is.na(energy_all) & energy_all > 0 & bp_all > 0
    bp_all <- bp_all[valid]
    energy_all <- energy_all[valid]
    protocol_all <- protocol_all[valid]
    
    # Check if we have any valid data points
    if(length(energy_all) == 0) {
      plot(NULL, NULL,
           main = paste(input$device, "-", input$vessel, "-", input$volume, "μL"),
           xlab = "Base Pair (bp)",
           ylab = "Energy (J)",
           xlim = c(100, 1500),
           ylim = c(0, 100))  # fallback range
      text(500, 50, "No valid data for this combination", cex = 1.2, col = "red")
      return()
    }
    
    # set up plot with valid data range
    plot(NULL, NULL,
         main = paste(input$device, "-", input$vessel, "-", input$volume, "μL"),
         xlab = "Base Pair (bp)",
         ylab = "Energy (J)",
         xlim = range(100, 1500),
         ylim = range(energy_all))  # Now safe to use since we checked length > 0
    
    # plot points based on protocol selection
    if(input$protocol_type == "Continuous") {
      idx <- protocol_all == "Continuous"
      if(any(idx)) {
        points(bp_all[idx], energy_all[idx], pch = 19, col = "blue")
      }
    } else if(input$protocol_type == "Pulsing") {
      idx <- protocol_all == "Pulsing"
      if(any(idx)) {
        points(bp_all[idx], energy_all[idx], pch = 17, col = "darkgreen")
      }
    }
    
    # add model if available
    if(!is.null(values$current_model)) {
      model <- values$current_model
      bp_seq <- seq(100, 1500, length.out = 100)
      pred <- predict_energy(model, bp_seq)
      
      # prediction interval
      polygon(c(bp_seq, rev(bp_seq)), 
              c(pred$pred_lower, rev(pred$pred_upper)),
              col = rgb(0.8, 0.8, 0.8, 0.3), border = NA)
      
      # confidence interval
      polygon(c(bp_seq, rev(bp_seq)), 
              c(pred$conf_lower, rev(pred$conf_upper)),
              col = rgb(0.5, 0.5, 0.5, 0.5), border = NA)
      
      # fitted line
      lines(bp_seq, pred$fit, col = "red", lwd = 2)
    }
    
    # target line
    abline(v = input$target_bp, col = "purple", lty = 2, lwd = 2)
    
    # legend
    legend_items <- c(input$protocol_type, "Fitted", "95% CI", "95% PI", "Target")
    legend_cols <- c(ifelse(input$protocol_type == "Pulsing", "darkgreen", "blue"), 
                     "red", rgb(0.5, 0.5, 0.5, 0.5), rgb(0.8, 0.8, 0.8, 0.3), "purple")
    legend_pch <- c(ifelse(input$protocol_type == "Pulsing", 17, 19), NA, NA, NA, NA)
    legend_lty <- c(NA, 1, NA, NA, 2)
    legend_fill <- c(NA, NA, rgb(0.5, 0.5, 0.5, 0.5), rgb(0.8, 0.8, 0.8, 0.3), NA)
    
    legend("topright", legend = legend_items, col = legend_cols, pch = legend_pch,
           lty = legend_lty, lwd = 2, fill = legend_fill)
    
    # add data count
    n_shown <- sum(protocol_all == input$protocol_type)
    mtext(paste("Showing", n_shown, "points"), side = 3, line = 0.5, cex = 0.8)
  })  
  # prediction output
  output$prediction <- renderPrint({
    req(values$current_model)
    
    pred <- predict_energy(values$current_model, input$target_bp)
    
    cat("Model:", attr(values$current_model, "model_type"), "\n")
    cat("Target BP:", input$target_bp, "bp\n\n")
    cat("Point Estimate:", round(pred$fit, 2), "J\n")
    cat("95% Confidence Interval: [", round(pred$conf_lower, 2), ",", round(pred$conf_upper, 2), "] J\n")
    cat("95% Prediction Interval: [", round(pred$pred_lower, 2), ",", round(pred$pred_upper, 2), "] J\n")
    
    if(pred$extrapolating) {
      cat("\n⚠️ WARNING: Extrapolating outside data range!\nPredictions may be unreliable.\n")
    }
  })
  
  # protocol output (updated)
  output$protocol <- renderPrint({
    req(values$current_model, values$data, input$volume)
    
    # get prediction
    pred <- predict_energy(values$current_model, input$target_bp)
    
    # filter data
    filtered <- filter_data(values$data, input$device, input$vessel, 
                            as.numeric(input$volume), input$protocol_type)
    
    # calculate protocol
    protocol <- calculate_protocol(pred$fit, filtered)
    
    if(!protocol$success) {
      cat(protocol$message)
      return()
    }
    
    cat("TARGET PROTOCOL PARAMETERS\n")
    cat("========================\n")
    cat("Target BP:", input$target_bp, "bp\n")
    cat("Target Energy:", round(pred$fit, 0), "J\n\n")
    
    cat("CALCULATED PROTOCOL\n")
    cat("-------------------\n")
    cat("Peak Incident Power:", round(protocol$pip, 1), "W\n")
    cat("Duration:", round(protocol$duration, 1), "s\n")
    cat("Duty Factor:", round(protocol$duty_factor, 0), "%\n")
    cat("Sample Volume:", input$volume, "μL\n")
    cat("Water Level: Placeholder \n")
    cat("Z-Height: Placeholder \n")
    cat("Dither: Placeholder \n")
    cat("Calculated Energy:", round(protocol$calculated_energy, 0), "J\n")
    
    # add pulsing parameters if applicable
    if(input$protocol_type == "Pulsing") {
      raw_data <- filtered$raw
      valid <- !is.na(raw_data[, COLS$energy])
      if(any(valid)) {
        # find closest match to get pulsing parameters
        energy_diff <- abs(as.numeric(raw_data[valid, COLS$energy]) - pred$fit)
        closest_idx <- which.min(energy_diff)
        
        cycles <- raw_data[valid, ][closest_idx, COLS$cycles]
        iterations <- raw_data[valid, ][closest_idx, COLS$iterations]
        
        cat("\nPULSING PARAMETERS\n")
        cat("------------------\n")
        cat("Cycles per Burst:", ifelse(is.na(cycles), "N/A", cycles), "\n")
        cat("Iterations:", ifelse(is.na(iterations), "N/A", iterations), "\n")
      }
    }
    
    if(length(protocol$warnings) > 0) {
      cat("\n⚠️", paste(protocol$warnings, collapse = "\n⚠️ "), "\n")
    }
    
    cat("\nVALIDATION\n")
    cat("----------\n")
    cat("Energy check:", round(protocol$pip, 1), "×", round(protocol$duration, 1), 
        "×", round(protocol$duty_factor, 0), "/100 =", round(protocol$calculated_energy, 0), "J\n")
  })
  
  # summary stats
  output$summary_stats <- renderPrint({
    req(values$all_models)
    
    models <- values$all_models
    successful <- models[models$status == "success", ]
    
    cat("MODEL SELECTION SUMMARY\n")
    cat("=======================\n")
    cat("Total combinations:", nrow(models), "\n")
    cat("Successful fits:", nrow(successful), "\n")
    cat("Failed fits:", nrow(models) - nrow(successful), "\n\n")
    
    if(nrow(successful) > 0) {
      cat("Successful Fits by Protocol:\n")
      protocol_table <- table(successful$protocol)
      for(p in names(protocol_table)) {
        cat("  ", p, ":", protocol_table[p], "\n")
      }
      
      cat("\nModel Type Distribution:\n")
      model_table <- table(successful$model_type)
      for(m in names(model_table)) {
        cat("  ", m, ":", model_table[m], 
            sprintf("(%.1f%%)\n", 100 * model_table[m] / nrow(successful)))
      }
      
      cat("\nModel Quality:\n")
      cat("  Average R²:", round(mean(successful$r_squared), 3), "\n")
      cat("  Average AIC:", round(mean(successful$aic), 1), "\n")
      
      # r-squared distribution
      cat("\nR² Distribution:\n")
      r2_breaks <- c(0, 0.7, 0.8, 0.9, 0.95, 1)
      r2_labels <- c("Poor (<0.7)", "Fair (0.7-0.8)", "Good (0.8-0.9)", 
                     "Very Good (0.9-0.95)", "Excellent (>0.95)")
      r2_cut <- cut(successful$r_squared, breaks = r2_breaks, labels = r2_labels)
      r2_table <- table(r2_cut)
      for(i in seq_along(r2_table)) {
        if(r2_table[i] > 0) {
          cat("  ", names(r2_table)[i], ":", r2_table[i], "\n")
        }
      }
    }
  })
  
  # equations table (updated)
  output$equations_table <- renderTable({
    req(values$all_models)
    
    # select columns to display (including volume)
    display_cols <- c("device", "vessel", "volume", "protocol", "model_type", 
                      "n_points", "n_pulsing", "n_continuous",
                      "aic", "r_squared", "equation", "status")
    
    # filter to columns that exist
    available_cols <- intersect(display_cols, names(values$all_models))
    values$all_models[, available_cols]
  })
  
  # download handler
  output$download_equations <- downloadHandler(
    filename = function() {
      paste0("dna_shearing_equations_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$all_models, file, row.names = FALSE)
    }
  )
  
  # example curves plot (updated)
  output$example_curves <- renderPlot({
    req(values$all_models)
    
    models_df <- values$all_models
    successful <- models_df[!is.na(models_df$model_type) & models_df$status == "success", ]
    
    if(nrow(successful) == 0) {
      plot(1, type = "n", xlab = "", ylab = "", main = "No successful models")
      return()
    }
    
    # get up to 4 different model types
    model_types <- unique(successful$model_type)
    n_types <- min(4, length(model_types))
    
    # set up 2x2 layout
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    
    for(i in 1:n_types) {
      # get example of this model type
      model_type <- model_types[i]
      example <- successful[successful$model_type == model_type, ][1, ]
      
      # get data for this device/vessel/volume/protocol
      if(example$protocol == "All") {
        filtered <- filter_data(values$data, example$device, example$vessel, example$volume)
      } else {
        filtered <- filter_data(values$data, example$device, example$vessel, example$volume, example$protocol)
      }
      
      if(filtered$n_points > 0) {
        # refit model for this example
        fit_result <- fit_best_model(filtered$bp, filtered$energy)
        
        if(fit_result$success) {
          # plot data
          plot(filtered$bp, filtered$energy,
               main = paste(example$device, "-", example$vessel, "-", example$volume, "μL\n", model_type, "model"),
               sub = paste("Protocol:", example$protocol),
               xlab = "Base Pair (bp)",
               ylab = "Energy (J)",
               pch = 19,
               col = "blue")
          
          # add fitted curve
          bp_seq <- seq(min(filtered$bp), max(filtered$bp), length.out = 100)
          pred <- predict_energy(fit_result$model, bp_seq)
          lines(bp_seq, pred$fit, col = "red", lwd = 2)
          
          # add confidence band
          polygon(c(bp_seq, rev(bp_seq)), 
                  c(pred$conf_lower, rev(pred$conf_upper)),
                  col = rgb(0.5, 0.5, 0.5, 0.2), border = NA)
          
          # add stats
          mtext(paste("AIC =", round(example$aic, 1), ", R² =", round(example$r_squared, 3)), 
                side = 3, line = -1.5, cex = 0.8)
          
          grid()
        }
      }
    }
    
    # reset layout
    par(mfrow = c(1, 1))
  })
  
  # download protocol handler (updated)
  output$download_protocol <- downloadHandler(
    filename = function() {
      req(input$device, input$vessel, input$volume, input$target_bp)
      paste0("DNA_Shearing_Protocol_", 
             gsub("[^A-Za-z0-9]", "_", input$device), "_",
             gsub("[^A-Za-z0-9]", "_", input$vessel), "_",
             input$volume, "uL_",
             input$target_bp, "bp_",
             Sys.Date(), ".txt")
    },
    
    content = function(file) {
      req(values$current_model, values$data, input$volume)
      
      # Get prediction and protocol data
      pred <- predict_energy(values$current_model, input$target_bp)
      filtered <- filter_data(values$data, input$device, input$vessel, 
                              as.numeric(input$volume), input$protocol_type)
      protocol <- calculate_protocol(pred$fit, filtered)
      
      if(!protocol$success) {
        writeLines(paste("Error:", protocol$message), file)
        return()
      }
      
      # Create formatted protocol text
      protocol_text <- c(
        "========================================",
        "     DNA SHEARING PROTOCOL",
        "========================================",
        "",
        paste("Generated:", Sys.time()),
        paste("Protocol Builder v1 - Covaris"),
        paste("Generated by:", Sys.getenv("USER", "User")),
        "",
        "EXPERIMENTAL PARAMETERS",
        "----------------------",
        paste("Device:", input$device),
        paste("Vessel:", input$vessel),
        paste("Sample Volume:", input$volume, "μL"),  # Added volume
        paste("Protocol Type:", input$protocol_type),
        paste("Target Fragment Size:", input$target_bp, "bp"),
        paste("Target Energy:", round(pred$fit, 0), "J"),
        "",
        "MODEL INFORMATION",
        "-----------------",
        paste("Model Type:", attr(values$current_model, "model_type")),
        paste("Model Equation:", format_equation(values$current_model)),
        paste("R-squared:", round(summary(values$current_model)$r.squared, 3)),
        paste("AIC:", round(AIC(values$current_model), 1)),
        paste("Training Data Points:", length(attr(values$current_model, "bp_data"))),
        paste("Data Range:", paste(round(attr(values$current_model, "data_range")), collapse = "-"), "bp"),
        ""
      )
      
      # Add prediction intervals
      protocol_text <- c(protocol_text,
                         "ENERGY PREDICTION",
                         "-----------------",
                         paste("Point Estimate:", round(pred$fit, 2), "J"),
                         paste("95% Confidence Interval: [", round(pred$conf_lower, 2), ",", 
                               round(pred$conf_upper, 2), "] J"),
                         paste("95% Prediction Interval: [", round(pred$pred_lower, 2), ",", 
                               round(pred$pred_upper, 2), "] J"),
                         ""
      )
      
      # Add extrapolation warning if needed
      if(pred$extrapolating) {
        protocol_text <- c(protocol_text,
                           "⚠️  WARNING: EXTRAPOLATION DETECTED",
                           "   Prediction is outside the range of training data.",
                           "   Results may be unreliable. Use with caution.",
                           ""
        )
      }
      
      # Add main protocol parameters
      protocol_text <- c(protocol_text,
                         "========================================",
                         "        RECOMMENDED PROTOCOL",
                         "========================================",
                         "",
                         "INSTRUMENT SETTINGS",
                         "-------------------",
                         paste("Peak Incident Power (PIP):", round(protocol$pip, 1), "W"),
                         paste("Duration:", round(protocol$duration, 1), "s"),
                         paste("Duty Factor:", round(protocol$duty_factor, 0), "%"),
                         paste("Sample Volume:", input$volume, "μL"),
                         paste("Water Level: [TO BE DETERMINED]"),
                         paste("Z-Height: [TO BE DETERMINED]"),
                         paste("Dither: [TO BE DETERMINED]"),
                         "",
                         "CALCULATED VALUES",
                         "-----------------",
                         paste("Calculated Energy:", round(protocol$calculated_energy, 0), "J"),
                         paste("Energy Check:", round(protocol$pip, 1), "×", round(protocol$duration, 1), 
                               "×", round(protocol$duty_factor, 0), "/100 =", 
                               round(protocol$calculated_energy, 0), "J"),
                         ""
      )
      
      # Add pulsing parameters if applicable
      if(input$protocol_type == "Pulsing") {
        raw_data <- filtered$raw
        valid <- !is.na(raw_data[, COLS$energy])
        if(any(valid)) {
          energy_diff <- abs(as.numeric(raw_data[valid, COLS$energy]) - pred$fit)
          closest_idx <- which.min(energy_diff)
          
          cycles <- raw_data[valid, ][closest_idx, COLS$cycles]
          iterations <- raw_data[valid, ][closest_idx, COLS$iterations]
          
          protocol_text <- c(protocol_text,
                             "PULSING PARAMETERS",
                             "------------------",
                             paste("Cycles per Burst:", ifelse(is.na(cycles), "N/A", cycles)),
                             paste("Iterations:", ifelse(is.na(iterations), "N/A", iterations)),
                             ""
          )
        }
      }
      
      # Add warnings if any
      if(length(protocol$warnings) > 0) {
        protocol_text <- c(protocol_text,
                           "⚠️  WARNINGS",
                           "-----------",
                           paste("⚠️ ", protocol$warnings),
                           ""
        )
      }
      
      # Add footer
      protocol_text <- c(protocol_text,
                         "========================================",
                         "             IMPORTANT NOTES",
                         "========================================",
                         "",
                         "1. Always verify water level and Z-height before starting",
                         "2. Ensure sample volume is appropriate for vessel type",
                         "3. Monitor temperature during long treatments",
                         "4. This protocol is based on predictive modeling",
                         "5. Validate results with analytical methods (gel, Bioanalyzer, etc.)",
                         "",
                         "For technical support, contact: bpfeiffer@covaris.com",
                         "",
                         paste("Protocol generated by DNA Shearing Protocol Builder v1"),
                         paste("File created:", Sys.time())
      )
      
      # Write to file
      writeLines(protocol_text, file)
    }
  )
}

# run app
shinyApp(ui = ui, server = server)
