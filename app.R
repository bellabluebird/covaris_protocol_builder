# dna shearing protocol builder v2
# built by bella pfeiffer (bpfeiffer@covaris.com)
# available here: https://bellabluebird.shinyapps.io/covaris_protocol_builder/
# last updated: 9/2/25

# load libraries
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}
pacman::p_load(shiny, rsconnect, DBI, RSQLite, dplyr, digest, DT)

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

# password setup
password <- Sys.getenv("COVARIS_PASSWORD", "covaris2025")
password_hash <- digest(password, algo = "sha256")

# ui definition
ui <- fluidPage(
  tags$head(
    tags$title("🧬 Covaris Protocol Builder")
  ),
  
  uiOutput("main_ui")
)

# server logic
server <- function(input, output, session) {
  
  # reactive values to store state
  values <- reactiveValues(
    authenticated = FALSE,
    data = NULL,
    current_model = NULL,
    all_models = NULL,
    data_loaded = FALSE,
    error_message = NULL,
    comp_source_energy = NULL,
    comp_vessel_options = NULL,
    comp_results_ready = FALSE
  )
  
  # helper functions
  
  # check password
  check_password <- function(pwd) {
    if(is.null(pwd) || pwd == "") return(FALSE)
    digest(pwd, algo = "sha256") == password_hash
  }
  
  # get database connection with error handling
  get_database_data <- function() {
    possible_paths <- c(
      "data/database.sqlite",
      "./data/database.sqlite", 
      "database.sqlite",
      "./database.sqlite"
    )
    
    for(path in possible_paths) {
      if(file.exists(path)) {
        tryCatch({
          con <- dbConnect(RSQLite::SQLite(), path)
          
          if(!"data_protocols" %in% dbListTables(con)) {
            dbDisconnect(con)
            next
          }
          
          df <- dbReadTable(con, "data_protocols")
          dbDisconnect(con)
          
          if(ncol(df) < 17) {
            next
          }
          
          return(list(success = TRUE, data = df))
        }, error = function(e) {
          if(exists("con")) dbDisconnect(con)
          next
        })
      }
    }
    
    return(list(success = FALSE, message = "Database file not found or invalid"))
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
    
    # fit power model (only model type for now)
    tryCatch({
      model <- lm(log(energy_j) ~ log(bp_mode), data = df_model)
      
      # add metadata
      attr(model, "model_type") <- "power"
      attr(model, "data_range") <- range(bp)
      attr(model, "bp_data") <- bp
      attr(model, "energy_data") <- energy
      
      list(success = TRUE, model = model, name = "power")
    }, error = function(e) {
      list(success = FALSE, message = e$message)
    })
  }
  
  # predict with intervals
  predict_energy <- function(model, bp_target) {
    pred_df <- data.frame(bp_mode = bp_target)
    
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
    sprintf("E = exp(%.3f) * bp^%.3f", coefs[1], coefs[2])
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
  
  # PROTOCOL COMPARISON HELPER FUNCTIONS
  
  # calculate energy from manual protocol input
  calculate_energy_from_protocol <- function(pip, duration, duty_factor) {
    pip * duration * (duty_factor / 100)
  }
  
  # find compatible volumes on target device
  find_compatible_volumes <- function(data, target_device, source_volume, protocol_type) {
    # get all vessel/volume combos for target device
    device_data <- data[data[, COLS$device] == target_device, ]
    
    if(nrow(device_data) == 0) {
      return(list(success = FALSE, message = "No data for target device"))
    }
    
    # get unique vessel/volume combinations with sufficient data
    vessel_volume_combos <- unique(device_data[, c(COLS$vessel, COLS$volume)])
    vessel_volume_combos$volume_numeric <- as.numeric(vessel_volume_combos[, 2])
    
    # calculate volume differences
    vessel_volume_combos$volume_diff <- abs(vessel_volume_combos$volume_numeric - source_volume)
    vessel_volume_combos$volume_ratio <- vessel_volume_combos$volume_numeric / source_volume
    
    # check data availability for each combo
    vessel_volume_combos$data_points <- 0
    vessel_volume_combos$has_protocol_type <- FALSE
    
    for(i in 1:nrow(vessel_volume_combos)) {
      filtered <- filter_data(data, target_device, 
                              vessel_volume_combos[i, 1], 
                              vessel_volume_combos[i, 2],
                              protocol_type)
      vessel_volume_combos$data_points[i] <- filtered$n_points
      vessel_volume_combos$has_protocol_type[i] <- filtered$n_points >= 3
    }
    
    # sort by volume similarity and data availability
    vessel_volume_combos <- vessel_volume_combos[order(
      !vessel_volume_combos$has_protocol_type,  # prioritize those with data
      vessel_volume_combos$volume_diff          # then by volume similarity
    ), ]
    
    list(
      success = TRUE,
      options = vessel_volume_combos,
      exact_match = any(vessel_volume_combos$volume_diff == 0),
      closest_volume = vessel_volume_combos$volume_numeric[1]
    )
  }
  
  # translate protocol between devices
  translate_protocol <- function(source_energy, target_device, target_vessel, 
                                 target_volume, protocol_type, data) {
    # filter data for target combination
    filtered <- filter_data(data, target_device, target_vessel, target_volume, protocol_type)
    
    if(filtered$n_points < 3) {
      return(list(
        success = FALSE, 
        message = "Insufficient data for target device/vessel/volume combination"
      ))
    }
    
    # fit model for target
    fit_result <- fit_best_model(filtered$bp, filtered$energy)
    
    if(!fit_result$success) {
      return(list(success = FALSE, message = "Failed to fit model for target"))
    }
    
    # reverse engineer BP from energy
    model <- fit_result$model
    coefs <- coef(model)
    
    # for power model: E = exp(a) * BP^b
    # so: BP = (E / exp(a))^(1/b)
    estimated_bp <- (source_energy / exp(coefs[1]))^(1/coefs[2])
    
    # calculate protocol for this energy
    protocol <- calculate_protocol(source_energy, filtered)
    
    if(!protocol$success) {
      return(list(success = FALSE, message = protocol$message))
    }
    
    # add additional info
    protocol$estimated_bp <- estimated_bp
    protocol$target_energy <- source_energy
    protocol$model_r2 <- summary(model)$r.squared
    protocol$data_points <- filtered$n_points
    
    protocol
  }
  
  # authentication handlers
  observeEvent(input$login_btn, {
    if(check_password(input$password)) {
      values$authenticated <- TRUE
      updateTextInput(session, "password", value = "")
    } else {
      updateTextInput(session, "password", value = "")
      showNotification("Incorrect password", type = "error")
    }
  })
  
  observeEvent(input$logout_btn, {
    values$authenticated <- FALSE
  })
  
  # main ui renderer
  output$main_ui <- renderUI({
    if(!values$authenticated) {
      div(style = "max-width: 400px; margin: 100px auto; padding: 30px; border: 1px solid #ddd; border-radius: 10px;",
          h2("🧬 Covaris Protocol Builder", style = "text-align: center;"),
          p("Please enter password:", style = "text-align: center; color: #666;"),
          br(),
          passwordInput("password", "Password:", width = "100%"),
          br(),
          actionButton("login_btn", "Login", class = "btn-primary", style = "width: 100%;")
      )
    } else {
      tagList(
        div(style = "position: absolute; top: 10px; right: 10px;",
            actionButton("logout_btn", "Logout", class = "btn-default btn-sm")
        ),
        
        titlePanel("Covaris Protocol Builder"),
        
        sidebarLayout(
          sidebarPanel(
            # show data loading status
            verbatimTextOutput("data_status"),
            hr(),
            
            # dynamic dropdowns
            selectInput("device", "Select Device:", 
                        choices = list("Loading data..." = "")),
            selectInput("vessel", "Select Vessel:", 
                        choices = list("Select device first" = "")),
            selectInput("volume", "Select Sample Volume:",
                        choices = list("Select vessel first" = "")),
            radioButtons("protocol_type", "Protocol Type:",
                         choices = c("Continuous" = "Continuous",
                                     "Pulsing (Majority Not Available)" = "Pulsing"),
                         selected = "Continuous"),
            
            # target selection
            sliderInput("target_bp", "Target Base Pair (bp):", 
                        min = 100, max = 1500, value = 500, step = 25),
            actionButton("calculate", "Calculate Protocol", class = "btn-primary"),
            
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
              
              # protocol comparison tab
              tabPanel("Protocol Comparison",
                       # Step 1: Source Protocol
                       fluidRow(
                         column(12,
                                wellPanel(
                                  h4("Step 1: Define Your Current Protocol"),
                                  
                                  fluidRow(
                                    column(3,
                                           selectInput("comp_source_device", "Device:",
                                                       choices = list("Loading data..." = ""),
                                                       width = "100%")),
                                    column(3,
                                           selectInput("comp_source_vessel", "Vessel:",
                                                       choices = list("Select device first" = ""),
                                                       width = "100%")),
                                    column(3,
                                           selectInput("comp_source_volume", "Volume:",
                                                       choices = list("Select vessel first" = ""),
                                                       width = "100%")),
                                    column(3,
                                           radioButtons("comp_protocol_type", "Type:",
                                                        choices = c("Continuous" = "Continuous",
                                                                    "Pulsing" = "Pulsing"),
                                                        selected = "Continuous",
                                                        inline = TRUE))
                                  ),
                                  
                                  hr(),
                                  
                                  fluidRow(
                                    column(6,
                                           radioButtons("comp_input_method", "How do you know your protocol?",
                                                        choices = c("I know my target BP" = "bp",
                                                                    "I know my exact settings" = "manual"),
                                                        selected = "bp"),
                                           conditionalPanel(
                                             condition = "input.comp_input_method == 'bp'",
                                             sliderInput("comp_target_bp", "Target Fragment Size:", 
                                                         min = 100, max = 1500, value = 500, 
                                                         step = 25, post = " bp", width = "100%")
                                           ),
                                           conditionalPanel(
                                             condition = "input.comp_input_method == 'manual'",
                                             fluidRow(
                                               column(4, numericInput("comp_manual_pip", "PIP (W):", 
                                                                      value = 175, width = "100%")),
                                               column(4, numericInput("comp_manual_duration", "Duration (s):", 
                                                                      value = 60, width = "100%")),
                                               column(4, numericInput("comp_manual_duty", "Duty (%):", 
                                                                      value = 10, width = "100%"))
                                             )
                                           )
                                    ),
                                    column(6,
                                           div(style = "background-color: #fff3cd; padding: 15px; border-radius: 5px;",
                                               h5("Quick Tip"),
                                               p("Select 'target BP' if you're designing a new protocol."),
                                               p("Select 'exact settings' if you want to replicate an existing protocol on a different device."),
                                               uiOutput("comp_energy_display")
                                           )
                                    )
                                  )
                                )
                         )
                       ),
                       
                       # Step 2: Target Device
                       fluidRow(
                         column(12,
                                wellPanel(
                                  h4("Step 2: Select Target Device"),
                                  
                                  fluidRow(
                                    column(4,
                                           selectInput("comp_target_device", "Target Device:",
                                                       choices = list("Select source first" = ""),
                                                       width = "100%")),
                                    column(8,
                                           div(style = "padding-top: 25px;",
                                               actionButton("comp_calculate", 
                                                            "Compare Protocols", 
                                                            class = "btn-primary btn-lg", 
                                                            style = "width: 100%;"))
                                    )
                                  )
                                )
                         )
                       ),
                       
                       # Results Section
                       conditionalPanel(
                         condition = "output.show_comparison_results",
                         
                         fluidRow(
                           column(12,
                                  wellPanel(
                                    h4("Comparison Results"),
                                    
                                    fluidRow(
                                      column(6,
                                             h5("Volume Compatibility"),
                                             uiOutput("comp_volume_status"),
                                             hr(),
                                             h6("Available Options:"),
                                             DT::dataTableOutput("comp_vessel_options_dt")
                                      ),
                                      column(6,
                                             h5("Recommended Protocol"),
                                             uiOutput("comp_protocol_display"),
                                             hr(),
                                             fluidRow(
                                               column(6,
                                                      downloadButton("download_comparison", 
                                                                     "Download Report",
                                                                     class = "btn-success", 
                                                                     style = "width: 100%;")),
                                               column(6,
                                                      actionButton("copy_protocol", 
                                                                   "Copy to Clipboard",
                                                                   class = "btn-info", 
                                                                   style = "width: 100%;"))
                                             )
                                      )
                                    ),
                                    
                                    fluidRow(
                                      column(12,
                                             h5("Side-by-Side Comparison"),
                                             DT::dataTableOutput("comp_summary_visual")
                                      )
                                    )
                                  )
                           )
                         )
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
    }
  })
  
  # load data from database
  observe({
    req(values$authenticated)
    
    result <- get_database_data()
    
    if(result$success) {
      values$data <- add_protocol_type(result$data)
      values$data_loaded <- TRUE
      values$error_message <- NULL
      
      # update device dropdown
      devices <- unique(values$data[, COLS$device])
      devices <- devices[!is.na(devices)]
      if(length(devices) > 0) {
        updateSelectInput(session, "device", choices = devices, selected = devices[1])
        # also update comparison source device dropdown
        updateSelectInput(session, "comp_source_device", choices = devices, selected = devices[1])
      } else {
        values$error_message <- "No valid devices found in database"
      }
    } else {
      values$error_message <- result$message
      values$data_loaded <- FALSE
    }
  })
  
  # data loading status output
  output$data_status <- renderText({
    req(values$authenticated)
    if(values$data_loaded) {
      paste("✓ Data loaded successfully:", nrow(values$data), "records")
    } else if(!is.null(values$error_message)) {
      paste("✗ Error:", values$error_message)
    } else {
      "Loading data..."
    }
  })
  
  # device change handler
  observeEvent(input$device, {
    req(values$data, input$device, values$data_loaded, values$authenticated)
    
    # get vessels for this device
    device_data <- values$data[values$data[, COLS$device] == input$device, ]
    vessels <- unique(device_data[, COLS$vessel])
    vessels <- vessels[!is.na(vessels)]
    
    if(length(vessels) > 0) {
      updateSelectInput(session, "vessel", choices = vessels, selected = vessels[1])
    }
  })
  
  # vessel change handler
  observeEvent(input$vessel, {
    req(values$data, input$device, input$vessel, values$data_loaded, values$authenticated)
    
    # get volumes for this device/vessel combination
    device_vessel_data <- values$data[values$data[, COLS$device] == input$device & 
                                        values$data[, COLS$vessel] == input$vessel, ]
    volumes <- sort(unique(as.numeric(device_vessel_data[, COLS$volume])))
    volumes <- volumes[!is.na(volumes)]
    
    if(length(volumes) > 0) {
      updateSelectInput(session, "volume", choices = volumes, selected = volumes[1])
    }
  })
  
  # calculate button handler
  observeEvent(input$calculate, {
    req(values$data, input$device, input$vessel, input$volume, values$data_loaded, values$authenticated)
    
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
  
  # PROTOCOL COMPARISON EVENT HANDLERS
  
  # source device change handler
  observeEvent(input$comp_source_device, {
    req(values$data, input$comp_source_device, values$data_loaded, values$authenticated)
    
    # get vessels for this device
    device_data <- values$data[values$data[, COLS$device] == input$comp_source_device, ]
    vessels <- unique(device_data[, COLS$vessel])
    vessels <- vessels[!is.na(vessels)]
    updateSelectInput(session, "comp_source_vessel", choices = vessels, selected = vessels[1])
    
    # update target device dropdown (exclude source device)
    all_devices <- unique(values$data[, COLS$device])
    other_devices <- all_devices[all_devices != input$comp_source_device]
    updateSelectInput(session, "comp_target_device", choices = other_devices, 
                      selected = if(length(other_devices) > 0) other_devices[1] else "")
  })
  
  # source vessel change handler
  observeEvent(input$comp_source_vessel, {
    req(values$data, input$comp_source_device, input$comp_source_vessel, values$data_loaded, values$authenticated)
    
    # get volumes for this device/vessel combination
    device_vessel_data <- values$data[values$data[, COLS$device] == input$comp_source_device & 
                                        values$data[, COLS$vessel] == input$comp_source_vessel, ]
    volumes <- sort(unique(as.numeric(device_vessel_data[, COLS$volume])))
    volumes <- volumes[!is.na(volumes)]
    updateSelectInput(session, "comp_source_volume", choices = volumes, selected = volumes[1])
  })
  
  # comparison calculation handler
  observeEvent(input$comp_calculate, {
    req(values$data, input$comp_source_device, input$comp_source_vessel, 
        input$comp_source_volume, input$comp_target_device, values$data_loaded, values$authenticated)
    
    # calculate source energy
    if(input$comp_input_method == "manual") {
      # direct calculation from manual input
      source_energy <- calculate_energy_from_protocol(
        input$comp_manual_pip,
        input$comp_manual_duration,
        input$comp_manual_duty
      )
      values$comp_source_energy <- source_energy
    } else {
      # calculate from BP target
      filtered <- filter_data(values$data, input$comp_source_device, 
                              input$comp_source_vessel, as.numeric(input$comp_source_volume), 
                              input$comp_protocol_type)
      
      if(filtered$n_points < 3) {
        showNotification("Insufficient data for source device/vessel/volume", type = "error")
        return()
      }
      
      # fit model and predict
      fit_result <- fit_best_model(filtered$bp, filtered$energy)
      if(!fit_result$success) {
        showNotification("Failed to fit model for source", type = "error")
        return()
      }
      
      pred <- predict_energy(fit_result$model, input$comp_target_bp)
      values$comp_source_energy <- pred$fit
    }
    
    # find compatible volumes on target device
    volume_options <- find_compatible_volumes(
      values$data, 
      input$comp_target_device,
      as.numeric(input$comp_source_volume),
      input$comp_protocol_type
    )
    
    if(!volume_options$success) {
      showNotification(volume_options$message, type = "error")
      return()
    }
    
    values$comp_vessel_options <- volume_options$options
    values$comp_results_ready <- TRUE
  })
  
  # calculate all models for export
  calculate_all_models <- reactive({
    req(values$data, values$data_loaded, values$authenticated)
    
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
          
          # initialize result entry
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
    
    do.call(rbind, results)
  })
  
  # PROTOCOL COMPARISON OUTPUT FUNCTIONS
  
  # show/hide results section
  output$show_comparison_results <- reactive({
    values$comp_results_ready
  })
  outputOptions(output, "show_comparison_results", suspendWhenHidden = FALSE)
  
  # energy display
  output$comp_energy_display <- renderUI({
    if(input$comp_input_method == "manual") {
      energy <- calculate_energy_from_protocol(
        input$comp_manual_pip,
        input$comp_manual_duration,
        input$comp_manual_duty
      )
      div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 15px;",
          h6("Calculated Energy:"),
          p(style = "font-size: 18px; margin: 0;", 
            strong(round(energy, 0), " J")))
    } else {
      if(!is.null(input$comp_target_bp)) {
        div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 15px;",
            h6("Target Fragment Size:"),
            p(style = "font-size: 18px; margin: 0;", 
              strong(input$comp_target_bp, " bp")))
      }
    }
  })
  
  # volume status display
  output$comp_volume_status <- renderUI({
    req(values$comp_vessel_options)
    
    options <- values$comp_vessel_options
    source_vol <- as.numeric(input$comp_source_volume)
    
    if(options[1, "volume_diff"] == 0) {
      div(
        p(style = "color: #28a745; font-weight: bold;", 
          "✓ Exact volume match found!"),
        p(paste("Using", options[1, "volume_numeric"], "μL in", options[1, 1]))
      )
    } else {
      div(
        p(style = "color: #ffc107; font-weight: bold;", 
          "⚠ No exact volume match"),
        p(paste("Source:", source_vol, "μL → Closest:", 
                options[1, "volume_numeric"], "μL")),
        p(paste("Volume ratio:", round(options[1, "volume_ratio"], 2))),
        p(style = "font-size: 12px; color: #6c757d;", 
          "Note: Energy calculations assume similar efficiency across volumes.")
      )
    }
  })
  
  # vessel options DataTable
  output$comp_vessel_options_dt <- DT::renderDataTable({
    req(values$comp_vessel_options)
    
    options <- values$comp_vessel_options
    display_df <- data.frame(
      Vessel = options[, 1],
      Volume = paste0(options$volume_numeric, " μL"),
      `Match` = ifelse(options$volume_diff == 0, "Exact", 
                       paste0("+", round(options$volume_diff, 1), " μL")),
      `Points` = options$data_points,
      Status = ifelse(options$has_protocol_type, 
                      '<span style="color: #28a745;">✓ Ready</span>', 
                      '<span style="color: #dc3545;">✗ No data</span>'),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(head(display_df, 5), 
                  options = list(
                    dom = 't',
                    pageLength = 5,
                    ordering = FALSE,
                    searching = FALSE,
                    paging = FALSE,
                    info = FALSE
                  ),
                  escape = FALSE,
                  rownames = FALSE)
  }, server = FALSE)
  
  # protocol display
  output$comp_protocol_display <- renderUI({
    req(values$comp_source_energy, values$comp_vessel_options)
    
    best_option <- values$comp_vessel_options[1, ]
    
    if(!best_option$has_protocol_type) {
      return(div(
        p(style = "color: #dc3545; font-weight: bold;", 
          "✗ Insufficient data"),
        p(paste("No protocol data for", input$comp_target_device, 
                "-", best_option[, 1], "-", best_option$volume_numeric, "μL"))
      ))
    }
    
    # translate protocol
    result <- translate_protocol(
      values$comp_source_energy,
      input$comp_target_device,
      best_option[, 1],
      best_option$volume_numeric,
      input$comp_protocol_type,
      values$data
    )
    
    if(!result$success) {
      return(div(
        p(style = "color: #dc3545; font-weight: bold;", 
          "✗ Translation failed"),
        p(result$message)
      ))
    }
    
    # create protocol display
    tagList(
      div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
          p(style = "color: #28a745; font-weight: bold; margin: 0;",
            "✓ Protocol successfully translated!")),
      
      fluidRow(
        column(6,
               p(strong("PIP:"), br(), 
                 span(style = "font-size: 18px;", round(result$pip, 1), " W"))),
        column(6,
               p(strong("Duration:"), br(), 
                 span(style = "font-size: 18px;", round(result$duration, 1), " s")))
      ),
      fluidRow(
        column(6,
               p(strong("Duty Factor:"), br(), 
                 span(style = "font-size: 18px;", round(result$duty_factor, 0), "%"))),
        column(6,
               p(strong("Energy:"), br(), 
                 span(style = "font-size: 18px;", round(result$calculated_energy, 0), " J")))
      ),
      
      hr(),
      
      div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px;",
          p(strong("Confidence Metrics:")),
          p(style = "margin: 5px 0;", 
            "Model R²: ", strong(round(result$model_r2, 3))),
          p(style = "margin: 5px 0;", 
            "Data points: ", strong(result$data_points)),
          p(style = "margin: 5px 0;", 
            "Est. BP: ", strong(round(result$estimated_bp, 0), " bp"))
      ),
      
      if(length(result$warnings) > 0) {
        div(style = "margin-top: 10px;",
            p(style = "color: #ffc107; font-weight: bold;", 
              "⚠ Warnings:"),
            tags$ul(
              lapply(result$warnings, function(w) tags$li(w))
            ))
      }
    )
  })
  
  # comparison summary visual
  output$comp_summary_visual <- DT::renderDataTable({
    req(values$comp_source_energy, values$comp_vessel_options)
    
    # create comparison data
    comparison_data <- data.frame(
      Parameter = character(),
      Source = character(),
      Target = character(),
      stringsAsFactors = FALSE
    )
    
    # basic info
    comparison_data <- rbind(comparison_data,
                             data.frame(
                               Parameter = "Device",
                               Source = input$comp_source_device,
                               Target = input$comp_target_device
                             ))
    
    comparison_data <- rbind(comparison_data,
                             data.frame(
                               Parameter = "Vessel",
                               Source = input$comp_source_vessel,
                               Target = values$comp_vessel_options[1, 1]
                             ))
    
    comparison_data <- rbind(comparison_data,
                             data.frame(
                               Parameter = "Volume",
                               Source = paste0(input$comp_source_volume, " μL"),
                               Target = paste0(values$comp_vessel_options[1, "volume_numeric"], " μL")
                             ))
    
    comparison_data <- rbind(comparison_data,
                             data.frame(
                               Parameter = "Protocol Type",
                               Source = input$comp_protocol_type,
                               Target = input$comp_protocol_type
                             ))
    
    # energy info
    comparison_data <- rbind(comparison_data,
                             data.frame(
                               Parameter = "Target Energy",
                               Source = paste0(round(values$comp_source_energy, 0), " J"),
                               Target = ""
                             ))
    
    # try to get protocol details
    best_option <- values$comp_vessel_options[1, ]
    if(best_option$has_protocol_type) {
      result <- translate_protocol(
        values$comp_source_energy,
        input$comp_target_device,
        best_option[, 1],
        best_option$volume_numeric,
        input$comp_protocol_type,
        values$data
      )
      
      if(result$success) {
        comparison_data[comparison_data$Parameter == "Target Energy", "Target"] <- 
          paste0(round(result$calculated_energy, 0), " J")
        
        # add protocol parameters
        comparison_data <- rbind(comparison_data,
                                 data.frame(
                                   Parameter = "PIP",
                                   Source = ifelse(input$comp_input_method == "manual",
                                                   paste0(input$comp_manual_pip, " W"),
                                                   "Calculated"),
                                   Target = paste0(round(result$pip, 1), " W")
                                 ))
        
        comparison_data <- rbind(comparison_data,
                                 data.frame(
                                   Parameter = "Duration",
                                   Source = ifelse(input$comp_input_method == "manual",
                                                   paste0(input$comp_manual_duration, " s"),
                                                   "Calculated"),
                                   Target = paste0(round(result$duration, 1), " s")
                                 ))
        
        comparison_data <- rbind(comparison_data,
                                 data.frame(
                                   Parameter = "Duty Factor",
                                   Source = ifelse(input$comp_input_method == "manual",
                                                   paste0(input$comp_manual_duty, "%"),
                                                   "Calculated"),
                                   Target = paste0(round(result$duty_factor, 0), "%")
                                 ))
      }
    }
    
    DT::datatable(comparison_data,
                  options = list(
                    dom = 't',
                    pageLength = 20,
                    ordering = FALSE,
                    searching = FALSE,
                    paging = FALSE,
                    info = FALSE
                  ),
                  rownames = FALSE)
  }, server = FALSE)
  
  # model summary
  output$model_summary <- renderPrint({
    req(values$current_model, values$authenticated)
    model <- values$current_model
    
    cat("Model:", attr(model, "model_type"), "\n")
    cat("R²:", round(summary(model)$r.squared, 3), "\n")
    cat("AIC:", round(AIC(model), 1), "\n")
    cat("Data points:", length(attr(model, "bp_data")), "\n")
    cat("BP range:", paste(round(attr(model, "data_range")), collapse = "-"), "\n")
  })
  
  # main plot
  output$energy_plot <- renderPlot({
    req(values$data, input$device, input$vessel, input$volume, values$data_loaded, values$authenticated)
    
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
    
    # check if we have any valid data points
    if(length(energy_all) == 0) {
      plot(NULL, NULL,
           main = paste(input$device, "-", input$vessel, "-", input$volume, "μL"),
           xlab = "Base Pair (bp)",
           ylab = "Energy (J)",
           xlim = c(100, 1500),
           ylim = c(0, 100))
      text(500, 50, "No valid data for this combination", cex = 1.2, col = "red")
      return()
    }
    
    # set up plot with valid data range
    plot(NULL, NULL,
         main = paste(input$device, "-", input$vessel, "-", input$volume, "μL"),
         xlab = "Base Pair (bp)",
         ylab = "Energy (J)",
         xlim = range(100, 1500),
         ylim = range(energy_all))
    
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
    req(values$current_model, values$authenticated)
    
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
  
  # protocol output
  output$protocol <- renderPrint({
    req(values$current_model, values$data, input$volume, values$authenticated)
    
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
    req(values$all_models, values$authenticated)
    
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
  
  # equations table
  output$equations_table <- renderTable({
    req(values$all_models, values$authenticated)
    
    # select columns to display
    display_cols <- c("device", "vessel", "volume", "protocol", "model_type", 
                      "n_points", "n_pulsing", "n_continuous",
                      "aic", "r_squared", "equation", "status")
    
    # filter to columns that exist
    available_cols <- intersect(display_cols, names(values$all_models))
    values$all_models[, available_cols]
  })
  
  # download equations handler
  output$download_equations <- downloadHandler(
    filename = function() {
      paste0("dna_shearing_equations_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$authenticated)
      write.csv(values$all_models, file, row.names = FALSE)
    }
  )
  
  # example curves plot
  output$example_curves <- renderPlot({
    req(values$all_models, values$authenticated)
    
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
  
  # download protocol handler
  output$download_protocol <- downloadHandler(
    filename = function() {
      req(input$device, input$vessel, input$volume, input$target_bp, values$authenticated)
      paste0("DNA_Shearing_Protocol_", 
             gsub("[^A-Za-z0-9]", "_", input$device), "_",
             gsub("[^A-Za-z0-9]", "_", input$vessel), "_",
             input$volume, "uL_",
             input$target_bp, "bp_",
             Sys.Date(), ".txt")
    },
    
    content = function(file) {
      req(values$current_model, values$data, input$volume, values$authenticated)
      
      # get prediction and protocol data
      pred <- predict_energy(values$current_model, input$target_bp)
      filtered <- filter_data(values$data, input$device, input$vessel, 
                              as.numeric(input$volume), input$protocol_type)
      protocol <- calculate_protocol(pred$fit, filtered)
      
      if(!protocol$success) {
        writeLines(paste("Error:", protocol$message), file)
        return()
      }
      
      # create formatted protocol text
      protocol_text <- c(
        "========================================",
        "     DNA SHEARING PROTOCOL",
        "========================================",
        "",
        paste("Generated:", Sys.time()),
        paste("Protocol Builder v2 - Covaris"),
        paste("Generated by:", Sys.getenv("USER", "User")),
        "",
        "EXPERIMENTAL PARAMETERS",
        "----------------------",
        paste("Device:", input$device),
        paste("Vessel:", input$vessel),
        paste("Sample Volume:", input$volume, "μL"),
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
      
      # add prediction intervals
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
      
      # add extrapolation warning if needed
      if(pred$extrapolating) {
        protocol_text <- c(protocol_text,
                           "⚠️  WARNING: EXTRAPOLATION DETECTED",
                           "   Prediction is outside the range of training data.",
                           "   Results may be unreliable. Use with caution.",
                           ""
        )
      }
      
      # add main protocol parameters
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
      
      # add pulsing parameters if applicable
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
      
      # add warnings if any
      if(length(protocol$warnings) > 0) {
        protocol_text <- c(protocol_text,
                           "⚠️  WARNINGS",
                           "-----------",
                           paste("⚠️ ", protocol$warnings),
                           ""
        )
      }
      
      # add footer
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
                         paste("Protocol generated by DNA Shearing Protocol Builder v2"),
                         paste("File created:", Sys.time())
      )
      
      # write to file
      writeLines(protocol_text, file)
    }
  )
  
  # PROTOCOL COMPARISON DOWNLOAD HANDLERS
  
  # download comparison report
  output$download_comparison <- downloadHandler(
    filename = function() {
      paste0("Protocol_Comparison_",
             input$comp_source_device, "_to_",
             input$comp_target_device, "_",
             Sys.Date(), ".txt")
    },
    content = function(file) {
      req(values$authenticated)
      # comprehensive comparison report
      report_lines <- c(
        "========================================",
        "     PROTOCOL COMPARISON REPORT",
        "========================================",
        "",
        paste("Generated:", Sys.time()),
        paste("Protocol Builder v2 - Covaris"),
        "",
        "SOURCE PROTOCOL",
        "---------------",
        paste("Device:", input$comp_source_device),
        paste("Vessel:", input$comp_source_vessel),
        paste("Volume:", input$comp_source_volume, "μL"),
        paste("Protocol Type:", input$comp_protocol_type),
        paste("Target Energy:", round(values$comp_source_energy, 0), "J"),
        ""
      )
      
      if(input$comp_input_method == "manual") {
        report_lines <- c(report_lines,
                          "Source Protocol Parameters:",
                          paste("  PIP:", input$comp_manual_pip, "W"),
                          paste("  Duration:", input$comp_manual_duration, "s"),
                          paste("  Duty Factor:", input$comp_manual_duty, "%"),
                          "")
      } else {
        report_lines <- c(report_lines,
                          paste("Target BP:", input$comp_target_bp, "bp"),
                          "")
      }
      
      # add target device info
      report_lines <- c(report_lines,
                        "TARGET DEVICE ANALYSIS",
                        "---------------------",
                        paste("Target Device:", input$comp_target_device),
                        "")
      
      # add vessel options
      if(!is.null(values$comp_vessel_options)) {
        options <- head(values$comp_vessel_options, 5)
        report_lines <- c(report_lines,
                          "Compatible Vessel/Volume Options:",
                          "")
        
        for(i in 1:nrow(options)) {
          report_lines <- c(report_lines,
                            paste0(i, ". ", options[i, 1], " - ", 
                                   options[i, "volume_numeric"], " μL",
                                   " (", ifelse(options[i, "has_protocol_type"], 
                                                "Available", "No data"), ")"))
        }
        
        report_lines <- c(report_lines, "")
        
        # add best recommendation
        best_option <- options[1, ]
        if(best_option$has_protocol_type) {
          result <- translate_protocol(
            values$comp_source_energy,
            input$comp_target_device,
            best_option[, 1],
            best_option$volume_numeric,
            input$comp_protocol_type,
            values$data
          )
          
          if(result$success) {
            report_lines <- c(report_lines,
                              "========================================",
                              "       RECOMMENDED PROTOCOL",
                              "========================================",
                              "",
                              paste("Device:", input$comp_target_device),
                              paste("Vessel:", best_option[, 1]),
                              paste("Volume:", best_option$volume_numeric, "μL"),
                              "",
                              "Protocol Parameters:",
                              paste("  Peak Incident Power:", round(result$pip, 1), "W"),
                              paste("  Duration:", round(result$duration, 1), "s"),
                              paste("  Duty Factor:", round(result$duty_factor, 0), "%"),
                              paste("  Water Level:", result$water_level),
                              paste("  Z-Height:", result$z_height),
                              paste("  Dither:", result$dither),
                              paste("  Calculated Energy:", round(result$calculated_energy, 0), "J"),
                              "",
                              "Validation Metrics:",
                              paste("  Estimated BP:", round(result$estimated_bp, 0), "bp"),
                              paste("  Model R²:", round(result$model_r2, 3)),
                              paste("  Training data points:", result$data_points),
                              "")
            
            if(length(result$warnings) > 0) {
              report_lines <- c(report_lines,
                                "⚠️ Warnings:",
                                paste("  ", result$warnings),
                                "")
            }
          }
        }
      }
      
      report_lines <- c(report_lines,
                        "========================================",
                        "             IMPORTANT NOTES",
                        "========================================",
                        "",
                        "1. Protocol translation assumes similar efficiency across devices",
                        "2. Always validate results experimentally",
                        "3. Volume differences may affect shearing efficiency",
                        "4. Consider vessel-specific factors (geometry, material)",
                        "",
                        "For technical support: bpfeiffer@covaris.com",
                        "")
      
      writeLines(report_lines, file)
    }
  )
  
  # copy to clipboard handler
  observeEvent(input$copy_protocol, {
    req(values$comp_source_energy, values$comp_vessel_options, values$authenticated)
    
    best_option <- values$comp_vessel_options[1, ]
    if(!best_option$has_protocol_type) {
      showNotification("No protocol available to copy", type = "error")
      return()
    }
    
    result <- translate_protocol(
      values$comp_source_energy,
      input$comp_target_device,
      best_option[, 1],
      best_option$volume_numeric,
      input$comp_protocol_type,
      values$data
    )
    
    if(!result$success) {
      showNotification("Protocol translation failed", type = "error")
      return()
    }
    
    # create simple protocol text for clipboard
    protocol_text <- paste(
      "Device:", input$comp_target_device,
      "\nVessel:", best_option[, 1],
      "\nVolume:", best_option$volume_numeric, "μL",
      "\nPIP:", round(result$pip, 1), "W",
      "\nDuration:", round(result$duration, 1), "s",
      "\nDuty Factor:", round(result$duty_factor, 0), "%",
      "\nWater Level:", result$water_level,
      "\nZ-Height:", result$z_height,
      "\nDither:", result$dither,
      "\nEnergy:", round(result$calculated_energy, 0), "J"
    )
    
    # TODO: clipboard functionality requires additional JS
    showNotification("Protocol details displayed - copy manually from the display", type = "info")
  })
}

# run app
shinyApp(ui = ui, server = server)
