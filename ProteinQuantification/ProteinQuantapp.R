#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(navbarPage("BCA Protein Analysis & Western Blot Prep",
                tabPanel("Analysis", fluid = TRUE, icon = icon('calculator'),
                sidebarLayout(sidebarPanel(
                  fileInput('file1', 'Choose xlsx file',
                            accept = c(".xlsx", ".csv")),
                  numericInput("protein",
                               "Desired Amount of Protein per Alliquot (\u03BCg)",
                               min = 1,
                               max = 10000,
                               value = 75),
                  numericInput("protein2",
                               "Second Desired Amount of Protein per Alliquot (\u03BCg)",
                               min = 1,
                               max = 10000,
                               value = 5),
                  numericInput("protein3",
                               "ThirdDesired Amount of Protein per Alliquot (\u03BCg)",
                               min = 1,
                               max = 10000,
                               value = 100),
                  numericInput("std_min",
                               "Minimum Standard Concentration",
                               min = 0,
                               max = 10000,
                               value = 0),
                  numericInput("std_max",
                               "Maximum Standard Concentration",
                               min = 0,
                               max = 10000,
                               value = 10),
                  numericInput("wells",
                               "Number of Standards",
                               min = 0,
                               max = 10000,
                               value = 6),
                  numericInput("steps",
                               "Steps",
                               min = 0,
                               max = 10000,
                               value = 2),
                  numericInput("wbpro",
                               "Protein in Final WB Preparation Solution (\u03BCg)",
                               min = 0,
                               max = 10000,
                               value = 90),
                  numericInput("wb",
                               "Volume of Final WB Preparation Solution (\u03BCL)",
                               min = 0,
                               max = 10000,
                               value = 90),
                  checkboxInput("ripabox", "Check if no added buffer is needed"),
                  numericInput("ripa",
                               "Volume Desired for RIPA (\u03BCL)",
                               min = 0,
                               max = 10000,
                               value = 5),
                  selectInput("lsb",
                               "Laemmli Buffer Concentration",
                               choices = c("Select","2X", "4X", "6X"))),
                  
                mainPanel(
                  tabsetPanel(
                    tabPanel("Volume of Sample for Desired Amount of Protein (\u03BCL)", tableOutput('calculated')),
                    tabPanel("Western Blot", tableOutput('western')),
                    tabPanel("Raw Data", tableOutput('raw')),
                    tabPanel("Sample Protein Concentration (\u03BCg/mL)", tableOutput('content')),
                    tabPanel("Linear Model", verbatimTextOutput('model')),
                    tabPanel("Sample Names",
                      textInput("sample1", "Sample 1 Name"),
                      textInput("sample2", "Sample 2 Name"),
                      textInput("sample3", "Sample 3 Name"),
                      textInput("sample4", "Sample 4 Name"),
                      textInput("sample5", "Sample 5 Name"),
                      textInput("sample6", "Sample 6 Name"),
                      textInput("sample7", "Sample 7 Name"),
                      textInput("sample8", "Sample 8 Name"),
                      textInput("sample9", "Sample 9 Name"),
                      textInput("sample10", "Sample 10 Name"),
                      textInput("sample11", "Sample 11 Name"),
                      textInput("sample12", "Sample 12 Name"),
                      textInput("sample13", "Sample 13 Name"),
                      textInput("sample14", "Sample 14 Name"),
                      textInput("sample15", "Sample 15 Name"),
                      textInput("sample16", "Sample 16 Name"),
                      textInput("sample17", "Sample 17 Name"),
                      textInput("sample18", "Sample 18 Name"))),
                  downloadButton('download_bca',"Download BCA Data"),
                  downloadButton('download_wb',"Download Western Blot Data"),
                  )
                )
                ),
                tabPanel("User Guide", fluid = TRUE, icon = icon('book-open'), 
                         fluidRow(
                  h2(HTML("<b>Instructions</b>")),
                  h5(HTML("<p>For version control and other questions, see my <a href='https://github.com/stdecker/ProteinQuantification-BCA-Western'>GitHub repo</a>.</p>")),
                  h4(HTML("<b>Uploading Data</b>")),
                  h5(HTML("Note: <b>Before uploading, data must be in the format described on my GitHub Repo (see above).</b>")),
                  h5("Upload data using the 'Choose xlsx File'"),
                  h4(HTML("<b>Adjusting Parameters</b>")),
                  h5("Adjust necessary parameters using the input selections given on the left of the page."),
                  h4(HTML("<b>Naming Samples</b>")),
                  h5("Samples can be named using the 'Sample Names' tab and typing in the names of each corresponding sample."),
                  h4(HTML("<b>Saving and Downloading Data</b>")),
                  h5("After data have been checked and verified, select the 'Download Final Data' button. This will prompt a new window where you can rename the CSV file and save it onto your local system. 
                     Hit save. The data should be stored in the selected location and the app may now be closed.")
                  )
)))

server <- function(input, output) {
  output$raw <- renderTable({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    table1 <- readxl::read_excel(inFile$datapath, col_names = T)
    colnames(table1) <- c("", "Well 1", "Well 2", "Well 3", "Well 4", "Well 5", "Well 6", "Well 7", "Well 8", "Well 9", "Well 10", "Well 11", "Well 12")                             
    table1
  }, hover = TRUE, striped = T, bordered = T, align = "c", caption = "Raw Values")
  
  output$content <- renderTable({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    table1 <- data.frame(readxl::read_excel(inFile$datapath,sheet=1, col_names = T))
    
    
    stdmeans_table <- data.frame(rowMeans(table1[2:4], na.rm = TRUE))
    
    std_table <- data.frame(stdmeans_table[c(1:input$wells),])
    
    colnames(std_table) <- c("Means")
    
    std_table$BSA <- seq(input$std_min, input$std_max, input$steps)
    
    calculated_table <- table1[-c(1)] - min(std_table[1])
    
    std_table[1] <- std_table[1] - min(std_table[1])
    
    lin_model <- lm(std_table$Means ~ std_table$BSA)
    
    Sample1 <- (mean(c(calculated_table[1,4], calculated_table[2,4], calculated_table[3,4]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample2 <- (mean(c(calculated_table[1,5], calculated_table[2,5], calculated_table[3,5]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample3 <- (mean(c(calculated_table[1,6], calculated_table[2,6], calculated_table[3,6]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample4 <- (mean(c(calculated_table[1,7], calculated_table[2,7], calculated_table[3,7]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample5 <- (mean(c(calculated_table[1,8], calculated_table[2,8], calculated_table[3,8]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample6 <- (mean(c(calculated_table[1,9], calculated_table[2,9], calculated_table[3,9]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample7 <- (mean(c(calculated_table[1,10], calculated_table[2,10], calculated_table[3,10]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample8 <- (mean(c(calculated_table[1,11], calculated_table[2,11], calculated_table[3,11]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample9 <- (mean(c(calculated_table[1,12], calculated_table[2,12], calculated_table[3,12]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample10 <- (mean(c(calculated_table[7,4], calculated_table[5,4], calculated_table[6,4]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample11 <- (mean(c(calculated_table[7,5], calculated_table[5,5], calculated_table[6,5]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample12 <- (mean(c(calculated_table[7,6], calculated_table[5,6], calculated_table[6,6]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample13 <- (mean(c(calculated_table[7,7], calculated_table[5,7], calculated_table[6,7]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample14 <- (mean(c(calculated_table[7,8], calculated_table[5,8], calculated_table[6,8]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample15 <- (mean(c(calculated_table[7,9], calculated_table[5,9], calculated_table[6,9]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample16 <- (mean(c(calculated_table[7,10], calculated_table[5,10], calculated_table[6,10]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample17 <- (mean(c(calculated_table[7,11], calculated_table[5,11], calculated_table[6,11]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample18 <- (mean(c(calculated_table[7,12], calculated_table[5,12], calculated_table[6,12]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    
    content_table <- cbind(Sample1, Sample2, Sample3, Sample4, Sample5, Sample6, Sample7, Sample8, Sample9, Sample10, Sample11, Sample12, Sample13, Sample14, Sample15, Sample16, Sample17, Sample18)
    
    colnames(content_table) <- c(paste(input$sample1),paste(input$sample2),paste(input$sample3),paste(input$sample4),paste(input$sample5),
                                 paste(input$sample6),paste(input$sample7),paste(input$sample8),paste(input$sample9),paste(input$sample10),
                                 paste(input$sample11),paste(input$sample12),paste(input$sample13),paste(input$sample14),paste(input$sample15),
                                 paste(input$sample16),paste(input$sample17),paste(input$sample18))
    
    content_table <- replace(content_table, which(content_table < 0), NA)
    
    t(content_table[, colSums(is.na(content_table)) != nrow(content_table)])
    
  }, hover = TRUE, striped = T, bordered = T, align = "c", caption = "Protein Content")
  
  
  output$calculated <- renderTable({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    table1 <- data.frame(readxl::read_excel(inFile$datapath,sheet=1, col_names = T))
    
    
    stdmeans_table <- data.frame(rowMeans(table1[2:4], na.rm = TRUE))
    
    std_table <- data.frame(stdmeans_table[c(1:input$wells),])
    
    colnames(std_table) <- c("Means")
    
    std_table$BSA <- seq(input$std_min, input$std_max, input$steps)
    
    calculated_table <- table1[-c(1)] - min(std_table[1])
    
    std_table[1] <- std_table[1] - min(std_table[1])
    
    lin_model <- lm(std_table$Means ~ std_table$BSA)
    
    Sample1 <- (mean(c(calculated_table[1,4], calculated_table[2,4], calculated_table[3,4]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample2 <- (mean(c(calculated_table[1,5], calculated_table[2,5], calculated_table[3,5]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample3 <- (mean(c(calculated_table[1,6], calculated_table[2,6], calculated_table[3,6]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample4 <- (mean(c(calculated_table[1,7], calculated_table[2,7], calculated_table[3,7]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample5 <- (mean(c(calculated_table[1,8], calculated_table[2,8], calculated_table[3,8]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample6 <- (mean(c(calculated_table[1,9], calculated_table[2,9], calculated_table[3,9]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample7 <- (mean(c(calculated_table[1,10], calculated_table[2,10], calculated_table[3,10]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample8 <- (mean(c(calculated_table[1,11], calculated_table[2,11], calculated_table[3,11]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample9 <- (mean(c(calculated_table[1,12], calculated_table[2,12], calculated_table[3,12]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample10 <- (mean(c(calculated_table[7,4], calculated_table[5,4], calculated_table[6,4]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample11 <- (mean(c(calculated_table[7,5], calculated_table[5,5], calculated_table[6,5]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample12 <- (mean(c(calculated_table[7,6], calculated_table[5,6], calculated_table[6,6]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample13 <- (mean(c(calculated_table[7,7], calculated_table[5,7], calculated_table[6,7]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample14 <- (mean(c(calculated_table[7,8], calculated_table[5,8], calculated_table[6,8]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample15 <- (mean(c(calculated_table[7,9], calculated_table[5,9], calculated_table[6,9]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample16 <- (mean(c(calculated_table[7,10], calculated_table[5,10], calculated_table[6,10]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample17 <- (mean(c(calculated_table[7,11], calculated_table[5,11], calculated_table[6,11]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample18 <- (mean(c(calculated_table[7,12], calculated_table[5,12], calculated_table[6,12]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    
    content_table <- cbind(Sample1, Sample2, Sample3, Sample4, Sample5, Sample6, Sample7, Sample8, Sample9, Sample10, Sample11, Sample12, Sample13, Sample14, Sample15, Sample16, Sample17, Sample18)
    
    
    bca_table <- 1/(content_table/input$protein)
    
    colnames(bca_table) <- c(paste(input$sample1),paste(input$sample2),paste(input$sample3),paste(input$sample4),paste(input$sample5),
                                 paste(input$sample6),paste(input$sample7),paste(input$sample8),paste(input$sample9),paste(input$sample10),
                                 paste(input$sample11),paste(input$sample12),paste(input$sample13),paste(input$sample14),paste(input$sample15),
                                 paste(input$sample16),paste(input$sample17),paste(input$sample18))
    
    bca_table <- replace(bca_table, which(bca_table < 0), NA)
    
    bca_table <- t(bca_table[, colSums(is.na(bca_table)) != nrow(bca_table)])
    
    bca1 <- merge(data.frame(Amount = paste(input$protein, "\u03BCg")), bca_table)
    
    bca_table2 <- 1/(content_table/input$protein2)
    
    colnames(bca_table2) <- c(paste(input$sample1),paste(input$sample2),paste(input$sample3),paste(input$sample4),paste(input$sample5),
                             paste(input$sample6),paste(input$sample7),paste(input$sample8),paste(input$sample9),paste(input$sample10),
                             paste(input$sample11),paste(input$sample12),paste(input$sample13),paste(input$sample14),paste(input$sample15),
                             paste(input$sample16),paste(input$sample17),paste(input$sample18))
    
    bca_table2 <- replace(bca_table2, which(bca_table2 < 0), NA)
    
    bca_table2 <- t(bca_table2[, colSums(is.na(bca_table2)) != nrow(bca_table2)])
    
    bca2 <- merge(data.frame(Amount = paste(input$protein2, "\u03BCg")), bca_table2)
    
    bca_table3 <- 1/(content_table/input$protein3)
    
    colnames(bca_table3) <- c(paste(input$sample1),paste(input$sample2),paste(input$sample3),paste(input$sample4),paste(input$sample5),
                              paste(input$sample6),paste(input$sample7),paste(input$sample8),paste(input$sample9),paste(input$sample10),
                              paste(input$sample11),paste(input$sample12),paste(input$sample13),paste(input$sample14),paste(input$sample15),
                              paste(input$sample16),paste(input$sample17),paste(input$sample18))
    
    bca_table3 <- replace(bca_table3, which(bca_table3 < 0), NA)
    
    bca_table3 <- t(bca_table3[, colSums(is.na(bca_table3)) != nrow(bca_table3)])
    
    bca3 <- merge(data.frame(Amount = paste(input$protein3, "\u03BCg")), bca_table3)
    
    bca_table_comb <- rbind(bca1, bca2, bca3)
    
    thedata <- reactive(bca_table_comb)
    
    
    output$download_bca <- downloadHandler(
      filename = function() {paste("BCA_analysis_",Sys.Date(),".csv", sep = "")}, 
      content = function(fname){
        write.csv(thedata(), fname, col.names = TRUE, row.names = FALSE)
      }
    )
    
    bca_table_comb
    
  }, hover = TRUE, striped = T, bordered = T, align = "c", caption = "Volume for Selected Protein Content (\u03BCL)")
  
  
  output$model <- renderPrint({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    
    table1 <- data.frame(readxl::read_excel(inFile$datapath,sheet=1, col_names = T))
    
    
    stdmeans_table <- data.frame(rowMeans(table1[2:4], na.rm = TRUE))
    
    std_table <- data.frame(stdmeans_table[c(1:input$wells),])
    
    colnames(std_table) <- c("Means")
    
    std_table$BSA <- seq(input$std_min, input$std_max, input$steps)
    
    calculated_table <- table1[-c(1)] - min(std_table[1])
    
    summary(lm(std_table$Means ~ std_table$BSA))})
  
  wb_data <- reactive({ 
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    table1 <- data.frame(readxl::read_excel(inFile$datapath,sheet=1, col_names = T))
    
    
    stdmeans_table <- data.frame(rowMeans(table1[2:4], na.rm = TRUE))
    
    std_table <- data.frame(stdmeans_table[c(1:input$wells),])
    
    colnames(std_table) <- c("Means")
    
    std_table$BSA <- seq(input$std_min, input$std_max, input$steps)
    
    calculated_table <- table1[-c(1)] - min(std_table[1])
    
    std_table[1] <- std_table[1] - min(std_table[1])
    
    lin_model <- lm(std_table$Means ~ std_table$BSA)
    
    Sample1 <- (mean(c(calculated_table[1,4], calculated_table[2,4], calculated_table[3,4]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample2 <- (mean(c(calculated_table[1,5], calculated_table[2,5], calculated_table[3,5]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample3 <- (mean(c(calculated_table[1,6], calculated_table[2,6], calculated_table[3,6]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample4 <- (mean(c(calculated_table[1,7], calculated_table[2,7], calculated_table[3,7]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample5 <- (mean(c(calculated_table[1,8], calculated_table[2,8], calculated_table[3,8]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample6 <- (mean(c(calculated_table[1,9], calculated_table[2,9], calculated_table[3,9]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample7 <- (mean(c(calculated_table[1,10], calculated_table[2,10], calculated_table[3,10]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample8 <- (mean(c(calculated_table[1,11], calculated_table[2,11], calculated_table[3,11]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample9 <- (mean(c(calculated_table[1,12], calculated_table[2,12], calculated_table[3,12]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample10 <- (mean(c(calculated_table[7,4], calculated_table[5,4], calculated_table[6,4]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample11 <- (mean(c(calculated_table[7,5], calculated_table[5,5], calculated_table[6,5]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample12 <- (mean(c(calculated_table[7,6], calculated_table[5,6], calculated_table[6,6]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample13 <- (mean(c(calculated_table[7,7], calculated_table[5,7], calculated_table[6,7]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample14 <- (mean(c(calculated_table[7,8], calculated_table[5,8], calculated_table[6,8]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample15 <- (mean(c(calculated_table[7,9], calculated_table[5,9], calculated_table[6,9]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample16 <- (mean(c(calculated_table[7,10], calculated_table[5,10], calculated_table[6,10]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample17 <- (mean(c(calculated_table[7,11], calculated_table[5,11], calculated_table[6,11]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    Sample18 <- (mean(c(calculated_table[7,12], calculated_table[5,12], calculated_table[6,12]), na.rm = T) - lin_model$coefficients[1])/lin_model$coefficients[2]
    
    content_table <- cbind(Sample1, Sample2, Sample3, Sample4, Sample5, Sample6, Sample7, Sample8, Sample9, Sample10, Sample11, Sample12, Sample13, Sample14, Sample15, Sample16, Sample17, Sample18)
    
    colnames(content_table) <- c(paste(input$sample1),paste(input$sample2),paste(input$sample3),paste(input$sample4),paste(input$sample5),
                                 paste(input$sample6),paste(input$sample7),paste(input$sample8),paste(input$sample9),paste(input$sample10),
                                 paste(input$sample11),paste(input$sample12),paste(input$sample13),paste(input$sample14),paste(input$sample15),
                                 paste(input$sample16),paste(input$sample17),paste(input$sample18))
    
    pro_concentrarion <- t(content_table)
    
    vol_sample <- (input$wbpro/pro_concentrarion)
    
    wb_table <- data.frame(pro_concentrarion, vol_sample)
    
    colnames(wb_table) <- c("Sample Protein Concentration", "Sample Volume")
    
    wb_table$`Sample ID` <- c(paste(input$sample1),paste(input$sample2),paste(input$sample3),paste(input$sample4),paste(input$sample5),
      paste(input$sample6),paste(input$sample7),paste(input$sample8),paste(input$sample9),paste(input$sample10),
      paste(input$sample11),paste(input$sample12),paste(input$sample13),paste(input$sample14),paste(input$sample15),
      paste(input$sample16),paste(input$sample17),paste(input$sample18))
    
    wb_table <- wb_table |> dplyr::select('Sample ID', everything())
    
    for(k in input$ripabox){
      
      if(k == FALSE){
        
        wb_table$`RIPA Volume` <- input$ripa
        
      for(i in input$lsb){
        if(i == "2X"){
          wb_table$`Buffer Volume` <- input$wb/(2/1) - wb_table$`Sample Volume` - wb_table$`RIPA Volume`
          lsb_vol <- input$wb/2
        }
        
        if(i == "4X"){
          wb_table$`Buffer Volume` <- input$wb/(4/3) - wb_table$`Sample Volume` - wb_table$`RIPA Volume`
          lsb_vol <- input$wb/4
        }
        
        if(i == "6X"){
          wb_table$`Buffer Volume` <- input$wb/(6/5) - wb_table$`Sample Volume` - wb_table$`RIPA Volume`
          lsb_vol <- input$wb/6
        }
      }
      
      }
      
      
      if(k == TRUE){
        for(i in input$lsb){
          if(i == "2X"){
            wb_table$`RIPA Volume` <- 0
            
            wb_table$`Buffer Volume` <- (input$wb/2) - wb_table$`Sample Volume`
            
            lsb_vol <- input$wb/2
            
         }
        
          if(i == "4X"){
            wb_table$`RIPA Volume` <- 0
            
            wb_table$`Buffer Volume` <- (input$wb*(3/4)) - wb_table$`Sample Volume`
            
            
            lsb_vol <- input$wb/4
            
          }
        
          if(i == "6X"){
            wb_table$`RIPA Volume` <- 0
            
            wb_table$`Buffer Volume` <- (input$wb*(5/6)) - wb_table$`Sample Volume`
            
            lsb_vol <- input$wb/6
          }
        
          
        
        }
      
    }
    }
        wb_table$`Laemmli Buffer Volume` <- lsb_vol
        
        wb_table$`Final Protein Concentration` <- input$wbpro/input$wb
 
        wb_table$`Total Volume` <- wb_table$`Sample Volume` + wb_table$`Buffer Volume` + wb_table$`RIPA Volume` + wb_table$`Laemmli Buffer Volume`

        wb_table$`Sample Protein Concentration` <- replace(wb_table$`Sample Protein Concentration`, which(wb_table$`Sample Protein Concentration` < 0), NA)
        
        colnames(wb_table) <- c('Sample ID', paste0('Sample Protein Concentration (μg', '/', 'μL)'), 'Sample Volume (μL)', 'Added Buffer Volume (μL)', 'Buffer Volume (μL)', 'Laemmli Buffer Volume (μL)',  paste0('Final Protein Concentration (μg', '/', 'μL)'), 'Total Volume (μL)')
        

   wb_table <- na.omit(wb_table)
   return(wb_table)
   


  })
  
  
  output$western <- renderTable({
    
    wb_data()
    
  }, hover = TRUE, striped = T, bordered = T, align = "c", caption = "Added buffer examples: RIPA, NaCl, etc.")
  
  
  output$download_wb <- downloadHandler(
    filename = function() {paste("WB_analysis_",Sys.Date(),".csv", sep = "")}, 
    content = function(fname){
      write.csv(wb_data(), fname, col.names = TRUE, row.names = FALSE)
    }
  )
  
}


# Run the application 
shinyApp(ui = ui, server = server)


