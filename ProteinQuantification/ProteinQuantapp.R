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
ui <- fluidPage(titlePanel("BCA Protein Analysis & Western Blot Prep"),
                sidebarLayout(sidebarPanel(
                  fileInput('file1', 'Choose xlsx file',
                            accept = c(".xlsx", ".csv")),
                  numericInput("protein",
                               "Desired Amount of Protein per Alliquot (\u03BCg)",
                               min = 1,
                               max = 10000,
                               value = 75),
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
                               "Amount of Protein Desired for Western Blot (\u03BCg)",
                               min = 0,
                               max = 10000,
                               value = 90),
                  numericInput("wb",
                               "Volume Desired for Western Blot (\u03BCL)",
                               min = 0,
                               max = 10000,
                               value = 90),
                  checkboxInput("ripabox", "Check if set volume of RIPA is needed"),
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
                    tabPanel( "Sample Names",
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
                )

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
    
    thedata <- reactive(bca_table)
    
    
    output$download_bca <- downloadHandler(
      filename = function() {paste("BCA_analysis_",Sys.Date(),".csv", sep = "")}, 
      content = function(fname){
        write.csv(thedata(), fname, col.names = TRUE, row.names = FALSE)
      }
    )
    
    bca_table
    
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
  
  output$western <- renderTable({
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
    
    wb_table <- data.frame(t(content_table))
    
    wb_table$vol_sample <- (input$wb/wb_table)
    
    colnames(wb_table) <- c("Sample Protein Concentration", "Sample Volume")
    
    wb_table$`Sample ID` <- c(paste(input$sample1),paste(input$sample2),paste(input$sample3),paste(input$sample4),paste(input$sample5),
      paste(input$sample6),paste(input$sample7),paste(input$sample8),paste(input$sample9),paste(input$sample10),
      paste(input$sample11),paste(input$sample12),paste(input$sample13),paste(input$sample14),paste(input$sample15),
      paste(input$sample16),paste(input$sample17),paste(input$sample18))
    
    wb_table <- wb_table |> dplyr::select('Sample ID', everything())
    
    wb_table$`Buffer Volume` <- (max(wb_table$`Sample Volume`) - wb_table$`Sample Volume`)/(input$wbpro/input$wb)
    
    wb_table$`Total Volume` <- input$wb
    
    for(k in input$ripabox){
      
      if(k == TRUE){
        
        wb_table$`RIPA Volume` <- input$ripa
      
      for(i in input$lsb){
        if(i == "2X"){
          lsb_vol <- wb_table$`Sample Volume`
        }
        
        if(i == "4X"){
          lsb_vol <- ((wb_table$`Sample Volume` + wb_table$`Buffer Volume` + wb_table$`RIPA Volume`)/3)
        }
        
        if(i == "6X"){
          lsb_vol <- ((wb_table$`Sample Volume` + wb_table$`Buffer Volume` + wb_table$`RIPA Volume`)/5)
        }
      }
      
      }
      
      
      if(k == FALSE){
        for(i in input$lsb){
          if(i == "2X"){
            wb_table$`RIPA Volume` <- (input$wb - wb_table$`Sample Volume` + wb_table$`Buffer Volume`)/2
            
            lsb_vol <- input$wb/2
            
         }
        
          if(i == "4X"){
            wb_table$`RIPA Volume` <- ((3/4)*input$wb) - wb_table$`Sample Volume` - wb_table$`Buffer Volume`
            
            lsb_vol <- (1/4)*input$wb
            
          }
        
          if(i == "6X"){
            wb_table$`RIPA Volume` <- ((5/6)*input$wb) - wb_table$`Sample Volume` - wb_table$`Buffer Volume`
            
            lsb_vol <- (1/6)*input$wb
          }
        
          
        
        }
      
    }
    }
        wb_table$`Laemmli Buffer Volume` <- lsb_vol

 
        wb_table$`Total Volume` <- wb_table$`Sample Volume` + wb_table$`Buffer Volume` + wb_table$`RIPA Volume` + wb_table$`Laemmli Buffer Volume`

        wb_table$`Sample Protein Concentration` <- replace(wb_table$`Sample Protein Concentration`, which(wb_table$`Sample Protein Concentration` < 0), NA)
        

    na.omit(wb_table)
    
  }, hover = TRUE, striped = T, bordered = T, align = "c", caption = "All volumes are in \u03BCL")
  
  
}


# Run the application 
shinyApp(ui = ui, server = server)

