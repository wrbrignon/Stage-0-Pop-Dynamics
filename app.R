## Check if shiny package is installed, if not install it.
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
library(popbio)

cv = .2  # coefficent f variation for starting uncertainty ranges
sexr <- 0.50  # sex ratio
s1 <- 0.019  # egg to smolt survival
s2 <- 0.014  # Smolt to 3year old adult survival
s3 <- 0.6  # 3 year old ocean survival
s4 <- 0.7  # 4 year old ocean survival
s5 <- 0.8  # 5 year old ocean survival
b3 <- 0.0  # propensity of age 3 females to spawn
b4 <- 0.8  # propensity of age 4 females to spawn
b5 <- 1 # propensity of age 5 females to spawn
m3 <- 3257  # number of eggs per female spawner of age 3
m4 <- 4095  # number of eggs per female spawner of age 4
m5 <- 5149  # number of eggs per female spawner of age 5

## Define the user interface (UI)
ui <- fluidPage(
  
  # Application Title
  titlePanel("Stage 0 Restoration: Population Dynamics Model"),
  
  fluidRow(
    column(2,
           wellPanel(
             checkboxGroupInput(inputId = "stoch",
                                label = "Stochasticity",
                                choices = list("Demographic" = 1),
                                selected = 0),
             h3("Survival"),
             numericInput(inputId = "s1",
                          label = "s1 = egg to smolt",
                          min = 0, 
                          max = 0.5, 
                          value = s1,
                          step = 0.001),
             numericInput(inputId = "s2",
                          label = "s2 = smolt to age 3 adult",
                          min = 0, 
                          max = 0.5, 
                          value = s2,
                          step = 0.001),
             numericInput(inputId = "s3",
                          label = "s3 = age 3 ocean",
                          min = 0, 
                          max = 1, 
                          value = s3,
                          step = 0.01),
             numericInput(inputId = "s4",
                          label = "s4 = age 4 ocean",
                          min = 0, 
                          max = 1, 
                          value = s4,
                          step = 0.01),
             numericInput(inputId = "s5",
                          label = "s5 = age 5 ocean",
                          min = 0, 
                          max = 1, 
                          value = s5,
                          step = 0.01),
             h3("Other"),
             sliderInput(inputId = "sexr",label = "sexr = Prop. of adults that are female",
                         value = 0.5,
                         min = 0.3,
                         max = 0.7,
                         step = 0.02),
           )),
    column(2,
           wellPanel(
             h3("Spawning probability"),
             numericInput(inputId = "b3",
                          label = "b3 = age 3 female",
                          min = 0, 
                          max = 1, 
                          value = b3,
                          step = 0.01),
             numericInput(inputId = "b4",
                          label = "b4 = age 4 female",
                          min = 0, 
                          max = 1, 
                          value = b4,
                          step = 0.01),
             numericInput(inputId = "b5",
                          label = "b5 = age 5 female",
                          min = 0, 
                          max = 1, 
                          value = b5,
                          step = 0.01),
             
             h3("Eggs/female"),
             sliderInput(inputId = "m3",
                         label = "m3 = age 3",
                         min = 1000, 
                         max =5000, 
                         value = m3,
                         step = 100),
             
             sliderInput(inputId = "m4",
                         label = "m4 = age 4",
                         min = 1000, 
                         max =7500, 
                         value = m4,
                         step = 100),
             
             sliderInput(inputId = "m5",
                         label = "m5 = age 5",
                         min = 1000, 
                         max =10000, 
                         value = m5,
                         step = 100)
           )),
    column(2,wellPanel(
      sliderInput(inputId = "Adults",
                  label = "Age 4 abundance t = 0",
                  min = 0,
                  max = 2000,
                  value = 100,
                  step = 10),
    )),
    column(2,wellPanel(
      sliderInput(inputId = "Years",
                  label = "Number of years",
                  value = 10,
                  min = 1,
                  max = 25,
                  step = 1),
    )),
    column(2,wellPanel(
      numericInput(inputId = "Simz",
                   label = "Number of simulations",
                   value = 1000,
                   min = 100, max = 10000, step = 100),
      )),
    #### Set up the outputs ####
    mainPanel(
      tabsetPanel(type = "tabs",
                  
                  
                  tabPanel("Population Projections",
                           plotOutput("PopnPlot"),
                           tableOutput("PopnTable"),
                           textOutput("Stoch")),
                  tabPanel("Matrix Model", 
                           column(3,
                                  h3("Parameter Values"),
                                  tableOutput("ParameterTable")),
                           
                           h3("Model Structure"),
                           tableOutput("MatrixStructure"), 
                           h3("Model Values"),
                           tableOutput("MatrixValues")),
                  tabPanel("Questions for the group",
                           textOutput("Questions"),
                           h5("Do we need a fry life stage?"),
                           h5("How best to link in the restoration?     Maybe a change in survivals to smolt stage?"),
                           h5("Maybe add densisty dependence in terms of spawning habitat and juvi rearing habitat?"))
      )
    )
  )
)



server <- function(input, output){
  # there are 2 of these, one with byrow = T and one with byrow = F
  matrix.reactive.byrow.T <- reactive(matrix(c(0,0,input$s1*input$b3*input$m3*input$sexr,input$s1*input$b4*input$m4*input$sexr,input$s1*input$b5*input$m5*input$sexr,
                                               input$s2,0,0,0,00,
                                               0,input$s3,0,0,0,
                                               0,0,(1-input$b3)*input$s4,0,0,
                                               0,0,0,(1-input$b4)*input$s5,0),
                                             nrow = 5, byrow = T))
  matrix.reactive.byrow.F <- reactive(matrix(c(0,0,input$s1*input$b3*input$m3*input$sexr,input$s1*input$b4*input$m4*input$sexr,input$s1*input$b5*input$m5*input$sexr,
                                               input$s2,0,0,0,00,
                                               0,input$s3,0,0,0,
                                               0,0,(1-input$b3)*input$s4,0,0,
                                               0,0,0,(1-input$b4)*input$s5,0),
                                             nrow = 5, byrow = F))
  
  output$MatrixStructure <- renderTable({
    A <- as.data.frame(matrix(c(".",".","s1*b3*m3*sexr","s1*b4*m4*sexr","s1*b5*m5*sexr",
                                "s2",".",".",".",".",
                                ".","s3",".",".",".",
                                ".",".","(1-b3)*s4",".",".",
                                ".",".",".","(1-b4)*s5","."),
                              nrow = 5, byrow = T)) 
    colnames(A) <- NULL
    A
  })
  #  matrix model for viewing
  output$MatrixValues <- renderTable({
    matrix.reactive.byrow.T()
  })
  #  parameter table
  output$ParameterTable <- renderTable({
    Parameter <- c("sexr","s1","s2","s3","s4","s5","b3","b4","b5","m3","m4","m5")
    Value <- c(input$sexr,input$s1,input$s2,input$s3,input$s4,input$s5,input$b3,input$b4,input$b5,input$m3,input$m4,input$m5)
    df <- data.frame(Parameter,Value)
    colnames(df) = NULL
    df
  })
  
  
  #  Developes the reactive population projection
  PopProjection <- reactive({
    # setup starting populations using teh SAD from adult 4 year old returns
    sad <- stable.stage(matrix.reactive.byrow.T())
    age5 <- input$Adults/sad[4]*sad[5]
    age3 <- input$Adults/sad[4]*sad[3]
    age2 <- input$Adults/sad[4]*sad[2]
    age1 <- input$Adults/sad[4]*sad[1]
    outs <- array(NA,dim = c(5,input$Years+1,3))
    popn <- matrix(NA,nrow = 5,ncol = input$Years+1, byrow = T)
    colnames(popn) = seq(0,input$Years,by = 1)
    popn[,1] <- round(c(age1,age2,age3,input$Adults,age5))
    
    ## Deterministic
    for (i in 1:input$Years){
      popn[,i+1] = round(popn[,i] %*% matrix.reactive.byrow.F(),digits = 0)
    }
    
    ## Demographic stochasticity
    if ("1" %in% input$stoch & length(input$stoch) == 1){
      popn.array <- array(NA,dim = c(5,input$Years+1,input$Simz))
      popn.array[,1,] <- round(c(age1,age2,age3,input$Adults,age5))
      
      ## Function for demographic stochasiticy
      '%ds%' <- function (popn,TM){
        N1 <- sum(rpois(length(popn),popn*TM[1,]))
        N2 <- sum(rbinom(length(popn),popn,TM[2,]))
        N3 <- sum(rbinom(length(popn),popn,TM[3,]))
        N4 <- sum(rbinom(length(popn),popn,TM[4,]))
        N5 <- sum(rbinom(length(popn),popn,TM[5,]))
        c(N1,N2,N3,N4,N5)
      }  
      
      for (j in 1:input$Simz){
        for (i in 1:input$Years){
          popn.array[,i+1,j] = popn.array[,i,j] %ds% matrix.reactive.byrow.T()
        }
      }
      popn <- apply(popn.array,c(1,2),mean)  ## mean popn
      high.popn <- apply(popn.array,c(1,2),max)  ## max popn of simz
      low.popn <- apply(popn.array,c(1,2),min)  ## Lower popn of simz
      outs[,,2] = high.popn
      outs[,,3] = low.popn
    }
    outs[,,1] = popn
    
    outs
  })
  
  #  Population projection plot codez
  output$PopnPlot <- renderPlot({
    spawning.abundance <- round((PopProjection()[3,,1] * input$b3) + (PopProjection()[4,,1] * input$b4) + (PopProjection()[5,,1] * input$b5),digits = 0)
    plot(spawning.abundance,type = "l", ylab = "Adult Abundance", xlab = "Year")
    UCL.spawning.abundance <- round((PopProjection()[3,,2] * input$b3) + (PopProjection()[4,,2] * input$b4) + (PopProjection()[5,,2] * input$b5),digits = 0)
    LCL.spawning.abundance <- round((PopProjection()[3,,3] * input$b3) + (PopProjection()[4,,3] * input$b4) + (PopProjection()[5,,3] * input$b5),digits = 0)
    lines(UCL.spawning.abundance,type = "l",col = "blue")
    lines(LCL.spawning.abundance,type = "l",col = "blue")
  })
  #  Population projection table codez
  output$PopnTable <- renderTable(PopProjection()[,,1])

}



shinyApp(ui = ui, server = server)
# 

# runApp("Stage 0 pop dy shiny.R")
