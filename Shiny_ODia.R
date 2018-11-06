library("deSolve")

require(ReacTran, quietly = TRUE) # Reaction-Transport models in aquatic systems 
require(marelac, quietly = TRUE)  # toolbox for aquatic sciences

grid <- setup.grid.1D(N = 100, L = 10, dx.1 = 0.01)

O2model <- function (t, O2, p) {
  with (as.list(p), {
    # Transport term (using ReacTran routines)
    O2tran <- tran.1D(C = O2, C.up = O2BW, D = D/(1-log(porosity**2)), VF = porosity, dx = grid)
    # Respiration
    O2cons <- minrate*(grid$x.mid < mindepth)    
    # the function returns the time derivative
    return( list(O2tran$dC - O2cons,
                 O2cons = O2cons,
                 O2flux = O2tran$flux.up,
                 Resp= sum(O2cons*porosity*grid$dx) ) ) 
  })
}

parms <- c(
  porosity = 0.8    , # -
  minrate  = 10    , # nmol O2/cm3/d - oxygen consumption rate
  mindepth = 5    , # cm            - depth below which minrate = 0
  O2BW     = 300  , # nmol/cm3       - bottom water oxygen concentration
  D        = as.numeric(diffcoeff(species="O2")*86400*1e4)     # cm2/d - molecular diffusion coefficient
)

IC <- rep.int(c(200),length(grid$x.mid))

# computes the steady-state solution
DefaultRun  <- steady.1D(y = IC, parms = parms, func = O2model, nspec = 1, names = "O2")

server <- function(input, output,session) {
  observeEvent(input$resetButton, {
    updateNumericInput(session, "porosity", value = 0.8)
    updateNumericInput(session, "minrate" , value = 10)
    updateNumericInput(session, "mindepth", value = 5)
    updateNumericInput(session, "O2BW"    , value = 300)
    updateNumericInput(session, "D"       , value = as.numeric(diffcoeff(species="O2")*86400*1e4))
  })
  
  output$model <- renderPlot({
    IC <- rep.int(200,length(grid$x.mid))
    Parms <- c(porosity=input$porosity, minrate=input$minrate,
                 mindepth=input$mindepth, O2BW=input$O2BW, D = input$D)
    # computes the steady-state solution
    out  <- steady.1D(y = IC, parms = Parms, func = O2model, nspec = 1, names = "O2")
    
    if (input$default) {
      plot(out, DefaultRun, xyswap = TRUE, xlab = "mmol/m3", ylab = "cm", grid = grid$x.mid)
      legend("bottomright", col = 1:2, legend = c("current", "default"), lty = 1:2)
    } else  
      plot(out, xyswap = TRUE, xlab = "mmol/m3", ylab = "cm", grid = grid$x.mid)
  })
  
  output$table <- renderTable({
    IC <- rep.int(200,length(grid$x.mid))
    Parms <- c(porosity=input$porosity, minrate=input$minrate,
               mindepth=input$mindepth, O2BW=input$O2BW, D = input$D)
    # computes the steady-state solution
    out  <- steady.1D(y = IC, parms = Parms, func = O2model, nspec = 1, names = "O2")
    
    if (input$default) {
      data.frame("Current"=c(out$O2flux,out$Resp),
                 "Default"=c(DefaultRun$O2flux,DefaultRun$Resp),
                 row.names = c("O2 Flux [nmol/cm2/d]",
                               "Total Resp [nmol/cm2/d]"))
    } else
      data.frame("Current"=c(out$O2flux,out$Resp),
                 row.names = c("O2 Flux [nmol/cm2/d]",
                               "Total Resp [nmol/cm2/d]"))
    },include.rownames=T)
}

ui <- fluidPage(
  headerPanel("O2 model"),
  sidebarLayout(
    sidebarPanel(
      h3("Parameters"),
      sliderInput("porosity", label = "porosity",
                   min = 0.3, max = 0.99,  value = 0.8, step = 0.05, width=100),
      sliderInput("minrate", label = "minrate",
                   min = 0.0, max = 100,  value = parms["minrate"], step = 0.5, width=100),
      sliderInput("mindepth", label = "mindepth",
                   min = 0.1, max = 10, value = parms["mindepth"], step = 0.1, width=100),
      sliderInput("O2BW", label = "O2BW",
                  min = 0, max = 600, value = parms["O2BW"], step = 10, width=100),
      sliderInput("D", label = "D",
                  min = 0.5, max = 2.5, value = parms["D"], step = 0.05, width=100),
      checkboxInput(inputId = "default",
                    label = strong("Add default run"),
                    value = TRUE),
      actionButton("resetButton",label="Reset Parameters"),
      br()   # ends without ','
      ),
    mainPanel(
      h3("Simulation results"),
      plotOutput("model"),
      tableOutput("table")
    )
  )
)

shinyApp(ui = ui, server = server)