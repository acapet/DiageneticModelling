library("deSolve")

require(ReacTran, quietly = TRUE) # Reaction-Transport models in aquatic systems 
require(marelac, quietly = TRUE)  # toolbox for aquatic sciences

grid <- setup.grid.1D(N = 100, L = 10, dx.1 = 0.01)

O2model <- function (t, O2, p) {
  with (as.list(p), {
    # Transport term (using ReacTran routines)
    O2tran <- tran.1D(C = O2, C.up = O2BW, D = D, VF = porosity, dx = grid)
    # Respiration
    O2cons <- minrate*(O2/(O2+ks))
    # the function returns the time derivative
    return( list(O2tran$dC - O2cons,
                 O2cons = O2cons,
                 O2flux = O2tran$flux.up,
                 Resp= sum(O2cons*porosity*grid$dx) ) ) 
  })
}

parms <- c(
  porosity = 0.8    , # -
  minrate  = 1000    , # mmol O2/m3/d - oxygen consumption rate
  ks       = 1      , # mmol O2/m3 - Half saturation constant
  O2BW     = 300  , # mmol/m3       - bottom water oxygen concentration
  D        = as.numeric(diffcoeff(species="O2", S=35.7, t=24)*86400*1e4)/(1-log(0.8*0.8))     # cm2/d - molecular diffusion coefficient
)

IC <- rep.int(c(200),length(grid$x.mid))

# computes the steady-state solution
DefaultRun  <- steady.1D(y = IC, parms = parms, func = O2model, nspec = 1, names = "O2", pos=TRUE)

server <- function(input, output,session) {
  observeEvent(input$resetButton, {
    updateNumericInput(session, "porosity", value = 0.8)
    updateNumericInput(session, "minrate", value = 1000)
    updateNumericInput(session, "ks", value = 1)
    updateNumericInput(session, "O2BW", value = 300)
  })
  
  output$model <- renderPlot({
    IC <- rep.int(200,length(grid$x.mid))
    Parms <- c(porosity=input$porosity,
               minrate=input$minrate,
               ks=input$ks,
               O2BW=input$O2BW,
               D = parms[["D"]])
    # computes the steady-state solution
    out  <- steady.1D(y = IC, parms = Parms, func = O2model, nspec = 1, names = "O2", pos=TRUE)
    
    if (input$default) {
      plot(out, DefaultRun, xyswap = TRUE, xlab = "mmol/m3", ylab = "cm",
           grid = grid$x.mid, ylim=c(2,0))
      legend("bottomright", col = 1:2, legend = c("current", "default"), lty = 1:2)
    } else{  
      plot(out, xyswap = TRUE, xlab = "mmol/m3", ylab = "cm",
           grid = grid$x.mid, ylim=c(2,0))
    }
  })
  
  output$table <- renderTable({
    IC <- rep.int(200,length(grid$x.mid))
    Parms <- c(porosity=input$porosity,
               minrate=input$minrate,
               ks=input$ks,
               O2BW=input$O2BW,
               D = parms[["D"]])
    # computes the steady-state solution
    out  <- steady.1D(y = IC, parms = Parms, func = O2model, nspec = 1, names = "O2", pos=TRUE)
    
    if (input$default) {
      data.frame("Current"=c(out$O2flux/100,out$Resp/100),
                 "Default"=c(DefaultRun$O2flux/100,DefaultRun$Resp/100),
                 row.names = c("O2 Flux [mmol/m2/d]",
                               "Total Resp [mmol/m2/d]"))
    } else
      data.frame("Current"=c(out$O2flux/100,out$Resp/100 ),
                 row.names = c("O2 Flux [mmol/m2/d]",
                               "Total Resp [mmol/m2/d]"))
    },include.rownames=T)
}

ui <- fluidPage(
  headerPanel("O2 model"),
  sidebarLayout(
    sidebarPanel(
      h3("Parameters"),
      sliderInput("porosity", label = "porosity",
                   min = 0.3, max = 0.99,  value = 0.8, step = 0.05, width=100),
      sliderInput("minrate", label = "minrate - mmol/m3/d",
                   min = 0.0, max = 80000,  value = parms["minrate"], step = 1000, width=100),
      sliderInput("ks", label = "ks - mmol/m3",
                   min = 0.01, max = 20, value = parms["ks"], step = 0.1, width=100),
      sliderInput("O2BW", label = "O2BW - mmol/m3",
                  min = 0, max = 600, value = parms["O2BW"], step = 10, width=100),
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