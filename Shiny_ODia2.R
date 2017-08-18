library("deSolve")

require(ReacTran, quietly = TRUE) # Reaction-Transport models in aquatic systems 
require(marelac, quietly = TRUE)  # toolbox for aquatic sciences

grid <- setup.grid.1D(N = 100, L = 10, dx.1 = 0.01)

H2Smodel <- function (t, S, p) {
  with (as.list(p), {
    O2<-S[1:100]
    H2S<-S[101:200]
    
    # Transport term (using ReacTran routines)
    O2tran <- tran.1D(C = O2, C.up = O2BW, D = DO2, VF = porosity, dx = grid)
    H2Stran <- tran.1D(C = H2S, C.up = H2SBW, D = DH2S, VF = porosity, dx = grid)
    
    # Respiration
    resp <- minrate*(grid$x.mid < mindepth)
    O2LIM<-O2/(O2+kO2lim)
    
    O2resp  <- resp*O2LIM
    H2Sprod <- resp-O2resp
    
    # H2S oxidation
    H2SOx<-rH2Sox*H2S*O2/(O2+ksO2H2Sox)
    
    dO2  <-  O2tran$dC - O2resp  - H2SOx
    dH2S <-  H2Stran$dC + H2Sprod - H2SOx
    
    return( list(cbind(dO2,dH2S),
                 H2SOx=H2SOx,
                 O2resp = O2resp,
                 O2flux = O2tran$flux.up,
                 Resp=sum(resp*porosity*grid$dx)) ) 
  })
}

parms <- c(
  porosity  = 0.8    , # -
  minrate   = 30     , # nmol O2/cm3/d - oxygen consumption rate
  mindepth  = 5      , # cm            - depth below which minrate = 0
  O2BW      = 300    , # nmol/cm3       - bottom water oxygen concentration
  DO2       = as.numeric(diffcoeff(species="O2")*86400*1e4)/(1-log(0.8*0.8)),     # cm2/d - molecular diffusion coefficient
  DH2S      = as.numeric(diffcoeff(species="H2S")*86400*1e4)/(1-log(0.8*0.8)),     # cm2/d - molecular diffusion coefficient
  H2SBW     = 0      , # nmol/cm3       - bottom water H2S concentration
  kO2lim    = .3     , # Oxygen limitation for oxic respiration
  rH2Sox    = 5   , # rate of H2S oxidation
  ksO2H2Sox = 10       # Oxygen limitation for H2S oxidation 
)

ICO2 <- rep.int(c(300),length(grid$x.mid))
ICH2S <- rep.int(c(1),length(grid$x.mid))

IC<-cbind(ICO2,ICH2S)

# computes the steady-state solution
DefaultRun  <- ode.1D(times=seq(0,100,.1),y = IC, parms = parms, func = H2Smodel, nspec = 2, names = c("O2","H2S"))
image(DefaultRun,legend = T,ylim=c(10,0),grid = grid$x.mid)

DefaultRun  <- steady.1D(y = IC, parms = parms, func = H2Smodel, nspec = 2, names = c("O2","H2S"))
plot(DefaultRun,xyswap=T)


server <- function(input, output,session) {
  observeEvent(input$resetButton, {
    updateNumericInput(session, "porosity", value = parms[["porosity"]])
    updateNumericInput(session, "minrate", value = parms[["minrate"]])
    updateNumericInput(session, "mindepth", value = parms[["mindepth"]])
    updateNumericInput(session, "O2BW", value = parms[["O2BW"]])
  })
  
  output$model <- renderPlot({
    IC <- cbind(rep.int(300,length(grid$x.mid)),rep.int(1,length(grid$x.mid)))
    Parms <- c(porosity=input$porosity,
               minrate=input$minrate,
               mindepth=input$mindepth,
               O2BW=input$O2BW,
               DO2 = parms[["DO2"]],
               DH2S      = parms[["DH2S"]]      ,     # cm2/d - molecular diffusion coefficient
               H2SBW     = parms[["H2SBW"]]      , # nmol/cm3       - bottom water H2S concentration
               kO2lim    = parms[["kO2lim"]]    , # Oxygen limitation for oxic respiration
               rH2Sox    = parms[["rH2Sox"]]    , # rate of H2S oxidation
               ksO2H2Sox = parms[["ksO2H2Sox"]]       # Oxygen limitation for H2S oxidation 
    )
    # computes the steady-state solution
    out  <- steady.1D(y = IC, parms = Parms, func = H2Smodel, nspec = 2, names = c("O2","H2S"),pos=T)
    
    if (input$default) {
      plot(out, DefaultRun , which=c('O2','H2S','O2resp','H2SOx'),
           xyswap = TRUE, xlab = c("mmol/m3","mmol/m3","mmol/m3/d","mmol/m3/d"),
           ylab = "cm", grid = grid$x.mid)
      legend("bottomright", col = 1:2, legend = c("current", "default"), lty = 1:2)
    } else  
      plot(out, which=c('O2','H2S','O2resp','H2SOx'),
           xyswap = TRUE, xlab =  c("mmol/m3","mmol/m3","mmol/m3/d","mmol/m3/d"),
           ylab = "cm", grid = grid$x.mid)
  })
  
  output$table <- renderTable({
    IC <- cbind(rep.int(200,length(grid$x.mid)),rep.int(2,length(grid$x.mid)))
    Parms <- c(porosity=input$porosity,
               minrate=input$minrate,
               mindepth=input$mindepth,
               O2BW=input$O2BW,
               DO2 = parms[["DO2"]],
               DH2S      = parms[["DH2S"]]      ,     # cm2/d - molecular diffusion coefficient
               H2SBW     = parms[["H2SBW"]]      , # nmol/cm3       - bottom water H2S concentration
               kO2lim    = parms[["kO2lim"]]    , # Oxygen limitation for oxic respiration
               rH2Sox    = parms[["rH2Sox"]]    , # rate of H2S oxidation
               ksO2H2Sox = parms[["ksO2H2Sox"]]       # Oxygen limitation for H2S oxidation 
    )
    # computes the steady-state solution
    out  <- steady.1D(y = IC, parms = Parms, func = H2Smodel, nspec = 2, names = c("O2","H2S"))
    
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
  headerPanel("O2 + H2S model"),
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