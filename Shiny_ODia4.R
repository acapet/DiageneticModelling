library("deSolve")

require(ReacTran, quietly = TRUE) # Reaction-Transport models in aquatic systems 
require(marelac, quietly = TRUE)  # toolbox for aquatic sciences

grid <- setup.grid.1D(N = 100, L = 10, dx.1 = 0.01)

TOCmodel <- function (t, S, p) {
  with (as.list(p), {
    O2<-S[1:100]
    ODU<-S[101:200]
    TOC<-S[201:300]
    
    # Transport term (using ReacTran routines)
    O2tran <- tran.1D(C = O2, C.up = O2BW, D = DO2, VF = porosity, dx = grid)
    ODUtran <- tran.1D(C = ODU, C.up = ODUBW, D = DODU, VF = porosity, dx = grid)
    TOCtran <- tran.1D(C = TOC, flux.up = OMflux , D = Dbioturb, VF = 1- porosity, dx = grid)
    
    # Respiration
    resp  <- TOC*minrate  # ! per volume of SOLIDS !
    respL <- resp*(1-porosity)/porosity # ! per volume of LIQUID !
    
    O2LIM <- O2/(O2+kO2lim)
    
    O2resp  <- respL*O2LIM
    ODUprod <- respL-O2resp
    
    # ODU oxidation
    OduOx<-rODUox*ODU*O2/(O2+ksO2oduox)
    
    #Irrigation 
    alpha<-alpha0*(grid$x.mid < irrigdepth)
    O2irrig  <- - alpha * (O2-O2BW)
    ODUirrig <- - alpha * (ODU-ODUBW)
    
    dO2  <-  O2tran$dC - O2resp  - OduOx + O2irrig
    dODU <-  ODUtran$dC + ODUprod - OduOx +ODUirrig
    dTOC <-  TOCtran$dC - resp
  
    return( list(cbind(dO2,dODU,dTOC),
                 OduOx=OduOx,
                 O2resp = O2resp,
                 O2flux = O2tran$flux.up,
                 ODUflux = ODUtran$flux.up,
                 O2irrigflux = sum(porosity*O2irrig*grid$dx),
                 ODUirrigflux = sum(porosity*ODUirrig*grid$dx),
                 Resp=sum(respL*porosity*grid$dx) ,
                 OrgC=TOC*12*1e-9/2.5*100 )) # Convert to percent of dry weight.
            
  })
}

parms <- c(
  porosity  = 0.8    , # -
  minrate   = 1/100   , # /d - mineralisation rate
  OMflux    = 1000    , # nmol/cm2/d Flux of Organic Matter
  O2BW      = 300    , # nmol/cm3       - bottom water oxygen concentration
  DO2       = as.numeric(diffcoeff(species="O2")*86400*1e4)/(1-log(0.8*0.8)),     # cm2/d - molecular diffusion coefficient
  DODU      = as.numeric(diffcoeff(species="H2S")*86400*1e4)/(1-log(0.8*0.8)),     # cm2/d - molecular diffusion coefficient
  ODUBW     = 0      , # nmol/cm3       - bottom water ODU concentration
  kO2lim    = .3     , # Oxygen limitation for oxic respiration
  rODUox    = 5      , # rate of ODU oxidation
  ksO2oduox = 10     , # Oxygen limitation for ODU oxidation 
  alpha0    = 0      , # Irrigation coefficient
  irrigdepth  = 3    , # Irrigation Depth
  Dbioturb   = 1/365      # cm2/d Bioturbation coeficient
)

ICO2 <- rep.int(c(300),length(grid$x.mid))
ICODU <- rep.int(c(1),length(grid$x.mid))
ICTOC <- rep.int(c(.03),length(grid$x.mid))

IC<-cbind(ICO2,ICODU,ICTOC)

# computes the steady-state solution
#DefaultRun  <- ode.1D(times=seq(0,100,.1),y = IC, parms = parms, func = TOCmodel, nspec = 2, names = c("O2","ODU"))
#image(DefaultRun,legend = T,ylim=c(10,0),grid = grid$x.mid)

# computes the steady-state solution
DefaultRun  <- steady.1D(y = IC, parms = parms, func = TOCmodel, nspec = 3, names = c("O2","ODU","TOC"),pos=T)
plot(DefaultRun,xyswap=T,grid=grid$x.mid)


server <- function(input, output,session) {
  observeEvent(input$resetButton, {
    updateNumericInput(session, "porosity", value = parms[["porosity"]])
    updateNumericInput(session, "minrate", value = parms[["minrate"]])
    updateNumericInput(session, "O2BW", value = parms[["O2BW"]])
    updateNumericInput(session, "alpha0", value = parms[["alpha0"]])
    updateNumericInput(session, "irrigdepth", value = parms[["irrigdepth"]])
    updateNumericInput(session, "OMflux", value = parms[["OMflux"]])
    updateNumericInput(session, "Dbioturb", value = parms[["Dbioturb"]])
  })
  
  output$model <- renderPlot({
    IC <- cbind(rep.int(300,length(grid$x.mid)),rep.int(1,length(grid$x.mid)),rep.int(1,length(grid$x.mid)))
    Parms <- c(porosity=input$porosity,
               minrate=input$minrate,
               O2BW=input$O2BW,
               DO2        = parms[["DO2"]],
               DODU       = parms[["DODU"]]      , # cm2/d - molecular diffusion coefficient
               ODUBW      = parms[["ODUBW"]]     , # nmol/cm3       - bottom water ODU concentration
               kO2lim     = parms[["kO2lim"]]    , # Oxygen limitation for oxic respiration
               rODUox     = parms[["rODUox"]]    , # rate of ODU oxidation
               ksO2oduox  = parms[["ksO2oduox"]] , # Oxygen limitation for ODU oxidation 
               alpha0     = input$alpha0,
               irrigdepth = input$irrigdepth,
               OMflux     = input$OMflux,
               Dbioturb   = input$Dbioturb
    )
    # computes the steady-state solution
    out  <- steady.1D(y = IC, parms = Parms, func = TOCmodel, nspec = 3, names = c("O2","ODU","TOC"),pos=T)
    
    if (input$default) {
      plot(out, DefaultRun , which=c('O2','ODU','OrgC'),
           xyswap = TRUE, xlab = c("mmol/m3(liquids)","mmol/m3(liquids)","% d.w."),
           ylab = "cm", grid = grid$x.mid)
      legend("bottomright", col = 1:2, legend = c("current", "default"), lty = 1:2)
    } else  
      plot(out, which=c('O2','ODU','OrgC'),
           xyswap = TRUE, xlab =  c("mmol/m3(liquids)","mmol/m3(liquids)","% d.w."),
           ylab = "cm", grid = grid$x.mid)
  })
  
  output$table <- renderTable({
    IC <- cbind(rep.int(300,length(grid$x.mid)),rep.int(1,length(grid$x.mid)),rep.int(1,length(grid$x.mid)))
    Parms <- c(porosity=input$porosity,
               minrate=input$minrate,
               O2BW=input$O2BW,
               DO2        = parms[["DO2"]],
               DODU       = parms[["DODU"]]      , # cm2/d - molecular diffusion coefficient
               ODUBW      = parms[["ODUBW"]]     , # nmol/cm3       - bottom water ODU concentration
               kO2lim     = parms[["kO2lim"]]    , # Oxygen limitation for oxic respiration
               rODUox     = parms[["rODUox"]]    , # rate of ODU oxidation
               ksO2oduox  = parms[["ksO2oduox"]] , # Oxygen limitation for ODU oxidation 
               alpha0     = input$alpha0,
               irrigdepth = input$irrigdepth,
               OMflux     = input$OMflux,
               Dbioturb   = input$Dbioturb
    )
    # computes the steady-state solution
    out  <- steady.1D(y = IC, parms = Parms, func = TOCmodel, nspec = 3, names = c("O2","ODU","TOC"),pos=T)
    
    if (input$default) {
      data.frame("Current"=c(out$O2flux,
                             out$ODUflux,
                             out$O2irrigflux,
                             out$ODUirrigflux,
                             sum(c(out$O2flux,-out$ODUflux,out$O2irrigflux,-out$ODUirrigflux)),
                             out$Resp),
                 "Default"=c(DefaultRun$O2flux,
                             DefaultRun$ODUflux,
                             DefaultRun$O2irrigflux,
                             DefaultRun$ODUirrigflux,
                             sum(c(DefaultRun$O2flux,-DefaultRun$ODUflux,DefaultRun$O2irrigflux,-DefaultRun$ODUirrigflux)),
                             DefaultRun$Resp),
                 row.names = c("O2 diffusive Flux [nmol/cm2/d]",
                               "ODU diffusive Flux [nmol/cm2/d]",
                               "O2 irrigative Flux [nmol/cm2/d]",
                               "ODU irrigative Flux [nmol/cm2/d]",
                               "sum",
                               "Total Respiration"))
    } else
    data.frame("Current"=c(out$O2flux,out$ODUflux,out$O2irrigflux,out$ODUirrigflux),
               row.names = c("O2 diffusive Flux [nmol/cm2/d]",
                             "ODU diffusive Flux [nmol/cm2/d]",
                             "O2 irrigative Flux [nmol/cm2/d]",
                             "ODU irrigative Flux [nmol/cm2/d]"))
    },include.rownames=T)
}

ui <- fluidPage(
  headerPanel("O2+ODU+irrigation+OM model"),
  sidebarLayout(
    sidebarPanel(
      h3("Parameters"),
      sliderInput("porosity", label = "porosity - []",
                   min = 0.3, max = 0.99,  value = 0.8, step = 0.05, width=100),
      sliderInput("minrate", label = "minrate - [/d]",
                   min = 0.0, max = 1/10,  value = parms["minrate"], step = 0.001, width=100),
      sliderInput("O2BW", label = "O2BW - [mmol/m3]",
                  min = 0, max = 600, value = parms["O2BW"], step = 10, width=100),
      sliderInput("alpha0", label = "alpha0 - [d]",
                  min = 0, max = 3, value = parms["alpha0"], step = .01, width=100),
      sliderInput("irrigdepth", label = "irrigdepth - [cm]",
                  min = 0, max = 7, value = parms["irrigdepth"], step = .5, width=100),
      sliderInput("OMflux", label = "OMflux - [nmol/cm2/d]",
                  min = 0, max = 5000, value = parms["OMflux"], step = 100, width=100),
      sliderInput("Dbioturb", label = "Dbioturb - [cm2/d]",
                  min = 0, max = .03, value = parms["Dbioturb"], step = .0003, width=100),
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