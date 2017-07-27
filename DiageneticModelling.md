Oxygen in porous non-permeable sediments
================
Arthur Capet, Marilaure Grégoire, Karline Soetaert
July 2017

Set-up
======

We start with a simple diffusive model for oxygen, considering a constant respiration rate above a certain depth.

For this we need

### 1. Librairies

``` r
require(ReacTran, quietly = TRUE) # Reaction-Transport models in aquatic systems 
require(marelac, quietly = TRUE)  # toolbox for aquatic sciences
```

    ## 
    ## Attaching package: 'marelac'

    ## The following objects are masked from 'package:oce':
    ## 
    ##     coriolis, gravity

### 2. Grid

We define a vertical 1D grid with 100 levels over 10cm. It is good to have finer spatial resolution in the upper cells, where steeper gradients occur.

``` r
grid <- setup.grid.1D(N = 100, L = 10, dx.1 = 0.01)
plot(grid$x.mid)
```

![](DiageneticModelling_files/figure-markdown_github/grid-1.png)

### 3. Model

The unique **state variable** is *O*<sub>2</sub>(*z*, *t*), the oxygen concentration in pore waters, given in *m**m**o**l* *m*<sub>*l**i**q**u**i**d*</sub><sup>−3</sup>. We consider **diffusion** and constant **respiration** above a certain depth.
$$\\frac{\\partial O\_2(z)}{\\partial t}=\\frac{\\partial O\_2(z)}{\\partial t}|\_{transport}+\\frac{\\partial O\_2(z)}{\\partial t}|\_{consumption}$$

$$ \\frac{\\partial O\_2(z)}{\\partial t} = \\frac{\\partial}{\\partial z}\[ D \\frac{\\partial C}{\\partial z }\]  - \\gamma (z) $$
 For boundary conditions, we impose :

-   constant concentration at the sediment-water interface
    *O*<sub>2</sub>|<sub>*z* = 0</sub> = *O*<sub>2, *B**W*</sub>
-   nul gradient at the lower cell
    $$ \\frac{\\partial O\_2}{\\partial z}|\_{z=\\text{deepest cell}}=0$$

This model is expressed as a **model function** that receives **current time** (t), **current state** (O2) and **parameters** (p) and returns **time derivatives**.

``` r
O2model <- function (t, O2, p) {
  with (as.list(p), {
    # Transport term (using ReacTran routines)
    O2tran <- tran.1D(C = O2, C.up = O2BW, D = D, VF = porosity, dx = grid)
    # Respiration
    O2cons <- minrate*(grid$x.mid < mindepth)    
    
    dO2<-O2tran$dC - O2cons
    # The model function returns the time derivatives as first argument.
    # The next arguments in the list are diagnostics that can be used for later visualisation.
    return( list(dO2, O2cons = O2cons) ) 
  })
}
```

### 4. Parameters

``` r
parms <- c(
  porosity = 0.8    , # -
  minrate  = 10     , # nmol O2/cm3/d - oxygen consumption rate
  mindepth = 5      , # cm            - depth below which minrate = 0
  O2BW     = 300    , # nmol/cm3      - bottom water oxygen concentration
  D        = as.numeric(diffcoeff(species="O2")*86400*1e4)/(1-log(0.8*0.8)) # cm2/d - molecular diffusion coefficient, corrected for tortuosity
)
```

Simulation
==========

To solve the model and plot the ouput we first need to define initial conditions. Then, we'll use the package `deSolve` to solve the model differential equations. For a dynamic solution (resolving the system evolution with time) we use the function `ode`. This function requires the time steps `times` at which we'll compute the new model states. We also specify that we model one species (*nspec*), and give a name to this species (*names*).

``` r
# Define initial conditions
IC <- rep(200, times = length(grid$x.mid))

# computes the dynamical solution for `parms`
out  <- ode.1D(times=seq(0,50,.1), y = IC, parms = parms, func = O2model, nspec = 1, names = "O2")

# plot the output
image(out, ylim=c(10,0), grid=grid$x.mid, legend=TRUE, ylab="Sediment Depth - [cm]")
```

<img src="DiageneticModelling_files/figure-markdown_github/runs-1.png" style="display: block; margin: auto;" />

The system converges towards an equilibrium solution.

Instead of computing the full dynamical evolution of the system, from the initial conditions until it reaches the steady state, it is possible to compute directly the steady-state solutions. The problem is to find *C* such that $\\frac{\\partial C}{\\partial t}=0$. For this, the package `rootSolve` provides a function answering our needs : `steady.1D`.

``` r
# computes the steady-state solution
out  <- steady.1D(y = IC, parms = parms, func = O2model, nspec = 1, names = "O2")
plot(out, xyswap = TRUE, xlab = "mmol/m3", ylab = "cm", grid = grid$x.mid)
```

<img src="DiageneticModelling_files/figure-markdown_github/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

------------------------------------------------------------------------

Exercise
========

**Q** Add the oxygen flux at the sediment-water interface as another diagnostic value to the list returned by the `O2model` (in addition to the time derivatives). The value can be found from the list returned by the function `tran.1D`: `O2tran$flux.up`.

**Q**: Explore how the the oxygen profile and the oxygen flux at the sediment water interface at steady state varies when a parameter is modified.

``` r
#Redefine the model function 'O2model' with an additional diagnostic (O2flux=..) in the list returned at the end of the function.

#Recompute the model solution, and check the value of the O2 flux at the interface. 

#Create a new vector with new parameter values: `parms2`

#Compute the new steady-state solution : out2

#Plot both solution on a graph   ( use `plot(out,out2,...) `)

#Compare the fluxes from out and out2
```

! WHEN YOU'RE DONE ! Execute the following code chunks in Rstudio, to explore interactively the role of each parameter.

( PS: you might have to execute : `install.packages('shiny', repos='https://cran.rstudio.com/')` ) ( PS2: Find a published version of the running Shiny App [here](https://acapet.shinyapps.io/DiageneticModelling/) )

``` r
library(shiny)
runApp('Shiny_ODia.R')
```

------------------------------------------------------------------------

O2 and H2S
==========

If you change mineralisation rates to high values, you'll notice that oxygen can get negative. This is because the oxygen consumption is independent of the oxygen concentration.

We should limit O2 consumption when oxygen is not available by using a Monod function, and define a new state variable `H2S`, that is produced whenever oxygen is not sufficient for organic matter mineralization.

Below you'll find the previous O2 model, and the updated H2S model.

The main changes are highlighted here:

-   Previously, the *State Vector* contained only one variable, `O2`, and had 100 elements: one value for each grid point. Now it contains two variables `O2` and `H2S`, concatenated into a longer vector (the size of the new *State Vector* will thus be 2\*100. We used `c` to concatenate the vectors).
-   H2S is a solute variable that is transported similarly to `O2`. So we had to define new parameters:
    -   `DH2S`, the diffusion coefficient, from package marelac.
    -   `H2SBW`, the H2S bottom water concentration.
-   The total respiration term `resp` is the same as we used before, but now only a part is respired with oxygen, while the remaining part produced reduced substances (think of it as if organic matter was oxidised by sulfate, thereby prH2Scing hydrogen sulfide).
-   Oxic respiration is limited by oxygen availability. The oxygen consumed for respiration is equal to the total respiration times a limiting function *O*<sub>2</sub>/(*O*<sub>2</sub> + *k*<sub>*O*<sub>2</sub>*l**i**m*</sub>).
-   The remaining part of `resp` is a source of `H2S`.
-   In addition, the `H2S` is oxidised (consuming oxygen) when it diffuses in the oxic regions of the sediments. This term can be expressed as *r*<sub>*H*2*S**o**x*</sub>.*H*2*S*.*O*<sub>2</sub>/(*O*<sub>2</sub> + *k*<sub>*O*<sub>2</sub>*H*2*S**o**x*</sub>) (it only occurs when there is oxygen) and is a loss for `H2S` and for `O2`.
-   We use `c` to concatenate dO2 and dH2S into one long vector to be returned by the model function.
-   To ensure model stability, we added the option `pos=T` as last argument to the solving function `steady.1d`. It forbids negative value, that could trigger instability.

``` r
#  O2 model 

#Grid
grid <- setup.grid.1D(N = 100, L = 10, dx.1 = 0.01)

#Model
O2model <- function (t, O2, p) {
  with (as.list(p), {
    # Transport term (using ReacTran routines)
    O2tran <- tran.1D(C = O2, C.up = O2BW, D = D, VF = porosity, dx = grid)
    # Respiration
    resp <- minrate*(grid$x.mid < mindepth)    
    # the function returns the time derivative
    dO2<-O2tran$dC - resp
    return( list(dO2, O2cons = resp) ) 
  })
}

# Parameters
parms <- c(
  porosity = 0.8    , # -
  minrate  = 10    , # nmol O2/cm3/d - oxygen consumption rate
  mindepth = 5    , # cm            - depth below which minrate = 0
  O2BW     = 300  , # nmol/cm3       - bottom water oxygen concentration
  D        = as.numeric(diffcoeff(species="O2")*86400*1e4)/(1-log(0.8*0.8))     # cm2/d - molecular diffusion coefficient
)

ICO2<-rep(200,length(grid$x.mid))
IC<-ICO2

out  <- steady.1D(y = IC, parms = parms, func = O2model, nspec = 1, names = "O2")
plot(out, xyswap = TRUE, xlab = "mmol/m3", ylab = "cm", grid = grid$x.mid)
```

![](DiageneticModelling_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
# This is the O2 + H2S model  

#Grid
grid <- setup.grid.1D(N = 100, L = 10, dx.1 = 0.01)

#Model
H2Smodel <- function (t, S, p) {
  with (as.list(p), {
    O2 <-S[1:100]
    H2S<-S[101:200]
    
    # Transport term (using ReacTran routines)
    O2tran <- tran.1D(C = O2, C.up = O2BW, D = DO2, VF = porosity, dx = grid)
    H2Stran <- tran.1D(C = H2S, C.up = H2SBW, D = DH2S, VF = porosity, dx = grid)
    
    # Respiration
    resp <- minrate*(grid$x.mid < mindepth)
    O2LIM<-O2/(O2+kO2lim)
    
    O2resp  <- resp*O2LIM
    H2Sprod <- resp-O2resp
    
    H2SOx<-rH2Sox*H2S*O2/(O2+ksO2H2Sox)
    
    dO2  <-  O2tran$dC - O2resp  - H2SOx
    dH2S <-  H2Stran$dC + H2Sprod - H2SOx
    
    return( list(c(dO2,dH2S), H2SOx=H2SOx, O2resp = O2resp, O2flux = O2tran$flux.up) ) 
  })
}

parms <- c(
  porosity  = 0.8    , # -
  minrate   = 40     , # nmol O2/cm3/d - oxygen consumption rate
  mindepth  = 5      , # cm            - depth below which minrate = 0
  O2BW      = 300    , # nmol/cm3       - bottom water oxygen concentration
  DO2       = as.numeric(diffcoeff(species="O2")*86400*1e4)/(1-log(0.8*0.8)),     # cm2/d - molecular diffusion coefficient
  DH2S      = as.numeric(diffcoeff(species="H2S")*86400*1e4)/(1-log(0.8*0.8)),     # cm2/d - molecular diffusion coefficient
  H2SBW     = 0      , # nmol/cm3       - bottom water H2S concentration
  kO2lim    = 0.3    , # Oxygen limitation for oxic respiration
  rH2Sox    = 5   , # rate of H2S oxidation
  ksO2H2Sox = 10       # Oxygen limitation for H2S oxidation 
)

ICO2 <- rep.int(c(300),length(grid$x.mid))
ICH2S <- rep.int(c(1),length(grid$x.mid))

IC <- cbind(ICO2,ICH2S)

# computes the dynamic solution
#DefaultRun  <- ode.1D(times=seq(0,100,.1),y = IC, parms = parms, func = H2Smodel, nspec = 2, names = c("O2","H2S"),pos=T)
#image(DefaultRun,legend = T,ylim=c(10,0),grid = grid$x.mid)

# computes the steady-state solution
DefaultRun  <- steady.1D(y = IC, parms = parms, func = H2Smodel, nspec = 2, names = c("O2","H2S"), pos=TRUE)
plot(DefaultRun,xyswap=T, grid=grid$x.mid)
```

![](DiageneticModelling_files/figure-markdown_github/unnamed-chunk-6-1.png)

For an interactive version, execute the following code chunks in Rstudio. ( PS: Find a published version of the running Shiny App [here](https://acapet.shinyapps.io/diageneticmodelling__oxygen_and_hydrogen_sulphide/) )

``` r
library(shiny)
runApp('Shiny_ODia2.R')
```

------------------------------------------------------------------------

Exercise
========

Bioirrigation can be represented in 1D models by allowing non-local exchange, ie. a term *α*(*C*(*z*)−*C*(0)) by which solutes are exchanged between any depth *z* and bottom waters, without having to diffuse through near-surface sediments.

**Q**: Add a biorrigation term to the equations, for both O2 and H2S.

-   Start by copying the O2 + H2S model.
-   Within the model function, create a vector of irrigation rates: *α*. It should have a fixed value *α*<sub>0</sub> above a given depth (e.g. a parameter `irrigdepth`), and 0 below (remember the expression we used to define the respiration rate in the first model).
    $$
    \\alpha(z) = \\begin{cases}
    \\alpha\_0 \\text{ for $z &lt; irrigDepth$} \\\\
    0 \\text{ for $z&gt; irrigDepth$} 
    \\end{cases}
    $$
-   Add a term *O*2*i**r**r**i**g*(*z*)=*α*(*z*)\[*O*<sub>2</sub>(*z*)−*O*<sub>2</sub>(0)\] to the time derivative of *O*<sub>2</sub>. Do the same for H2S.
-   Provide the following diagnostics :

    -   Diffusive flux, provided from the transport terms.
    -   Irrigative flux, computed by integrating the irrigation terms $\\int\\limits\_\\infty^0\\phi(z)\\alpha(z)\[ O\_2(z) - O\_2(0)\]dz$ (*d**z* can be obtained from `grid$dx`)
-   How much oxygen is consumed by oxic mineralization, and how much for H2S oxidation ? Is that ratio dependent of the irrigation rate ?

For an interactive version, execute the following code chunks in Rstudio. <https://acapet.shinyapps.io/diageneticmodelling__oxygen_and_hydrogen_sulphide/>

``` r
library(shiny)
runApp('Shiny_ODia3.R')
```

------------------------------------------------------------------------

O2 + H2S + Organic Matter
=========================

Let us now consider more realistic respiration rates by introducing a solid state variable for organic matter.

The equations of the system looks like this :

-   Organic Matter
    $$ \\frac{\\partial TOC(z)}{\\partial t} = \\underbrace{\\frac{\\partial}{\\partial z}\[ D\_b \\frac{\\partial TOC}{\\partial z }\]}\_{bioturbation}  - \\underbrace{\\gamma TOC}\_{resp}$$

-   Oxygen

$$ \\frac{\\partial O\_2(z)}{\\partial t} = \\underbrace{\\frac{\\partial}{\\partial z}\[ D\_{O\_2} \\frac{\\partial O\_2(z)}{\\partial z }\]}\_{diffusion} 
- \\underbrace{\\gamma TOC.\\frac{1-\\phi}{\\phi}.\\frac{O\_2}{O\_2+k\_{O\_2lim}}}\_{oxic ~respiration} \\\\
- \\underbrace{\\alpha (O\_2(z)-O\_2(0))}\_{irrigation}
- \\underbrace{r\_{H2Sox}.H2S.\\frac{O\_2}{O\_2+k\_{O\_2H2SOx}}}\_{H2S~oxidation}$$

-   H2S

$$ \\frac{\\partial H2S(z)}{\\partial t} = \\underbrace{\\frac{\\partial}{\\partial z}\[ D\_{H2S} \\frac{\\partial H2S(z)}{\\partial z }\]}\_{diffusion} 
- \\underbrace{\\gamma TOC.\\frac{1-\\phi}{\\phi}.(1-\\frac{O\_2}{O\_2+k\_{O\_2lim}})}\_{anoxic ~respiration} \\\\
- \\underbrace{\\alpha (H2S(z)-H2S(0))}\_{irrigation} 
- \\underbrace{r\_{H2Sox}.H2S.\\frac{O\_2}{O\_2+k\_{O\_2H2SOx}}}\_{H2S~oxidation}$$

-   For the boundary conditions of the solid fractions we impose `OMflux`: the sedimenting organic matter input.

``` r
grid <- setup.grid.1D(N = 100, L = 10, dx.1 = 0.01)

TOCmodel <- function (t, S, p) {
  with (as.list(p), {
    O2<-S[1:100]
    H2S<-S[101:200]
    TOC<-S[201:300]
    
    # Transport term (using ReacTran routines)
    O2tran <- tran.1D(C = O2, C.up = O2BW, D = DO2, VF = porosity, dx = grid)
    H2Stran <- tran.1D(C = H2S, C.up = H2SBW, D = DH2S, VF = porosity, dx = grid)
    TOCtran <- tran.1D(C = TOC, flux.up = OMflux , D = Dbioturb, VF = 1- porosity, dx = grid)
    
    # Respiration
    resp  <- TOC*minrate  # ! per volume of SOLIDS !
    respL <- resp*(1-porosity)/porosity # ! per volume of LIQUID !
    
    O2LIM <- O2/(O2+kO2lim)
    
    O2resp  <- respL*O2LIM
    H2Sprod <- respL-O2resp
    
    # H2S oxidation
    H2SOx<-rH2Sox*H2S*O2/(O2+ksO2H2Sox)
    
    #Irrigation 
    alpha    <- alpha0*(grid$x.mid < irrigdepth)
    O2irrig  <- - alpha * (O2-O2BW)
    H2Sirrig <- - alpha * (H2S-H2SBW)
    
    dO2  <-  O2tran$dC - O2resp  - H2SOx + O2irrig
    dH2S <-  H2Stran$dC + H2Sprod - H2SOx +H2Sirrig
    dTOC <-  TOCtran$dC - resp
  
    return( list(cbind(dO2,dH2S,dTOC),
                 H2SOx=H2SOx,
                 O2resp = O2resp,
                 O2flux = O2tran$flux.up,
                 H2Sflux = H2Stran$flux.up,
                 O2irrigflux = sum(porosity*O2irrig*grid$dx),
                 H2Sirrigflux = sum(porosity*H2Sirrig*grid$dx),
                 OrgC=TOC*12*1e-9/2.5*100) # Convert to percent of dry weight.
            ) 
  })
}

parms <- c(
  porosity  = 0.8    , # -
  minrate   = 1/100   , # /d - mineralisation rate
  OMflux    = 1000    , # nmol/cm2/d Flux of Organic Matter
  O2BW      = 300    , # nmol/cm3       - bottom water oxygen concentration
  DO2       = as.numeric(diffcoeff(species="O2")*86400*1e4)/(1-log(0.8*0.8)),     # cm2/d - molecular diffusion coefficient
  DH2S      = as.numeric(diffcoeff(species="H2S")*86400*1e4)/(1-log(0.8*0.8)),     # cm2/d - molecular diffusion coefficient
  H2SBW     = 0      , # nmol/cm3       - bottom water H2S concentration
  kO2lim    = .3     , # Oxygen limitation for oxic respiration
  rH2Sox    = 5      , # rate of H2S oxidation
  ksO2H2Sox = 10     , # Oxygen limitation for H2S oxidation 
  alpha0    = 0      , # Irrigation coefficient
  irrigdepth  = 3    , # Irrigation Depth
  Dbioturb   = 1/365      # cm2/d Bioturbation coeficient
)

ICO2 <- rep.int(c(300),length(grid$x.mid))
ICH2S <- rep.int(c(1),length(grid$x.mid))
ICTOC <- rep.int(c(.03),length(grid$x.mid))

IC<-cbind(ICO2,ICH2S,ICTOC)

# computes the steady-state solution
#DefaultRun  <- ode.1D(times=seq(0,100,.1),y = IC, parms = parms, func = TOCmodel, nspec = 2, names = c("O2","H2S"))
#image(DefaultRun,legend = T,ylim=c(10,0),grid = grid$x.mid)

# computes the steady-state solution
DefaultRun  <- steady.1D(y = IC, parms = parms, func = TOCmodel, nspec = 3, names = c("O2","H2S","TOC"),pos=T)
plot(DefaultRun,xyswap=T,grid=grid$x.mid)
```

![](DiageneticModelling_files/figure-markdown_github/unnamed-chunk-10-1.png)

[Interactive version:](https://acapet.shinyapps.io/diageneticmodelling_oxygen_sulphide_bioirrigation_orgc/)

``` r
library(shiny)
runApp('Shiny_ODia4.R')
```

------------------------------------------------------------------------

Exercise
========

Q1 : Add a solute variable DIC to the previous model.

Q2 : Add ammonium and nitrate.

-   Ammonium is produced from OM mineralisation (e.g. consider a C/N ratio of 10).
-   Nitrification can be expressed as `Nitri = rnit*NH3*O2/(O2+ksO2nitri)`
-   As before, compute first the total respiration on the basis of the TOC variable.
-   Total resp is then partitionned by computing first the relative importance of the differetn pathways :

    -   Oxicminlim = O2/(O2+ksO2oxic)
    -   Denitrilim = (1-O2/(O2+kinO2denit)) \* NO3/(NO3+ksNO3denit), ie. denitrification is limited by nitrate and inhibited by oxygen.
    -   Anoxiclim = (1-O2/(O2+kinO2anox))\*(1-NO3/(NO3+kinNO3anox)), anoxic mineralisation is inhibited by both nitrate and oxygen.
    -   Rescale = 1/(Oxicminlim+Denitrilim+Anoxiclim), we normalize, so that the three mineralisation pathways sums to 1.
-   The oxydants consumption is then provided by

    OxicMin = totresp*Oxicminlim*Rescale ! Oxic mineralisation Denitrific = totresp*Denitrilim*Rescale ! Denitrification AnoxicMin = totresp*Anoxiclim *Rescale ! Anoxic mineralisation, ie. H2S prH2Sction