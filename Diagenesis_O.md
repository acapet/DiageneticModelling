Oxygen Diagenesis
================
Karline Soetaert and Arthur Capet
July 2017

Set-up
======

We start with a simple oxygen diagenetic model. We need

### 1. to load librairies R-package ReacTran.

``` r
require(ReacTran, quietly = TRUE) # transport models in aquatic systems 
require(marelac, quietly = TRUE)  # toolbox for aquatic sciences
```

    ## 
    ## Attaching package: 'marelac'

    ## The following objects are masked from 'package:oce':
    ## 
    ##     coriolis, gravity

### 2. to set a spatial framework

Here, a vertical 1D grid with a 100 levels over 10cm.

``` r
grid <- setup.grid.1D(N = 100, L = 10)
```

### 3. Model equations

We consider **diffusion** and constant **respiration** above a certain depth.
$$\\frac{\\partial O\_2(z)}{\\partial t}=\\frac{\\partial O\_2(z)}{\\partial t}|\_{transport}-\\gamma(z)$$

*O*<sub>2</sub> is a concentration in *m**m**o**l* *m*<sub>*l**i**q**u**i**d*</sub><sup>−3</sup> and we have to consider porosity *ϕ*.
$$ \\frac{\\partial O\_2(z)}{\\partial t} =  \\frac{1}{\\phi}\\frac{\\partial}{\\partial z}\[\\phi D \\frac{\\partial C}{\\partial z }\]  - \\gamma (z) $$
 Boundary conditons are imposed as:

-   constant concentration at the sediment-water interface
    *O*<sub>2</sub>|<sub>*z* = 0</sub> = *O*<sub>2, *b**o**t**t**o**m*</sub>
-   nul gradient at the lower cell
    $$ \\frac{\\partial O\_2}{\\partial z}|\_{z=deep}=0$$

``` r
O2model <- function (t, O2, p) {
  with (as.list(p), {
    # Transport term (using ReacTran routines)
    O2tran <- tran.1D(C = O2, C.up = O2BW, D = D, VF = porosity, dx = grid)
    # Respiration
    O2cons <- minrate*(grid$x.mid < mindepth)    
    # the function returns the time derivative
    list(O2tran$dC - O2cons, O2cons = O2cons)
  })
}
```

### 4. to assign a value to the model parameters

``` r
parms <- c(
  porosity = 0.8    , # -
  minrate  = 1    , # nmol O2/cm3/d - oxygen consumption rate
  mindepth = 5    , # cm            - depth below which minrate = 0
  O2BW     = 300  , # nmol/cm3       - bottom water oxygen concentration
  D        = 1      # cm2/d         - molecular diffusion coefficient
)
```

Simulation
==========

Let us now compute the steady state solutions for three set of parameters (considering different respiration rates) and plot the ouputs.

``` r
parms2 <- parms3 <- parms
parms2["minrate"] <- 2
parms3["minrate"] <- 4

# random initial conditions
IC<-runif(length(grid$x.mid))

# computes the steay-state solution
out  <- steady.1D(y = IC, parms = parms, func = O2model, nspec = 1, names = "O2")
out2 <- steady.1D(y = IC, parms = parms2, func = O2model, nspec = 1, names = "O2")
out3 <- steady.1D(y = IC, parms = parms3, func = O2model, nspec = 1, names = "O2")

plot(out, out2, out3, xyswap = TRUE, xlab = "mmol/m3", ylab = "cm", grid = grid$x.mid)
```

<img src="Diagenesis_O_files/figure-markdown_github/runs-1.png" style="display: block; margin: auto;" />

------------------------------------------------------------------------

Specific Outputs
================

So far the variable ´out´ contains the steady-state solution and the diagnostic ´O2cons´

``` r
str(out)
```

    ## List of 2
    ##  $ y     : num [1:100, 1] 300 299 299 298 298 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : NULL
    ##   .. ..$ : chr "O2"
    ##  $ O2cons: num [1:100] 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "precis")= num [1:3] 6.51e+02 5.03e-03 2.13e-09
    ##  - attr(*, "steady")= logi TRUE
    ##  - attr(*, "class")= chr [1:3] "steady1D" "rootSolve" "list"
    ##  - attr(*, "dimens")= num 100
    ##  - attr(*, "nspec")= num 1
    ##  - attr(*, "ynames")= chr "O2"

The model declaration is modified to add a diagnostic : the oxygen flux at the sediment-water interface ´O2flux´.

``` r
O2model <- function (t, O2, p) {
  with (as.list(p), {
    O2tran <- tran.1D(C = O2, C.up = O2BW, D = D, VF = porosity, dx = grid)
    O2cons <- minrate*(grid$x.mid < mindepth)    
    list(O2tran$dC - O2cons, O2cons = O2cons,O2flux=O2tran$flux.up)
  })
}
```

``` r
out  <- steady.1D(y = IC, parms = parms, func = O2model, nspec = 1, names = "O2")
out2 <- steady.1D(y = IC, parms = parms2, func = O2model, nspec = 1, names = "O2")
out3 <- steady.1D(y = IC, parms = parms3, func = O2model, nspec = 1, names = "O2")

plot(x = c(parms["minrate"] , parms2["minrate"] , parms3["minrate"]), y = c(out$O2flux , out2$O2flux, out3$O2flux),
     xlab="Mineralization rate - [nmol O2/cm3/d]" , ylab="Oxygen Flux @ SWI - [nmol/cm2/d]")
```

<img src="Diagenesis_O_files/figure-markdown_github/unnamed-chunk-4-1.png" style="display: block; margin: auto;" /> At steady-state, the flux at the interface equals the total amount of oxygen consumed in the sediment. Here, consumption takes place only in the liquid, so we have to consider porosity in the integral.
$$ \\int\\limits\_\\infty^0 \\mathbf{\\phi} \\gamma  ~ dz =\\phi  \\gamma\_0 L $$
 -------------------

SensRange
=========

The library FME is a toolbox for exploring model characteristics.

``` r
library(FME)
```

    ## Loading required package: coda

We will only use the function sensRange, that allows to visualize the variability of model outputs when a parameter is allowed to vary within a given range.

´sensRange´ needs

-   a data frame to precise which parameter should vary and within which range.

``` r
parRange <- data.frame(min=1,max=10)
rownames(parRange)<-"minrate"
```

-   A function that reecievs a set of parameter and returns the response variable to be investigated (here we define two functions, for ´O2´ and ´O2flux´)

``` r
solveO2<-function(pars){
  out  <- steady.1D(y = IC, parms = pars, func = O2model, nspec = 1, names = "O2")
  return(data.frame(Depth=grid$x.mid , O2=out$y))
}

solveO2f<-function(pars){
  out  <- steady.1D(y = IC, parms = pars, func = O2model, nspec = 1, names = "O2")
  return(data.frame(Depth=0,O2flux=out$O2flux))
}


sR<-sensRange(func = solveO2, sensvar = "O2", parms=parms, 
              parRange = parRange, num=50, dist='grid')

sRf<-sensRange(func = solveO2f, sensvar = "O2flux", parms=parms,
               parRange = parRange, num=50, dist='grid')

plot(summary(sR),xyswap = TRUE, xlab = "nmol/cm3", ylab = "cm")
plot(summary(sRf), xlab = "nmol/cm2/d",  main="O2 flux")
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-7-1.png)![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-7-2.png)

Exercices
=========

1. Explore the range of model response to variations of other parameters
------------------------------------------------------------------------

*Hint* : Just adapt parRange in the upper code. [Start Here](Diagenesis_O_Ex1.Rmd)

2. Read measured data (Flux and/or Profiles), include them to the plot
----------------------------------------------------------------------

*Hint* : use read.table to read the data. Add them to the previous plot using ´line´ and ´points´

Example:

``` r
measuredO2flux = 20

O2prof<-matrix(ncol=2,  data=c(0.5,1,2,3,4,285,280,270,260,255))
colnames(O2prof)<-c("depth","O2")

plot(summary(sRf), xlab = "nmol/cm2/d", main="O2 flux")
points(x= 25,y=1, col="red", cex=5 )

plot(summary(sR),xyswap = TRUE, xlab = "nmol/cm3", ylab = "cm")
lines(y=O2prof[,'depth'] ,x=O2prof[,'O2'], col="red" )
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-8-1.png)![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-8-2.png) Can you try parameters value that match your measurements ?

``` r
plot(out, xyswap = TRUE, xlab = "mmol/m3", ylab = "cm", grid = grid$x.mid)
lines(y=c(0.5,1,2,3,4) ,x=c(295,280,270,260,255), col="red" )

barplot( height = c(out$O2flux, measuredO2flux ), ylab="Oxygen Flux @ SWI - [nmol/cm2/d]", names.arg = c("Model","Measured"),col=c('black','red'))
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-9-1.png)![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-9-2.png)

Same with the other parameters, also testing different distributions

\`\`\`r library(FME)

parRange &lt;- data.frame(min= c( .5 , 1 , 1 , 200 , .1 ), max= c( 1 , 10 , 10 , 400 , 10 )) rownames(parRange)&lt;-names(parms) \`\`\`

|          |    min|  max|
|----------|------:|----:|
| porosity |    0.5|    1|
| minrate  |    1.0|   10|
| mindepth |    1.0|   10|
| O2BW     |  200.0|  400|
| D        |    0.1|   10|

``` r
par(mfrow=c(2,2))

parRangelocal<-parRange[1,]

sR<-sensRange(func = solveO2, sensvar = "O2", parms=parms, parRange = parRangelocal,
              num=100, dist='grid')
sRf <- sensRange(func = solveO2f, sensvar = "O2flux", parms=parms, parRange = parRangelocal,
              num=100, dist='grid')
plot(summary(sR),xyswap = TRUE, xlab = "mmol/m3", ylab = "cm")
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
plot(summary(sRf),xyswap = TRUE, xlab = "nmol/m2/d", ylab = "cm")
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-12-2.png)

``` r
hist(sR[,1], xlab = rownames(parRangelocal))
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-12-3.png)

------------------------------------------------------------------------

``` r
parRangelocal<-parRange[2,]

sR<-sensRange(func = solveO2, sensvar = "O2", parms=parms, parRange = parRangelocal,
              num=100, dist='grid')
sRf <- sensRange(func = solveO2f, sensvar = "O2flux", parms=parms, parRange = parRangelocal,
              num=100, dist='grid')
plot(summary(sR),xyswap = TRUE, xlab = "mmol/m3", ylab = "cm")
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
plot(summary(sRf),xyswap = TRUE, xlab = "nmol/m2/d", ylab = "cm")
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-13-2.png)

``` r
hist(sR[,1], xlab = rownames(parRangelocal))
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-13-3.png)

------------------------------------------------------------------------

``` r
parRangelocal<-parRange[3,]
sR<-sensRange(func = solveO2, sensvar = "O2", parms=parms, parRange = parRangelocal,
              num=100, dist='grid')
sRf <- sensRange(func = solveO2f, sensvar = "O2flux", parms=parms, parRange = parRangelocal,
              num=100, dist='grid')
plot(summary(sR),xyswap = TRUE, xlab = "mmol/m3", ylab = "cm")
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
plot(summary(sRf),xyswap = TRUE, xlab = "nmol/m2/d", ylab = "cm")
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-14-2.png)

``` r
hist(sR[,1], xlab = rownames(parRangelocal))
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-14-3.png)

------------------------------------------------------------------------

``` r
parRangelocal<-parRange[4,]

sR<-sensRange(func = solveO2, sensvar = "O2", parms=parms, parRange = parRangelocal,
              num=100, dist='grid')
sRf <- sensRange(func = solveO2f, sensvar = "O2flux", parms=parms, parRange = parRangelocal,
              num=100, dist='grid')
plot(summary(sR),xyswap = TRUE, xlab = "mmol/m3", ylab = "cm")
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
plot(summary(sRf),xyswap = TRUE, xlab = "nmol/m2/d", ylab = "cm")
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-15-2.png)

``` r
hist(sR[,1], xlab = rownames(parRangelocal))
```

![](Diagenesis_O_files/figure-markdown_github/unnamed-chunk-15-3.png)

------------------------------------------------------------------------
