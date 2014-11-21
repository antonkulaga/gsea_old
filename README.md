GSEA
====
This is a temporal project for Gene Sets Enrichment Analyses. 
The data in examples belongs to Alexey Moskalev, see paper Mining Gene Expression Data for Pollutants (Dioxin, Toluene, Formaldehyde) and Low Dose of Gamma-Irradiation
A part of code belong to developers of R (bioconductor) package Piano, see http://www.bioconductor.org/packages/release/bioc/html/piano.html 

Getting started
---------------

You can run this project locally. The easiest way to do this is to run RStudio.
For this you need to have all packages mentioned in server.R in addition to Shiny package. 
In RStudio you 
 - load Shiny library library('Shiny')
 - set working directory to the one with ui.R and server.R with command setwd("full path here")
 - initiate Shiny app with runApp() command.
 

TROUBLESHOOTING:
----------------
 
 If you are using Ubuntu/Mint you may also need to install some linux packages (look into errors and see what is missing, sometimes you have to google for them), like:
 ```
 sudo apt-get install  libxml2-dev
 sudo apt-get install  gfortran
 sudo apt-get install r-cran-rgl 
 ```
You can also have an issue with permissions "installed directory not writable, cannot update packages". One of the possible
ways (but not the safe one) to fix it is to get to know where R is installed by:

```
whereis R 
```
and then apply to them

```
sudo chown -R 777 /pathto/R/folder
```
