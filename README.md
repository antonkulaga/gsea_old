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