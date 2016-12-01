# ProteinleaderClustering
A shiny app to perform clustering on protein Molecular Dynamics (MD) trajectories 

If you are interested in the science of this app, want to use it but don't understand how, contact me and I will provide a video-demo. If you know what clustering on protein MD trajectories is, then using the app should be straight forward. You will find additional explanatory tips within each step of the app. You do not need to worry about memory overheads in case your trajectory files are huge (they will be read frame by frame).

It currently analyzes dcd trajectory files. Xtc option will be available in the future (see below for converting xtc to dcd)

Acknowledgements:
Includes code from the bio3d package implememnted by B. Grant for reading the binary dcd files.

=================== TECHNICALITIES =================================================

DEPENDENCIES:
Make sure you have the following installed

- R version 3.3.1 or newer 
- Rstudio version 0.99.903 or newer

R packages
- shiny v0.13.2
- shinyjs v0.7  
- shinyBS v0.61   
- igraph v1.0.1 
- bio3d v2.2-4
- ggplot2 v2.1.0 

(You can check which version of R you got by typing sessionInfo() from within R. That way you also get the respective info for all the packages you have loaded in the current R session. If newer R or package versions give you trouble, contact me)

HOW TO RUN IT:
Make sure you have all the provided files in the same folder, including the subfolder www/ with its current content. Open the ui.R function with Rstudio and press the Run App button on the top right (or type runApp() in the R terminal below). Make sure you use the app within Rstudio to have all the features working properly!

NOTES
Known issue 1: The browse button, in the first tab of the app, works only within Rstudio. Outside Rstudio you will have to type the path to the files explicitly in the file text windows (i.e. when browsing for dcd and pdb files). The reason for the issue is that the code needs only the path for the trajectory file in order to read it and perform the clustering step by step (dcd files are typically huge). However, the current "file upload" options available in shiny directly load the whole files into memory, which is not an option here. As a workaround we use the function "file.choose", in r-base, to browse to the location of the file. The "file.choose" function works with a browse button only when used with Rstudio. Outside Rstudio paths should be given explicitly.  

Known issue 2 (only for linux users): When on linux, use the app in an external browser rather than the Rstudio viewer (i.e. runApp() and subsequently click on button "Open in Browser" on the top left of the Rstudio viewer) as the download buttons don't work in the viewer on linux. 

Note 3: The app currently works only with .dcd (NAMD) trajectory files. Gromacs users may use mdconvert to convert trajectories from xtc to dcd.

Note 4 (optional automatic rendering): The render.sh script will not work by default on any machine. This script is run when the "render" button in Step 2 of the app is pressed. The script calls VMD as well as the tachyon render, that usually comes by default with VMD. So if you want to use the automatic rendering of Step 2, you need to 1) have VMD in your machine and 2) open the render.sh script and change the VMD and tachyon paths to point to the correct locations in your system. If you don't want automatic rendering, it is ok, use VMD manually by using the leaders PDB file that is also provided in step 2.

Contact me for any issue.

