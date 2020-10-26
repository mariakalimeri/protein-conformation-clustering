# Protein conformational clustering

Source code for a shiny app to perform clustering on protein conformations as 
sampled by Molecular Dynamics (MD) simulations. 

If you are interested in the science of this app, want to use it but don't know 
how, contact me and I will provide two video-demos. If you know what is the 
purpose of clustering on protein MD trajectories is, then using the app should 
be straight forward, assuming you know how to fire up the shiny app GUI 
(see below). While using the app, there are explanatory tips within each step 
along the process. 

You do not need to worry about memory overheads in case your trajectory files 
are huge, frames are read one at a time. The code currently analyzes `.dcd` 
files. `.xtc` file option may become available in the future but for now see 
below for converting `xtc` to `dcd`.

This code was written in the context of a collaboration with [Dr. Fabio Sterpone, 
LBT/IBPC, CNRS, Paris (France)](https://sites.google.com/site/sterponefabio/) 


## Credits
Includes code for reading binary `dcd` files from the 
[bio3d](http://thegrantlab.org/bio3d/) package implememnted by B. Grant.

## Dependencies 

Make sure you have the following:

- R version 3.3.1 or newer 
- RStudio version 0.99.903 or newer

As well as the following R packages:

- shiny v0.13.2
- shinyjs v0.7  
- shinyBS v0.61   
- igraph v1.0.1 
- bio3d v2.2-4
- ggplot2 v2.1.0 

(If newer R or package versions give you trouble, contact me)

## Usage 

From within the R console and the source code's root directory, type

```r
shiny::runApp(appDir = ".")
```

(or if outside the root directory, specify the path to the shiny app source code
as the `appDir` argument.)


## Known Issues

1. The browse button, in the first tab of the app, works optimally only when you 
run the app "in Window" or in "Viewer Pane" within RStudio. When run outside 
RStudio, i.e. in a browser, the pop-up window doesn't quite pop-up. 
In this case, you will have to type the path to the files explicitly in 
the file text entries (i.e. when browsing for a `dcd` and a `pdb` file). 
The reason for the issue: the code needs only the path for the binary trajectory 
file in order to read it and perform the clustering step by step (`dcd` files are 
typically big). However, the current "file upload" options available in shiny 
directly load the whole files into memory, which is not an option here. As a 
workaround we use the function "file.choose", in r-base, to browse to the 
location of the file. The "file.choose" function works with a browse button 
only when used with Rstudio. Outside Rstudio paths should be given explicitly.  

2. (Only for linux users): When on linux, you actually have to use an external 
browser rather than the RStudio window or viewer pane, as the download buttons 
don't work in the viewer on linux. 

3. The app currently works only with `.dcd` binary files, i.e. NAMD trajectory 
files. Gromacs users may use [`mdconvert`](http://mdtraj.org/1.9.0/mdconvert.html) 
to convert trajectories from `xtc` to `dcd`.

4. (Optional protein comformation rendering): The `render.sh` script will not 
work by default on any machine. This script is called when the "render" button 
in "Step 2" is clicked. The script calls [`VMD`](https://www.ks.uiuc.edu/Research/vmd/) 
as well as the tachyon render, that usually comes by default with `VMD`. So if 
you want to use the rendering in Step 2, you need to 

1. have `VMD` installed and 
2. open the `render.sh` script and change the `VMD` and `tachyon` paths to point 
to the correct locations in your system. 

If you don't want or can't render protein conformation in the app rendering, 
you can use `VMD` manually by using the cluster leaders `pdb` file that is also 
provided in step 2.

Contact me for any issue.

