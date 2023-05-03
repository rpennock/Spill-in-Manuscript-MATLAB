# Spill-in-Manuscript-MATLAB
MATLAB codes used for analysis in Nat Comm shared AMPAR manuscript

The codes deposited here were used to analyze line scan images in our manuscript titled 'Afferent convergence onto a shared population
of interneuron AMPA receptors' in Nature Communications. Below is a description of the work flow and a brief description of each script.

1) AlignImages.m
    -This script is used to align individual line scan sweeps spatially so that accurate analysis of the amplitude and spread of 
     fluorescence signals can be performed. Make sure that the width (spatial dimension) and length (temporal dimension) of all
     images is equal or there will be an error. Images from a channel detecting a fill dye (e.g. Alexa 594) and a fluorescent indicator
     (e.g. Fluo5F) are needed for each sweep.
    -Make sure images are named sequentially so that cropped images produced by this script are in the correct order (e.g. 010123_Cell1_Fluoo5F_001, *_002, etc.).
     Indicate sweeps where failures occurred (e.g. 010123_Cell1_Fluo5F_F003). This will allow the next script in the work flow to separate
     sweeps with a signal from sweeps with failures.
    -The script will detect all images for each channel then detect peaks in the fluorescene of the dye channel. The user will then be 
     prompted to select a peak to crop the image around. There is an option to manually change the values for the peaks if needed if there 
     is poor detection on a small number of sweeps.
    -Cropped images will be placed in a new folder using the file names indicated.
    
2) CalculateTransients.m
    -The cropped images produced in (1) will now be used to calculate deltaG/R and deltaG/G. The sweeps with a signal present must be indicated
     before running the script. If sweeps with failures were named as suggested above this will be simple. If not each sweep with a signal will need
     to be indicated individually. The onset of any provided stimulus and the sampling rate of the line scans being analyzed will also need to be
     indicated.
    -ftWidth variable defines the number of pixels in the spatial dimension used to calculate dG/R and dG/G.
    -Filter parameters and scale bars can be adjusted as needed.
    -Upon running script an image of dG/R or dG/G will appear to allow the visualization of the signals before selecting the center of the signals.
     Dashed lines are placed at regular intervals to help with accurate centering of the signal when calculating dG/R or dG/G. Multiple runs
     may be required in some cases to get this correct.
    -Sweeps containing a signal will be averaged and plotted as a trace of dG/R or dG/G as a function of time. A similar trace will be overlaid
     for sweeps where failures occurred.
    -Raw data for these traces will be saved as .txt files in a new folder titled 'dGoverR' and 'dGoverG'
     
3) CaTransientsFx.m
    -After calculating dG/R running the next script fits each line scan (i.e. each row of the image) with a Gaussian function to measure the spread
     of the fluorescence signal over time.
    -The '*_dGoverR.mat' file from (2) will need to be copied into the script. All the other parameters were saved during (2) and will be loaded as
     long as the same source folder is being used.
    -The pixel size will need to be acquired from the image metadata and entered.
    -Filter parameters can be adjusted as needed
    -A sample of fluorescence profiles post stimulus will be plotted before fitting with Gaussians
    -There is an option to discard line scans post stimulus but before a signal is present, leaving only scans that can be fit to a Gaussian
    -All scans post stimulus for the time interval indicated (timePostStim) will be fit to Gaussians (**long time intervals can take a while
     to run since every row is plotted)
    -A file of the Gaussian fit parameters as a function of time will be saved as a .txt file
    
4) CaTPeakToPeak.m
    -When Ca2+ transients arising from multiple point sources are present in an image, this script allows all peaks to be fit with Gaussian functions
     and measures the distance between the peaks.
     
