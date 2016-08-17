GUI_Suite2P_Main.m

1) Edit "Line 151: toolbox_path = 'C:\Users\florin\Desktop\software\Suite2P_CSHL';" to match location of your local code folder
2) Edit file paths based on your preferences in "Line 63:71". Easiest way is find your host name, and have host specific settings. 
   You can also set all kinds of defaults here.
    ... to find your hostname type "!hostname" in command line
3) compile the deconvolution file by entering the following in cmd line. You might need to download a compiler - google it.
   "mex -largeArrayDims SpikeDetection/deconvL0.c"
    or "mex -largeArrayDims SpikeDetection/deconvL0.cpp"