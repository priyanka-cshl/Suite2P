GUI_Suite2P_Main.m

1) Edit "toolbox_path = 'C:\Users\florin\Desktop\software\Suite2P_CSHL';" to match location of your local code folder

2) Edit file paths based on your preferences in "Line 63:71". Easiest way is find your host name, and have host specific settings. 
   You can also set all kinds of defaults here.
    ... to find your hostname type "!hostname" in command line
    ------------------
    Filepaths - notes
    ------------------
    handles.filepaths.Data(1): Root storage path where all your data is located - can be a local server
        The assumed folder heirarchy is RootStorage\MouseName\Date\ExperimentFolder
        The key restrictions are:
            ExperimentFolder Name should be a number eg. '1', '2' ..    
            You can select multiple Experiment folders within one session - these will be treated as data from the same FOV - common registration and ROI extraction
            You can also specify multiple sessions - these will be analysed independently, but  useful for batch processing.

    handles.filepaths.Data(2) = Local folder - files will read/written from this location frequently - better if this drive is an SSD

    handles.filepaths.Data(3) = Save path for analyzed data

    handles.filepaths.Data(4) = Save path for Registered Tiffs. You can uncheck saving of registered tiffs in the GUI if you don't want to save the registered Tiffs,
                                and save only the registration offsets

3) Compile the deconvolution file by entering the following in cmd line. You might need to download a compiler - google it.
   "mex -largeArrayDims SpikeDetection/deconvL0.c"
    or "mex -largeArrayDims SpikeDetection/deconvL0.cpp"
    If you don't manage to compile this - uncheck deconvolve spikes, change cluster mode to 'standard' and neuropil subtraction to 'none' in the GUI.

4) GPU ... 