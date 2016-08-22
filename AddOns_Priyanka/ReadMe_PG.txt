To Run GUI_Suite2P_Main.m (Registration, ROI extraction, Spike deconvolution)

1) Edit file paths in My_Filepaths.m (Suite2P\AddOns_Priyanka\GUI_Suite2P_Main\My_FilePaths.m) based on your preferences.
If you will run the same code on many machines, you can have host specific settings in the same file. To find your hostname type "!hostname" in command line
   You can also set all kinds of defaults here.
    
    ------------------
    Filepaths - notes
    ------------------
    handles.filepaths.Data(1): Root storage path where all your data is located - can be a local server
        The assumed folder heirarchy is RootStorage\MouseName\Date\ExperimentFolder
        You can select multiple Experiment folders within one session - these will be treated as data from the same FOV - common registration and ROI extraction
        You can also specify multiple sessions - these will be analysed independently, but  useful for batch processing.

    handles.filepaths.Data(2) = Local folder - files will read/written from this location frequently - better if this drive is an SSD

    handles.filepaths.Data(3) = Save path for analyzed data

    handles.filepaths.Data(4) = Save path for Registered Tiffs. You can uncheck saving of registered tiffs in the GUI if you don't want to save the registered Tiffs,
                                and save only the registration offsets

    Edit "toolbox_path = 'C:\Users\florin\Desktop\software\Suite2P_CSHL';" to match location of your local code folder

3) Compile the deconvolution file by entering the following in cmd line. You might need to download a compiler - Go to Home/AddOns. 
    In AddOn explorer, search mingw - install the compiler. You will need to login to your mathworks account.
    Once the compiler is installed, execute the command below from within the Suite2P directory
   "mex -largeArrayDims SpikeDetection/deconvL0.c"
    or "mex -largeArrayDims SpikeDetection/deconvL0.cpp"
    If you don't manage to compile this - uncheck deconvolve spikes, change cluster mode to 'standard' and neuropil subtraction to 'none' in the GUI.

4) GPU ... uncheck use GPU if you don't have a GPU. 
    Updated documentation coming soon.



To Run GUI_Suite2P_ROIs.m (Post extraction pruning of ROIs).
    Updated documentation coming soon.
