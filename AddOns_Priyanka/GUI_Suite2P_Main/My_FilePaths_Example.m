% function to define local paths for GUI_Suite2P_Main

function [handles] = My_FilePaths(handles)
if strcmp(handles.computername,'marbprec')
    handles.filepaths.Data(1) = {'Z:\photoncerber_data'}; % Root storage
    handles.filepaths.Data(2) = {'C:\Users\florin\Desktop\DATA_processed\'}; % local read/write folder
    handles.filepaths.Data(3) = {'C:\Users\florin\Desktop\DATA_processed\F'}; % Save path
    handles.filepaths.Data(4) = {'C:\Users\florin\Desktop\DATA_processed\R'}; % Registered Tiffs
    handles.cluster_settings.Data(1:2) = [500; 200];
    handles.cluster_settings.Data(4) = 250;
    handles.session_list.String = {''};
    
    % toolbox
    handles.ops0.toolbox_path = 'C:\Users\florin\Desktop\software\Suite2P_CSHL';
else
    handles.ops0.toolbox_path = [];
end
end