function [] = Save_my_traces(h)

% get filename from user
timestamp = ['_',datestr(now,'yyyy-mm-dd'),'_',datestr(now,'HH:MM')];
timestamp = strrep(timestamp,':','-');
file1 = [h.dat.filename(1:end-4),timestamp,'.mat'];
file2 = [h.dat.filename(1:end-4),timestamp,'_condensed.mat'];

disp('saving session ...')
% basic params
Output.Lx = h.dat.cl.Lx; % image size
Output.Ly = h.dat.cl.Ly; % image size
Output.total_rois = sum(h.dat.cl.iscell); % number of rois chosen
Output.directories = h.dat.ops.expts;

% to get the indices of chosen rois
Output.rois_indices = find(h.dat.cl.iscell)';

% to get the image with roi pixels indexed by actual cluster index
Output.ROImasks.all = h.dat.res.iclust;
chosen_roi_pixels = ismember(h.dat.res.iclust, Output.rois_indices);
Output.ROImasks.chosen = h.dat.res.iclust .* chosen_roi_pixels;

% to get the variance explained image
h.dat.img0.V = reshape(h.dat.res.M, h.dat.cl.Ly, h.dat.cl.Lx) / h.dat.cl.MeanM;
Output.ROImasks.varianceexp = h.dat.img0.V .* chosen_roi_pixels;

% mean image
Output.meanimage.green = h.dat.mimg(:,:,2);
Output.flattenedimage.green = h.dat.mimg_proc(:,:,2);
if size(h.dat.mimg,3)>2
    Output.meanimage.red = h.dat.mimg(:,:,3);
    Output.flattenedimage.red = h.dat.mimg_proc(:,:,3);
end

% fluorescence traces
start_frame = 1;
end_frame =  0;

for j = 1:size(h.dat.Fcell,2) % each recording session
    Output.total_frames(j) = size(h.dat.Fcell{j},2);
    end_frame = end_frame + Output.total_frames(j);
    Output.traces(j).roi.raw = h.dat.Fcell{j}(Output.rois_indices,:);
    Output.traces(j).neuropil.raw = h.dat.FcellNeu{j}(Output.rois_indices,:);

    % neuropil subtraction coefficients and events
    count = 0;
    for i = 1:size(Output.rois_indices,1)
        if h.dat.cl.isroi(Output.rois_indices(i,1))
            count = size(find(h.dat.cl.isroi(1:Output.rois_indices(i,1))),2);
            Output.traces(j).neuropil.coeffs(i,:) = h.dat.cl.dcell{count}.B(2:3);
            
            % events
            my_spike_times = zeros(1,sum(Output.total_frames));
            my_spike_times(1,h.dat.cl.dcell{count}.st) = h.dat.cl.dcell{count}.c;
            Output.traces(j).events(i,1:Output.total_frames(j)) = my_spike_times(1,start_frame:end_frame);
        else
            Output.traces.neuropil.coeffs(i,:) = [NaN NaN];
            Output.traces(j).events(i,1:Output.total_frames(j)) = NaN;
        end
    end
    start_frame = end_frame + 1;
end
% save the condensed traces
save(file2, 'Output');

% save session parameters as well
h.dat = rmfield(h.dat,'Fcell');
h.dat = rmfield(h.dat,'FcellNeu');
h.dat = rmfield(h.dat,'F');
dat = h.dat;
save(file1, 'dat');
disp(['saved files: '])
disp(file1);
disp(file2);
end