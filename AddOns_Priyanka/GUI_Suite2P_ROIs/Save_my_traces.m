function [] = Save_my_traces(h)

% get filename from user
[file1, path] = uiputfile([h.dat.filename(1:end-4) '_proc.mat'],'Save Session As');
[~,file2] = fileparts(file1);
file2 = [file2,'_condensed.mat'];

disp('saving session ...')
% basic params
Output.Lx = h.dat.cl.Lx; % image size
Output.Ly = h.dat.cl.Ly; % image size
Output.total_rois = sum(h.dat.cl.iscell); % number of rois chosen
Output.total_frames = size(h.dat.F.trace,2);

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
Output.traces.roi.raw = h.dat.Fcell{1}(Output.rois_indices,:);
Output.traces.roi.filtered = 0*Output.traces.roi.raw;
Output.traces.neuropil.raw = h.dat.FcellNeu{1}(Output.rois_indices,:);
Output.traces.neuropil.filtered = 0*Output.traces.neuropil.raw;

for i = 1:Output.total_rois
    temp = Output.traces.roi.raw(i,:);
    Output.traces.roi.filtered(i,:) = my_conv_local(medfilt1(temp, 3), 3);
    temp = Output.traces.neuropil.raw(i,:);
    Output.traces.neuropil.filtered(i,:) = my_conv_local(medfilt1(temp, 3), 3);
end

% neuropil subtraction coefficients and events
count = 0;
for i = 1:size(Output.rois_indices,1)
    if h.dat.cl.isroi(Output.rois_indices(i,1))
        count = size(find(h.dat.cl.isroi(1:Output.rois_indices(i,1))),2);
        Output.traces.neuropil.coeffs(i,:) = h.dat.cl.dcell{count}.B(2:3);
        
        % events
        my_spike_times = zeros(1,Output.total_frames);
        my_spike_times(1,h.dat.cl.dcell{count}.st) = h.dat.cl.dcell{count}.c;
        Output.traces.events(i,:) = my_spike_times;
    else
        Output.traces.neuropil.coeffs(i,:) = [NaN NaN];
        Output.traces.events(i,1:Output.total_frames) = NaN;
    end
end


% save the condensed traces
save(fullfile(path,file2), 'Output');

% save the whole session as well
h.dat.F.trace = [];
dat = h.dat;
save(fullfile(path, file1), 'dat');
disp(['saved files: ', file1, ' & ', file2]);
end