
function redraw_meanimg_overlay(h, whichcase)
if nargin < 2
    whichcase = 0;
end

if h.dat.procmap
    I = h.dat.mimg_proc(:,:,h.dat.map);
else    
    I = h.dat.mimg(:,:,h.dat.map);
end

mu = median(I(:));
sd1 = mean(abs(I(I<mu) - mu));
sd2 = mean(abs(I(I>mu) - mu));

axes(h.axes2); imagesc(I, mu + 5*[-sd1 sd2]);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
colormap('gray')
axes(h.axes3); imagesc(I, mu + 5*[-sd1 sd2]);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
colormap('gray')

% ROI image
allROIs = hsv2rgb(cat(3, h.dat.img1.H, h.dat.img1.Sat, h.dat.img1.V));
allROIs = min(allROIs, 1);
allROIs = rgb2gray(allROIs);

allROIs_excluded = hsv2rgb(cat(3, h.dat.img2.H, h.dat.img2.Sat, h.dat.img2.V));
allROIs_excluded = min(allROIs_excluded, 1);
allROIs_excluded = rgb2gray(allROIs_excluded);

% create a mask
imgMaskColor = [0 0.8 0.8];
tmp = ones(size(I));
imgMask = cat(3, imgMaskColor(1)*tmp, imgMaskColor(2)*tmp, imgMaskColor(3)*tmp);

Mask_multiplier_1 = str2num(h.mask_multiplier_1.String);
Mask_multiplier_2 = str2num(h.mask_multiplier_2.String);

% display masked ROIs
axes(h.axes2);
hold on
hg2 = image(imgMask);
hg2.AlphaData = (Mask_multiplier_1)*allROIs;
hold off
axis off

if whichcase == 1
    axes(h.axes3);
    hold on
    hg3 = image(imgMask);
    hg3.AlphaData = (Mask_multiplier_1)*allROIs_excluded;
    hold off
    axis off
end

if whichcase == 2 % compare 2 sets of ROIs
    % create a new mask
    newROIs = getfield(h.ROIs_compare,char(fieldnames(h.ROIs_compare)));
    % get offsets to align images
    rows_to_cut = str2num(h.row_offset.String);
    columns_to_cut = str2num(h.column_offset.String);
    my_newROIs = 0*newROIs; 
    my_newROIs(1:end-rows_to_cut,1:end-columns_to_cut) = newROIs(rows_to_cut+1:end,columns_to_cut+1:end);
    newROIs = my_newROIs;
    imgMaskColor = [0.8 0 0];
    tmp = ones(size(newROIs));
    imgMask = cat(3, imgMaskColor(1)*tmp, imgMaskColor(2)*tmp, imgMaskColor(3)*tmp);
    axes(h.axes2);
    hold on
    hg4 = image(imgMask);   
    hg4.AlphaData = (Mask_multiplier_2)*newROIs;
    hold off
    axis off
end

Sat = h.dat.img1.Sat;
imgMaskColor = [0 0.8 0];
tmp = ones(size(Sat));
imgMask = cat(3, imgMaskColor(1)*tmp, imgMaskColor(2)*tmp, imgMaskColor(3)*tmp);
axes(h.axes2);
hold on
hg5 = image(imgMask);
hg5.AlphaData = (Mask_multiplier_1)*~Sat;
hold off
axis off
        
drawnow






