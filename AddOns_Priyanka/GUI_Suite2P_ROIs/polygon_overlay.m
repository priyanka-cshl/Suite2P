
function polygon_overlay(h)

if h.dat.procmap
    I = h.dat.mimg_proc(:,:,h.dat.map);
else
    I = h.dat.mimg(:,:,h.dat.map);
end

% create a mask
imgMaskColor = [0 0.8 0.8];
tmp = ones(size(I));
imgMask = cat(3, imgMaskColor(1)*tmp, imgMaskColor(2)*tmp, imgMaskColor(3)*tmp);

Mask_multiplier_1 = str2num(h.mask_multiplier_1.String);

axes(h.axes2);
hold on
hg6 = image(imgMask);
my_polygon = poly2mask(h.polygon(:,1), h.polygon(:,2), size(I,1), size(I,2));
hg6.AlphaData = (Mask_multiplier_1)*my_polygon;
hold off
axis off

drawnow






