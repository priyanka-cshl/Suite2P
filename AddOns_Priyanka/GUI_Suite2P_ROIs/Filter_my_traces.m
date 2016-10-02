function [filtered] = Filter_my_traces(raw);
filtered = 0*raw;
for i = 1:size(raw,1)
    temp = raw(i,:);
    filtered(i,:) = my_conv_local(medfilt1(temp, 3), 3);
end
end