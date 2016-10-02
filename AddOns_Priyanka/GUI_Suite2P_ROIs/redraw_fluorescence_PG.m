function redraw_fluorescence_PG(h)
axes(h.axes4)
hold off
MyCase = h.fluorescence_display_settings.SelectedObject.String;
% plot raw F trace
switch MyCase
    case {'Raw', 'Raw, neuropil'} % plot raw trace
        temp = double(h.dat.F.trace(h.dat.F.ichosen,:));
        if h.show_unfiltered.Value
            plot(temp,'color',0.5*[1 1 1]);
            hold on
        end
        plot(my_conv_local(medfilt1(temp, 3), 3),'y');
        axis tight
        hold on
        if isequal (MyCase, 'Raw, neuropil') % plot also neuropil
            if isfield(h.dat.F, 'neurop')
                temp = double(h.dat.F.neurop(h.dat.F.ichosen,:));
                if h.show_unfiltered.Value
                    plot(temp,'color',0.5*[1 1 1]);
                end
                plot(my_conv_local(medfilt1(temp, 3), 3),'c');
            end
        end
    case {'neuropil subtracted','neuropil subtracted, ceiled'}
        if isfield(h.dat.F, 'neurop')
            temp = double(h.dat.F.trace(h.dat.F.ichosen,:)) - double(h.dat.F.neurop(h.dat.F.ichosen,:));
            if h.show_unfiltered.Value
                plot(temp,'color',0.5*[1 1 1]);
                hold on
            end
            plot(my_conv_local(medfilt1(temp, 3), 3),'g');
            axis tight
        end
end
set(gca,'TickDir','out','TickLength',[0.002 0.005],'Fontsize',8);

% plot a separator if more than one recording session was concatenated
if size(h.dat.Fcell,2) > 1
    line_coord = 0;
    for i = 1:size(h.dat.Fcell,2)-1
        line_coord = line_coord + size(h.dat.Fcell{i},2);
        linehandle = line([line_coord line_coord],get(gca,'YLim'),'Color',0.5*[1 1 1],'LineStyle','--');
        uistack(linehandle,'bottom');
    end
end

if h.show_events.Value
    cellID = h.dat.cl.isroi(h.dat.F.ichosen)*...
        size(find(h.dat.cl.isroi(1:h.dat.F.ichosen)),2);
    if cellID
        axes(h.axes6)
        hold off
        set(gca,'box','off','XTick',[],'YTick',[]);
        my_spike_times = zeros(1,size(h.dat.F.trace,2));
        my_spike_times(1,h.dat.cl.dcell{cellID}.st) = h.dat.cl.dcell{cellID}.c;
        plot(my_spike_times,'r')
        axis tight
        set(gca,'XTick',[],'TickDir','out','TickLength',[0.002 0.005],'Fontsize',8);
        % plot a separator if more than one recording session was concatenated
        if size(h.dat.Fcell,2) > 1
            line_coord = 0;
            for i = 1:size(h.dat.Fcell,2)-1
                line_coord = line_coord + size(h.dat.Fcell{i},2);
                linehandle = line([line_coord line_coord],get(gca,'YLim'),'Color',0.5*[1 1 1],'LineStyle','--');
                uistack(linehandle,'bottom');
            end
        end
    end
end
end
