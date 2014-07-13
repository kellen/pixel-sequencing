function [  ] = graph_threshold( xvals, num_objects, precision, labeltext, filename )
    fh = figure('Visible','off');
    
    C = linspecer(6,'qualitative');
    C = C(2:4,:);
    set(gcf,'Color','black')
    line(xvals, num_objects,'Color',C(1,:),'LineWidth',3)
    haxes1 = gca; % handle to axes
    set(haxes1,'XColor','white', ...
        'YColor',C(1,:),'Color','none','LineWidth',2)
    set(get(haxes1,'YLabel'),'String','# objects');
    set(haxes1, 'YTickLabel', num2str(get(haxes1,'YTick')','%d'))
    haxes1_pos = get(haxes1,'Position'); 
    haxes2 = axes('Position',haxes1_pos,...
                  'YAxisLocation','right',...
                  'XColor','white','YColor',C(2,:),'Color','none','LineWidth',2);
    set(get(haxes2,'YLabel'),'String','precision');
    ylim(haxes2, [0 1]);
    set(get(haxes2,'XLabel'),'String',[labeltext]);
    line(xvals, precision,'Parent',haxes2,'Color',C(2,:),'LineWidth',3);
    
    set(fh, 'InvertHardCopy', 'off');
    
    sz = [0 0 7.5 6.25];
    %sz = [0 0 3.75 3.125];
    sz = [0 0 5.625 4.6875];
    set(fh,'PaperUnits','inches','PaperPosition',sz)
    print(fh,'-dtiff', '-r100', filename);
    %saveas(fh,filename,'tif')
end

