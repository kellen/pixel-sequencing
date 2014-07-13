
% for presentation



show_thresholds('09-sequence', 1, '12-thresholds');

show_registration('06-mips-full', '07-registered', '04-reorganized', 1, 'XX-presentation');
write_image(roi(imread('04-reorganized/slideA/1/1/G.tif')), 'XX-presentation/tophat-before.tif');
write_image(roi(imread('08-tophats-registered/slideA/1/1/G.tif')), 'XX-presentation/tophat-after.tif');

config = get_config();
S = load(sprintf('%s/%s', '09-sequence', sprintf(config('mat_sequence'), 1)));

C = linspecer(size(unique(S.seq), 1), 'sequential');
L = uint8(zeros(size(S.seq)));
labels = unique(S.seq);
% ugh this is wrong, but whatever.
for labidx=1:numel(labels)
    % 2d mask
    mask = S.seq == labels(labidx);
    L(mask) = labidx;
end
seqim = ind2rgb(L,C);





write_image(roi(label2rgb(S.seq, 'lines', 'k')), 'XX-presentation/seq.tif');


write_image(roi(uint8(S.avg)),  'XX-presentation/average.tif');
write_image(roi(uint8(double(label2rgb(S.seq, 'lines', 'k')) .* repmat((S.avg / 255), [1 1 3]))), 'XX-presentation/average-seq.tif');


write_image(roi(uint8(max(S.int, [], 3))),  'XX-presentation/maxintensity.tif');

write_image(roi(uint8(double(label2rgb(S.seq, 'lines', 'k')) .* repmat((max(S.int, [], 3) / 255), [1 1 3]))), 'XX-presentation/maxintensity-seq.tif');


write_image(roi(S.quality),  'XX-presentation/quality.tif');
write_image(roi(uint8(double(label2rgb(S.seq, 'lines', 'k')) .* repmat(S.quality, [1 1 3]))), 'XX-presentation/quality-seq.tif');



ignore_apriori = (S.seq == 1111) | (S.seq == 2222) | (S.seq == 3333) | (S.seq == 4444);
write_image(roi(ignore_apriori), 'XX-presentation/ignore-apriori.tif');
s = S.seq;
s(ignore_apriori) = 0;
write_image(roi(label2rgb(s, 'lines', 'k')), 'XX-presentation/ignore-apriori-seq.tif');


athresh = 40;
ignore = S.avg < athresh;
write_image(roi(ignore), 'XX-presentation/ignore-average.tif');
s = S.seq;
s(ignore) = 0;
write_image(roi(label2rgb(s, 'lines', 'k')), 'XX-presentation/ignore-average-seq.tif');


%
% how much each thing says to ignore
%

    config = get_config();
    position = 1;
    fmt = '%s/%s/';
    dofn = sprintf(fmt, do_input_dir, sprintf(config('img_do'), position));
    do = im2bw(imfilter(imread(dofn), fspecial('gaussian')), 0.10);

    ignore_apriori = (S.seq == 1111) | (S.seq == 2222) | (S.seq == 3333) | (S.seq == 4444);
    ignore_size = (filtered_size == 0);
    ignore_size_big = (filtered_size_big ~= 0);
    ignore_avg = (S.avg < avg_threshold);
    ignore_int = (max(S.int, [], 3) < int_threshold);
    ignore_quality = (S.quality < quality_threshold);

    ignore = logical(zeros(size(S.seq)));
    ignore(:,:,1) = ignore_apriori;
    ignore(:,:,2) = ignore_size;
    ignore(:,:,3) = ignore_avg;
    ignore(:,:,4) = ignore_int;
    ignore(:,:,5) = ignore_quality;
    ignore(:,:,6) = ignore_size_big;
    avg_ignore = mean(ignore, 3);
    
    fnpatt = sprintf(fmt, fn_path, '/ignore.tif');
    write_image(avg_ignore, fnpatt);
    figure,imshow(avg_ignore,[]);

%
% graphs
lowx = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];
lowp = [0.19387 0.19387 0.22291 0.23042 0.26359 0.31546 0.37483 0.43608 0.49103 0.54183 0.58363 0.62049 0.65186 0.67375 0.69777 0.72203 0.73756 0.74885 0.76546 0.772 0.78188 0.78553 0.7956 0.79357 0.79402 0.78947 0.78039 0.77262 0.76654 0.76629 0.76959];
lown = [21659 21659 17518 16704 13521 10461 8153 6563 5519 4721 4155 3718 3355 3059 2776 2547 2351 2174 2021 1886 1733 1590 1453 1337 1204 1083 979 884 801 706 638 ];
graph_threshold( lowx, lown, lowp, 'lower object size threshold (pixels)', 'XX-presentation/graph-low.tif');

upn = [0 17389 19886 21011 21424 21552 21600 21632 21641 21655 21655 21655 21655 21659 21659];
upp = [0 0.10093 0.14221 0.17638 0.18867 0.19223 0.19347 0.19397 0.19403 0.1939 0.1939 0.1939 0.1939 0.19387 0.19387];
upx = [0 10 20 30 40 50 60 70 80 90 100 110 120 130 140];
graph_threshold( upx, upn, upp, 'upper object size threshold (pixels)', 'XX-presentation/graph-upper.tif');

avgx = [0 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 ];
avgp = [0.19387 0.28922 0.39475 0.505 0.60801 0.68711 0.75652 0.80482 0.83956 0.87411 0.89682 0.90821 0.92169 0.93453 0.94568 0.9531 0.963 0.96817 0.97344 0.97406 0.98045 ];
avgn = [21659 11939 7924 5697 4393 3570 3027 2654 2356 2105 1919 1754 1609 1451 1307 1194 1081 974 866 771 665 ];
graph_threshold( avgx, avgn, avgp, 'average sequence intensity threshold', 'XX-presentation/graph-avg.tif');

maxn = [21659 9437 5370 3841 3032 2554 2217 1922 1670 1468 1235 1021 780 577 391 218 93 ];
maxp = [0.19387 0.34937 0.54581 0.68654 0.77836 0.83359 0.86919 0.89698 0.91617 0.92847 0.94332 0.9569 0.96282 0.974 0.97442 0.97706 0.95699 ];
maxx = [0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 ];
graph_threshold( maxx, maxn, maxp, 'maximum sequence intensity threshold', 'XX-presentation/graph-max.tif');

qtyx = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 ];
qtyp = [0.19387 0.19387 0.19387 0.20852 0.33584 0.44505 0.68122 0.74286 0.71205 0.53356 0.26056 ];
qtyn = [21659 21659 21659 19682 10886 7179 3620 2380 1469 581 284 ];
graph_threshold( qtyx, qtyn, qtyp, 'quality threshold', 'XX-presentation/graph-qual.tif');


%posx = [1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16];
posx = [1     2     3     4     6     7     8     9    10    11    15    16];

%cpp = [0.97167 0.98234 0.9897  0.97622 0 0.97931 0.98208 0.98736 0.98501 0.97693 0.94904 0 0 0 0.93028 0.97282];
%ppp = [0.97167 0.98101 0.98796 0.97013 0 0.97703 0.98446 0.98921 0.98024 0.98282 0.93902 0 0 0 0.96588 0.96676];
cpp = [0.97167 0.98234 0.9897  0.97622 0.97931 0.98208 0.98736 0.98501 0.97693 0.94904 0.93028 0.97282];
ppp = [0.97167 0.98101 0.98796 0.97013 0.97703 0.98446 0.98921 0.98024 0.98282 0.93902 0.96588 0.96676];

%cpvalid = [1818 3783 3172 780 0 2367 2028 1953 1905 720 298 0 0 0 507 3007];
%ppvalid = [1818 3513 2871 617 0 2042 1900 1742 1637 572 231 0 0 0 368 2647];
cpvalid = [1818 3783 3172 780 2367 2028 1953 1905 720 298 507 3007];
ppvalid = [1818 3513 2871 617 2042 1900 1742 1637 572 231 368 2647];

%cperror = [53 68 33 19 0 50 37 25 29 17 16 0 0 0 38 84];
%pperror = [53 68 35 19 0 48 30 19 33 10 15 0 0 0 13 91];
cperror = [53 68 33 19 50 37 25 29 17 16 38 84];
pperror = [53 68 35 19 48 30 19 33 10 15 13 91];


lw = 2;

    fh = figure; %('Visible','off');
    C = linspecer(12);
    set(gcf,'Color','black')
    haxes1 = gca; % handle to axes
    line(posx, cpvalid, 'Parent',haxes1, 'Color',C(1,:),'LineWidth',lw);
    line(posx, ppvalid, 'Parent',haxes1,'Color',C(2,:),'LineWidth',lw);
    line(posx, cperror, 'Parent',haxes1,'Color',C(10,:),'LineWidth',lw);
     line(posx, pperror, 'Parent',haxes1,'Color',C(9,:),'LineWidth',lw);
    set(haxes1,'XColor','white', 'YColor',C(1,:),'Color','none','LineWidth',2)
    set(get(haxes1,'YLabel'),'String','# objects');
    %set(haxes1, 'YTickLabel', num2str(get(haxes1,'YTick')','%d'))
    set(get(haxes1,'XLabel'),'String','position');
    xlim(haxes1,[1 16])
    ylim(haxes1,[0 max(cpvalid)]);
    set(fh, 'InvertHardCopy', 'off');
    sz = [0 0 5.625 4.6875];
    set(fh,'PaperUnits','inches','PaperPosition',sz)
    print(fh,'-dtiff', '-r100', 'XX-presentation/graph-comparision-objs.tif');
    
    
    fh = figure; %('Visible','off');
    C = linspecer(12);
    set(gcf,'Color','black')
    haxes1 = gca; % handle to axes
    line(posx, cpp, 'Parent',haxes1, 'Color',C(1,:),'LineWidth',lw);
    line(posx, ppp, 'Parent',haxes1,'Color',C(2,:),'LineWidth',lw);
    set(haxes1,'XColor','white', 'YColor',C(1,:),'Color','none','LineWidth',2)
    set(get(haxes1,'YLabel'),'String','precision');
    %set(haxes1, 'YTickLabel', num2str(get(haxes1,'YTick')','%f'))
    set(get(haxes1,'XLabel'),'String','position');
    xlim(haxes1,[1 16])
    ylim(haxes1, [0 1]);
    set(fh, 'InvertHardCopy', 'off');
    sz = [0 0 5.625 4.6875];
    set(fh,'PaperUnits','inches','PaperPosition',sz)
    print(fh,'-dtiff', '-r100', 'XX-presentation/graph-comparison-precision.tif');
    
    
    
    
    haxes1_pos = get(haxes1,'Position'); 
    haxes2 = axes('Position',haxes1_pos,...
                  'YAxisLocation','right',...
                  'XColor','white','YColor',C(2,:),'Color','none','LineWidth',2);
    set(get(haxes2,'YLabel'),'String','precision');
    ylim(haxes2, [0 1]);
    
    

    line(posx, cpp,'Parent',haxes2,'Color',C(5,:),'LineWidth',lw);
    line(posx, ppp,'Parent',haxes2,'Color',C(6,:),'LineWidth',lw);

    
    set(fh, 'InvertHardCopy', 'off');
    
    sz = [0 0 7.5 6.25];
    %sz = [0 0 3.75 3.125];
    sz = [0 0 5.625 4.6875];
    set(fh,'PaperUnits','inches','PaperPosition',sz)
    print(fh,'-dtiff', '-r100', filename);
    %saveas(fh,filename,'tif')





