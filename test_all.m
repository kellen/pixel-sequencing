function [] = test_all(input_dir, position)
    config = get_config();
    
    S = load(sprintf('%s/%s', input_dir, sprintf(config('mat_sequence'), position)));
    
    ignore_apriori = (S.seq == 1111) | (S.seq == 2222) | (S.seq == 3333) | (S.seq == 4444);
    
    do = im2bw(imread('07-registered\slideA\1\do.tif'));
    
    
    do = im2bw(imfilter(imread('07-registered\slideA\1\do.tif'), fspecial('gaussian')), 0.10);
    uthresh = 80;
    lthresh = 5;
    athresh = 25; %44;
    ithresh = 44; % 44
    
    qthresh = 0.1;
    s = apply(S, do, uthresh, lthresh, athresh, ithresh, qthresh);
    figure,imshow(label2rgb(s, 'lines', 'k'));

    keyboard;
    
        uthresh = 80;
    lthresh = 5;
    athresh = 25; %44;
    ithresh = 40; % 44
    for qthresh = 0:0.05:1.0
        s = apply(S, do, uthresh, lthresh, athresh, ithresh, qthresh);
        figure,imshow(label2rgb(s, 'lines', 'k'));
    end
    
    keyboard;
    
    
    
    %{ 
    % ideal?
    uthresh = 60;
    lthresh = 5;
    athresh = 16;
    ithresh = 20;
    qthresh = 0.4;

    
    [blobpositions, bloblabels, blobquality] = cp_centroids('13-cp-output', 1);
    valid = config('valid');
    for qthresh = 0.4:0.05:0.6
        bl = bloblabels(blobquality > qthresh);
        valid_blob = ismember(bl, valid);
        blob_precision = sum(valid_blob) / size(bl, 1);
        disp(['qthresh: ' num2str(qthresh) ', num: ' num2str(size(bl,1)) ', precision: ' num2str(blob_precision)]);
    end
    
    uthresh = 60;
    lthresh = 5;
    athresh = 40;
    ithresh = 50;
    qthresh = 0.55;
    for qthresh = 0.4:0.05:0.6
        for lthresh=0:10
            disp(['qthresh: ' num2str(qthresh) ', lthresh: ' num2str(lthresh)]);
            apply(S, do, uthresh, lthresh, athresh, ithresh, qthresh );
        end
    end
    
    
    
    % conservative thresholds
    uthresh = 100;
    lthresh = 5;
    athresh = 16;
    ithresh = 20;
    qthresh = 0.4;
    disp('Conservative');
    summarize(S, do, uthresh, lthresh, athresh, ithresh, qthresh );


    % aggressive thresholds
    uthresh = 60;
    lthresh = 9;
    athresh = 30;
    ithresh = 50;
    qthresh = 0.55;
    disp('Aggressive');
    summarize(S, do, uthresh, lthresh, athresh, ithresh, qthresh );

        
    s = filter_size(s, lthresh);
    figure,imshow(label2rgb(s, 'lines', 'k'));
    [p, n] = calc_precision(s);
    
    keyboard;
        %}
    
    %{
    for uthresh = config('show_steps_size_upper')
        s = S.seq;
        im = filter_size(s, uthresh);
        rm = ~(im == 0);
        s(ignore_apriori) = 0;
        s(rm) = 0;
        [p, n] = calc_precision(s);
        disp([ ...
            'Num: ' num2str(n) ...
            ' Precision: ' num2str(p)  ...
            ' Upp: ' num2str(uthresh) ...
            ]);
    end
    %}
    
    for lthresh = config('show_steps_size')
        s = S.seq;
        s = filter_size(s, lthresh);
        s(ignore_apriori) = 0;
        [p, n] = calc_precision(s);
        disp([ ...
            'Num: ' num2str(n) ...
            ' Precision: ' num2str(p)  ...
            ' Low: ' num2str(lthresh) ...
            ]);
    end
    
    keyboard;
    
    for athresh = config('show_steps_average')
        s = S.seq;
        ignore = (S.avg < athresh);
        s(ignore) = 0;
        s(ignore_apriori) = 0;
        [p, n] = calc_precision(s);
        disp([ ...
            'Num: ' num2str(n) ...
            ' Precision: ' num2str(p)  ...
            ' Avg: ' num2str(athresh) ...
            ]);
    end

    for ithresh = config('show_steps_intensity')
        s = S.seq;
        ignore = (max(S.int, [], 3) < ithresh);
        s(ignore) = 0;
        s(ignore_apriori) = 0;
        [p, n] = calc_precision(s);
        disp([ ...
            'Num: ' num2str(n) ...
            ' Precision: ' num2str(p)  ...
            ' Max: ' num2str(ithresh) ...
            ]);
    end
    
    for qthresh = config('show_steps_quality')
        s = S.seq;
        ignore = (S.quality < qthresh);
        s(ignore) = 0;
        s(ignore_apriori) = 0;
        [p, n] = calc_precision(s);
        disp([ ...
            'Num: ' num2str(n) ...
            ' Precision: ' num2str(p)  ...
            ' Qty: ' num2str(qthresh) ...
            ]);
    end
    
    
end

function [s] = apply(S, do, uthresh, lthresh, athresh, ithresh, qthresh )
    s = S.seq;
    ignore_apriori = (s == 1111) | (s == 2222) | (s == 3333) | (s == 4444);

    filtered_size_big = filter_size(s, uthresh);
    ignore_size_big = (filtered_size_big ~= 0);

    ignore_avg = (S.avg < athresh);
    ignore_int = (max(S.int, [], 3) < ithresh);
    ignore_quality = (S.quality < qthresh);

    s(ignore_apriori) = 0;
    s(~do) = 0;
    s(ignore_size_big) = 0;
    s(ignore_quality) = 0;
    s(ignore_avg) = 0;
    s(ignore_int) = 0;
    s = filter_size(s, lthresh);
    [p, n] = calc_precision(s);
    disp(['Final, num: ' num2str(n) ', precision: ' num2str(p)]);
end

function [] = summarize(S, do, uthresh, lthresh, athresh, ithresh, qthresh )
    s = S.seq;
    ignore_apriori = (s == 1111) | (s == 2222) | (s == 3333) | (s == 4444);

    filtered_size_big = filter_size(s, uthresh);
    ignore_size_big = (filtered_size_big ~= 0);
    
    ignore_avg = (S.avg < athresh);
    ignore_int = (max(S.int, [], 3) < ithresh);
    ignore_quality = (S.quality < qthresh);

    s(ignore_apriori) = 0;
    [p, n] = calc_precision(s);
    disp(['Apriori, num: ' num2str(n) ', precision: ' num2str(p)]);
    
    s(~do) = 0;
    [p, n] = calc_precision(s);
    disp(['DO, num: ' num2str(n) ', precision: ' num2str(p)]);
    
    s(ignore_size_big) = 0;
    [p, n] = calc_precision(s);
    disp(['Big, num: ' num2str(n) ', precision: ' num2str(p)]);
    
    s(ignore_quality) = 0;
    [p, n] = calc_precision(s);
    disp(['Quality, num: ' num2str(n) ', precision: ' num2str(p)]);

    s(ignore_avg) = 0;
    [p, n] = calc_precision(s);
    disp(['Avg intensity, num: ' num2str(n) ', precision: ' num2str(p)]);
    
    s(ignore_int) = 0;
    [p, n] = calc_precision(s);
    disp(['Max intensity, num: ' num2str(n) ', precision: ' num2str(p)]);
    
    s = filter_size(s, lthresh);
    [p, n] = calc_precision(s);
    disp(['Small, num: ' num2str(n) ', precision: ' num2str(p)]);

end


function [cent_precision,numcents] = calc_precision(seq)
    config = get_config();
    % get the list of labels, excluding 0 which belongs to no class
    labels = unique(seq);
    labels = labels(labels ~= 0);
    
    valid = config('valid');

    centpositions = [];
    centlabels = [];

    se = strel('disk', 1);
    for i=1:numel(labels);
        % get a binary image of the current label
        bw = seq == labels(i);
        % dilate the binary image to retain single pixels after watershed
        bw = imdilate(bw, se);
        % do watershed
        L = wshed(bw);
        % get the centroids
        s = regionprops(L, 'Centroid');
        cents = cat(1, s.Centroid);
        
        if ~isempty(cents)
            % wshed() removes the background label for us, but we have to
            % get rid of the centroid for this now non-existent label
            % which registers as a NaN
            cents = cents(isfinite(cents(:, 1)), :);
            if ~isempty(cents)
                centpositions = cat(1, centpositions, cents);
                centlabels = cat(1, centlabels, labels(i) * ones([size(cents, 1) 1]));
            end
        %else
            %disp(['Removed all instances of class ' num2str(labels(i))]);
        end
    end
    
    valid_cent = ismember(centlabels, valid);
    cent_precision = sum(valid_cent) / size(centlabels, 1);
    numcents = size(centlabels,1);
end


%{
test_all('09-sequence', 1)
Num: 17389      1 Precision: 0.10093 Upp: 10
Num: 19886      1 Precision: 0.14221 Upp: 20
Num: 21011      1 Precision: 0.17638 Upp: 30
Num: 21424      1 Precision: 0.18867 Upp: 40
Num: 21552      1 Precision: 0.19223 Upp: 50
Num: 21600      1 Precision: 0.19347 Upp: 60
Num: 21632      1 Precision: 0.19397 Upp: 70
Num: 21641      1 Precision: 0.19403 Upp: 80
Num: 21655      1 Precision: 0.1939 Upp: 90
Num: 21655      1 Precision: 0.1939 Upp: 100
Num: 21655 Precision: 0.1939 Upp: 110
Num: 21655 Precision: 0.1939 Upp: 120
Num: 21659 Precision: 0.19387 Upp: 130
Num: 21659 Precision: 0.19387 Upp: 140

Num: 21659      1 Precision: 0.19387 Low: 0
Num: 21659      1 Precision: 0.19387 Low: 1
Num: 17518      1 Precision: 0.22291 Low: 2
Num: 16704      1 Precision: 0.23042 Low: 3
Num: 13521      1 Precision: 0.26359 Low: 4
Num: 10461      1 Precision: 0.31546 Low: 5
Num: 8153     1 Precision: 0.37483 Low: 6
Num: 6563     1 Precision: 0.43608 Low: 7
Num: 5519     1 Precision: 0.49103 Low: 8
Num: 4721     1 Precision: 0.54183 Low: 9
Num: 4155     1 Precision: 0.58363 Low: 10
Num: 3718     1 Precision: 0.62049 Low: 11
Num: 3355     1 Precision: 0.65186 Low: 12
Num: 3059     1 Precision: 0.67375 Low: 13
Num: 2776     1 Precision: 0.69777 Low: 14
Num: 2547     1 Precision: 0.72203 Low: 15
Num: 2351 Precision: 0.73756 Low: 16
Num: 2174 Precision: 0.74885 Low: 17
Num: 2021 Precision: 0.76546 Low: 18
Num: 1886 Precision: 0.772 Low: 19
Num: 1733 Precision: 0.78188 Low: 20
Num: 1590 Precision: 0.78553 Low: 21
Num: 1453 Precision: 0.7956 Low: 22
Num: 1337 Precision: 0.79357 Low: 23
Num: 1204 Precision: 0.79402 Low: 24
Num: 1083 Precision: 0.78947 Low: 25
Num: 979 Precision: 0.78039 Low: 26
Num: 884 Precision: 0.77262 Low: 27
Num: 801 Precision: 0.76654 Low: 28
Num: 706 Precision: 0.76629 Low: 29
Num: 638 Precision: 0.76959 Low: 30

Num: 21659      1 Precision: 0.19387 Avg: 0
Num: 11939      1 Precision: 0.28922 Avg: 4
Num: 7924     1 Precision: 0.39475 Avg: 8
Num: 5697     1 Precision: 0.505 Avg: 12
Num: 4393     1 Precision: 0.60801 Avg: 16
Num: 3570     1 Precision: 0.68711 Avg: 20
Num: 3027     1 Precision: 0.75652 Avg: 24
Num: 2654     1 Precision: 0.80482 Avg: 28
Num: 2356     1 Precision: 0.83956 Avg: 32
Num: 2105     1 Precision: 0.87411 Avg: 36
Num: 1919     1 Precision: 0.89682 Avg: 40
Num: 1754     1 Precision: 0.90821 Avg: 44
Num: 1609     1 Precision: 0.92169 Avg: 48
Num: 1451     1 Precision: 0.93453 Avg: 52
Num: 1307     1 Precision: 0.94568 Avg: 56
Num: 1194     1 Precision: 0.9531 Avg: 60
Num: 1081     1 Precision: 0.963 Avg: 64
Num: 974    1 Precision: 0.96817 Avg: 68
Num: 866    1 Precision: 0.97344 Avg: 72
Num: 771    1 Precision: 0.97406 Avg: 76
Num: 665    1 Precision: 0.98045 Avg: 80

Num: 21659      1 Precision: 0.19387 Max: 0
Num: 9437     1 Precision: 0.34937 Max: 10
Num: 5370     1 Precision: 0.54581 Max: 20
Num: 3841     1 Precision: 0.68654 Max: 30
Num: 3032     1 Precision: 0.77836 Max: 40
Num: 2554     1 Precision: 0.83359 Max: 50
Num: 2217     1 Precision: 0.86919 Max: 60
Num: 1922     1 Precision: 0.89698 Max: 70
Num: 1670     1 Precision: 0.91617 Max: 80
Num: 1468     1 Precision: 0.92847 Max: 90
Num: 1235     1 Precision: 0.94332 Max: 100
Num: 1021     1 Precision: 0.9569 Max: 110
Num: 780    1 Precision: 0.96282 Max: 120
Num: 577    1 Precision: 0.974 Max: 130
Num: 391    1 Precision: 0.97442 Max: 140
Num: 218    1 Precision: 0.97706 Max: 150
Num: 93   1 Precision: 0.95699 Max: 160

Num: 21659      1 Precision: 0.19387 Qty: 0
Num: 21659      1 Precision: 0.19387 Qty: 0.1
Num: 21659      1 Precision: 0.19387 Qty: 0.2
Num: 19682      1 Precision: 0.20852 Qty: 0.3
Num: 10886      1 Precision: 0.33584 Qty: 0.4
Num: 7179     1 Precision: 0.44505 Qty: 0.5
Num: 3620     1 Precision: 0.68122 Qty: 0.6
Num: 2380     1 Precision: 0.74286 Qty: 0.7
Num: 1469     1 Precision: 0.71205 Qty: 0.8
Num: 581    1 Precision: 0.53356 Qty: 0.9
Num: 284    1 Precision: 0.26056 Qty: 1
%}