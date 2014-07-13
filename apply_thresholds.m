function [s] = apply_thresholds(S, do, uthresh, lthresh, athresh, ithresh, qthresh )
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