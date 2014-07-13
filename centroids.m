function [] = centroids(input_dir, output_dir) 
    config = get_config();
    write_config(output_dir);
    positions = config('positions'); 
    fmt = '%s/%s';
    
    for posidx = 1:numel(positions)
        position = positions(posidx);
        S = load(sprintf(fmt, input_dir, sprintf(config('mat_sequence'), position)));
        seq = S.seq;
        
        % get the list of labels, excluding 0 which belongs to no class
        labels = unique(seq);
        labels = labels(labels ~= 0);

        centpositions = [];
        centlabels = [];
        centareas = [];
        
        se = strel('disk', 1);
        for i=1:numel(labels);
            % get a binary image of the current label
            bw = seq == labels(i);
            % dilate the binary image to retain single pixels after watershed
            bw = imdilate(bw, se);
            % do watershed
            L = wshed(bw);
            % get the centroids
            s = regionprops(L, 'Centroid', 'Area');
            areas = cat(1, s.Area);
            cents = cat(1, s.Centroid);
            % wshed() removes the background label for us, but we have to
            % get rid of the centroid for this now non-existent label
            % which registers as a NaN
            if ~isempty(cents)
                cents = cents(isfinite(cents(:, 1)), :);
                areas = areas(isfinite(cents(:, 1)), :);

                if ~isempty(cents)
                    centpositions = cat(1, centpositions, cents);
                    centlabels = cat(1, centlabels, labels(i) * ones([size(cents, 1) 1]));
                else
                    disp(['Removed all instances of class ' num2str(labels(i))]);
                end
            end
        end
        
        fn = sprintf(fmt, output_dir, sprintf(config('mat_centroids'), position));
        mkdir_basename(fn);
        save(fn, 'centpositions', 'centlabels', 'centareas');
    end
end