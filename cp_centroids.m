function [blobpositions, bloblabels, quality] = cp_centroids(input_dir, position)
    config = get_config();
    cycles = config('cycles');
    types = config('types');
    key = config('label_cellprofiler');
    
    S = load(sprintf('%s/%s', input_dir, config('mat_cellprofiler')));
    m = S.handles.Measurements.blobs;
    
    % position data exists in these for each image; ignore all but the first
    % since these are the same in every image.
    blobpositions = cat(2, ...
        m.Location_Center_X{position}, ...
        m.Location_Center_Y{position});
    % shift to matlab-based positions (i.e start at 1 rather than 0)
    blobpositions = blobpositions + 1;

    intensity = cell([1 numel(cycles)]);
    for cycleidx=1:numel(cycles)
        cycle = cycles(cycleidx);
        intensity{cycleidx} = [];
        
        for typeidx=1:numel(types)  
            intensity{cycleidx} = cat(2, ...
                intensity{cycleidx}, ...
                m.(sprintf(key, cycle, typeidx)){position} ...
            );
        end
    end
    
    % accumulate the positions into the sequence format
    bloblabels = zeros([size(blobpositions, 1) 1]);
    quality = inf([size(blobpositions, 1) 1]);
    shift = 10 ^ (numel(cycles) - 1);
    for i=1:numel(cycles)
        total = sum(intensity{i}, 2);
        [maxes, indexes] = max(intensity{i}, [], 2);
        bloblabels = bloblabels + (indexes .* shift);
        shift = shift / 10;
        
        quality = min(double(maxes) ./ total, quality);
    end
end