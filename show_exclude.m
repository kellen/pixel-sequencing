function [ ] = show_exclude( input_dir, do_input_dir, output_dir, position )
    config = get_config();
    write_config(output_dir);
    positions = config('positions');
    key = config('thresholds_key');
    map = config(key);
    avg_threshold = map('avg');
    int_threshold = map('int');
    quality_threshold = map('quality');
    
    size_lower_threshold = map('size_lower');
    size_upper_threshold = map('size_upper');
    
    fmt = '%s/%s';
    
        
    S = load(sprintf(fmt, input_dir, sprintf(config('mat_sequence'), position)));
    dofn = sprintf(fmt, do_input_dir, sprintf(config('img_do'), position));
    do = im2bw(imfilter(imread(dofn), fspecial('gaussian')), 0.10);

    ignore_apriori = (S.seq == 1111) | (S.seq == 2222) | (S.seq == 3333) | (S.seq == 4444);
    ignore_avg = (S.avg < avg_threshold);
    ignore_int = (max(S.int, [], 3) < int_threshold);
    ignore_quality = (S.quality < quality_threshold);
    ignore_size_big = (filter_size(S.seq, size_upper_threshold) ~= 0);
    ignore_size_small = (filter_size(S.seq, size_lower_threshold) == 0);
    
    ignore = logical(zeros(size(S.seq)));
    ignore(:,:,1) = ignore_apriori;
    ignore(:,:,2) = ignore_size_small;
    ignore(:,:,3) = ignore_avg;
    ignore(:,:,4) = ignore_int;
    ignore(:,:,5) = ignore_quality;
    ignore(:,:,6) = ignore_size_big;
    avg_ignore = mean(ignore, 3);
    
    write_image(roi(avg_ignore), sprintf(fmt, output_dir, '/exclude-average-contribution.tif'));
    figure,imshow(avg_ignore,[]);
end

