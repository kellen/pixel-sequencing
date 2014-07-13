function [] = filtered(input_dir, do_input_dir, output_dir) 
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
    
    for posidx = 1:numel(positions)
        position = positions(posidx);
        
        S = load(sprintf(fmt, input_dir, sprintf(config('mat_sequence'), position)));
        dofn = sprintf(fmt, do_input_dir, sprintf(config('img_do'), position));
        do = im2bw(imfilter(imread(dofn), fspecial('gaussian')), 0.10);
        
        [seq,errors] = apply(S, do, size_upper_threshold, size_lower_threshold, ...
            avg_threshold, int_threshold, quality_threshold );
        
        fn = sprintf(fmt, output_dir, sprintf(config('img_filtered'), position));
        write_image(label2rgb(seq, 'lines', 'k'), fn);
        
        fn = sprintf(fmt, output_dir, sprintf(config('mat_sequence'), position));
        mkdir_basename(fn);
        save(fn, 'seq', 'errors', 'avg_threshold', 'int_threshold', 'quality_threshold', 'size_lower_threshold', 'size_upper_threshold');
    end
end

function [seq, errors] = apply(S, do, size_upper_threshold, size_lower_threshold, ...
            avg_threshold, int_threshold, quality_threshold)
    config = get_config();
    valid = config('valid');
        ignore_apriori = (S.seq == 1111) | (S.seq == 2222) | (S.seq == 3333) | (S.seq == 4444);
        ignore_avg = (S.avg < avg_threshold);
        ignore_int = (max(S.int, [], 3) < int_threshold);
        ignore_quality = (S.quality < quality_threshold);

        filtered_size_big = filter_size(S.seq, size_upper_threshold);
        ignore_size_big = (filtered_size_big ~= 0);
        
        ignore = ignore_apriori | ignore_avg | ignore_int | ignore_quality | ignore_size_big;
        seq = S.seq;
        seq(~do) = 0;
        seq(ignore) = 0;
        [seq] = filter_size(seq, size_lower_threshold);
        
        errors = ~ismember(seq, valid);
        errors(seq == 0) = 0;
end

