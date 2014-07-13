function [] = show_thresholds(input_dir, position, output_dir)
    config = get_config();
    key = config('thresholds_key');
    map = config(key);
    avg_threshold = map('avg');
    int_threshold = map('int');
    quality_threshold = map('quality');
    size_threshold = map('size');
    size_threshold_lower = 5;
    
    S = load(sprintf('%s/%s', input_dir, sprintf(config('mat_sequence'), position)));
    
    fmt = '%s/%s';
    fn_path = sprintf(fmt, output_dir, sprintf(config('dir_position'), position));
    
    disp('Average intensity:');
    fnpatt = sprintf(fmt, fn_path, 'avg/average-%d.tif');
    show_steps(S.seq, S.avg, config('show_steps_average'), fnpatt);
    
    disp('Quality:');
    fnpatt = sprintf(fmt, fn_path, 'quality/quality-%f.tif');
    show_steps(S.seq, S.quality, config('show_steps_quality'), fnpatt);
    
    disp('Max intensity:');
    fnpatt = sprintf(fmt, fn_path, 'int/int-%d.tif');
    show_steps(S.seq, max(S.int, [], 3), config('show_steps_intensity'), fnpatt);
    
    disp('Size lower:');
    fnpatt = sprintf(fmt, fn_path, 'size-lower/size-%d.tif');
    show_size_threshold(S.seq, fnpatt);
    
    disp('Size upper:');
    fnpatt = sprintf(fmt, fn_path, 'size-upper/size-%d.tif');
    show_size_threshold_upper(S.seq, fnpatt);
    
    filtered_size_big = filter_size(S.seq, size_threshold);
    filtered_size = filter_size(S.seq, size_threshold_lower);

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

    remove = max(ignore, [], 3);
    filtered_size(remove) = 0;
    
    fnpatt = sprintf(fmt, fn_path, '/final.tif');
    write_image(label2rgb(filtered_size, 'lines', 'k'), fnpatt);
end

function [] = show_size_threshold(seq, fnpatt)
    config = get_config();
    for threshold = config('show_steps_size')
        img = label2rgb(filter_size(seq, threshold), 'lines', 'k');
        write_image(img, sprintf(fnpatt, threshold));
    end
end

function [] = show_size_threshold_upper(seq, fnpatt)
    config = get_config();
    for threshold = config('show_steps_size_upper')
        im = filter_size(seq, threshold);
        rm = ~(im == 0);
        seq(rm) = 0;
        img = label2rgb(seq, 'lines', 'k');
        write_image(img, sprintf(fnpatt, threshold));
    end
end

function [] = show_steps(seq, criteria, steps, fnpatt)
    for threshold = steps
        ignore = (criteria < threshold);
        seq(ignore) = 0;
        
        img = label2rgb(seq, 'lines', 'k');
        write_image(img, sprintf(fnpatt, threshold));
    end
end

