function [] = tophat(input_dir, output_dir)
    config = get_config();
    write_config(output_dir);
    positions = config('positions');
    types = config('types');
    cycles = config('cycles');
    
    se = strel('disk', config('tophat_size'));
    fmt = '%s/%s';
    for posidx = 1:numel(positions)
        position = positions(posidx);
        pos = read_position(input_dir, position);
        images = pos('cycles');
        for cycleidx = 1:numel(images)
            cycle = cycles(cycleidx);
            image = images{cycleidx};
            for channel=1:size(image, 3)
                type = types{channel};
                fn = sprintf(fmt, output_dir, sprintf(config('img_cycle'), position, cycle, type));
                th = imtophat(image(:,:,channel), se);
                write_image(th,fn);
            end
        end
        % top hat the DO as well
        image = pos('do');
        th = imtophat(image, se);
        fn = sprintf(fmt, output_dir, sprintf(config('img_do'), position));
        write_image(th, fn);
    end
end