function [] = sequence(input_dir, output_dir)
    config = get_config();
    write_config(output_dir);
    positions = config('positions');
    valid = config('valid');
    
    for posidx = 1:numel(positions)
        position = positions(posidx);
        pos = read_position(input_dir, position);
        
        images = pos('cycles');
        
        [x, y, z] = size(images{1});
    
        % seq is the computed "best" sequence for a given pixel
        seq = zeros([x y]);

        % quality is the minimum quality for all bases in the sequence
        % where base quality is the % of intensity values for all channels
        quality = inf([x y]);

        % the average intensity value for the entire "best" sequence
        avg = zeros([x y]);

        % the actual intensity values for the entire "best" sequence
        int = zeros([x y z]);

        % accumulates the positions into the same format as ID_list_BCpanel
        shift = 10 ^ (numel(images) - 1);
        for i=1:numel(images)
            total = sum(images{i}, 3);
            [maxes, indexes] = max(images{i}, [], 3);

            int(:, :, i) = maxes;
            avg = avg + double(maxes);
            seq = seq + (indexes .* shift);
            quality = min(double(maxes) ./ total, quality);

            shift = shift / 10;
        end
        % set pixels with zero intensity in all channels to 0
        quality(quality == Inf) = 0;
        avg = avg ./ numel(images);
        errors = ~ismember(seq, valid);

        fn = sprintf('%s/%s', output_dir, sprintf(config('mat_sequence'), position));
        mkdir_basename(fn);
        save(fn, 'seq', 'quality', 'avg', 'int', 'errors');
    end
end




