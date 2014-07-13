function [] = mip(input_dir, output_dir)
    config = get_config();
    write_config(output_dir);
    positions = config('positions');
    cycles = config('cycles');
    
    for posidx = 1:numel(positions)
        position = positions(posidx);
        pos = read_position(input_dir, position);
        cycleimgs = pos('cycles');
        for cycleidx=1:numel(cycles)
            cycle = cycles(cycleidx);
            i = cycleimgs{cycleidx};
            fn = sprintf('%s/%s', output_dir, sprintf(config('img_mip'), position, cycle));
            write_image(max(i, [], 3), fn);
        end
    end
end