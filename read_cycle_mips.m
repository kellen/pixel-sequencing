function [mips] = read_cycle_mips(input_dir, cycle)
    config = get_config();
    positions = config('positions');
    
    mips = cell([numel(positions) 1]);
    for posidx = 1:numel(positions) 
        position = positions(posidx);
        fn = sprintf('%s/%s', input_dir, sprintf(config('img_mip'), position, cycle));
        mips{posidx} = imread(fn);
    end
end