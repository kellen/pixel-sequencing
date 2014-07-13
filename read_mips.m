function [mips] = read_mips(input_dir, position) 
    config = get_config();
    cycles = config('cycles');
    
    mips = cell([numel(cycles) 1]);
    for cycleidx = 1:numel(cycles)
        cycle = cycles(cycleidx);
        fn = sprintf('%s/%s', input_dir, sprintf(config('img_mip'), position, cycle));
        mips{cycleidx} = imread(fn);
    end
end