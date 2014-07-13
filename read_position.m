function [pos] = read_position(input_dir, position)
    config = get_config();
    types = config('types');
    cycles = config('cycles');
    
    pos = containers.Map;
    c = cell([numel(cycles) 1]);
    fmt = '%s/%s';
    for cycleidx=1:numel(cycles)
        cycle = cycles(cycleidx);
        for typeidx=1:numel(types)
            type = types{typeidx};
            fn = sprintf(fmt, input_dir, sprintf(config('img_cycle'), position, cycle, type));
            c{cycleidx}(:,:,typeidx) = imread(fn);
        end
    end
    pos('cycles') = c;
    fn = sprintf(fmt, input_dir, sprintf(config('img_do'), position));
    pos('do') = imread(fn);
end