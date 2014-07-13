function [] = reorganize(input_dir, output_dir)

    config = get_config();
    write_config(output_dir);
    prefix = config('prefix');
    positions = config('positions');
    cycles = config('input_cycles');
    types = config('types');
    channels = config('channels');
    do = config('do_dirname');
    do_channel = config('do_channel');
    
    % FIXME use config instead...
    fmt = '%s/%s';
    fmt_indir = '%s_%d_%s';
    fmt_outdir = '%s/%d/%d';
    fmt_infile = '%s_%d_%s_c%d.tif';
    fmt_outfile = '%s.tif';
    fmt_do_outdir = '%s/%d';
    
    for posidx = 1:numel(positions)
        position = positions(posidx);
        % copy the DO
        indir = sprintf(fmt, input_dir, sprintf(fmt_indir, prefix, position, do));
        outdir = sprintf(fmt, output_dir, sprintf(fmt_do_outdir, prefix, position));
        [s,mess,~] = mkdir(outdir);
        if ~s
            error(['Could not create ' outdir ' Message:' mess]);
        end
        channel = do_channel;
        infile = sprintf(fmt, indir, sprintf(fmt_infile, prefix, position, do, channel));
        outfile = sprintf(fmt, outdir, sprintf(fmt_outfile, do));
        copyfile(infile, outfile);
        
        % copy each of the cycles
        for cycleidx = 1:numel(cycles)
            cycle = cycles{cycleidx};
            indir = sprintf(fmt, input_dir, sprintf(fmt_indir, prefix, position, cycle));
            outdir = sprintf(fmt, output_dir, sprintf(fmt_outdir, prefix, position, cycleidx));
            [s,mess,~] = mkdir(outdir);
            if ~s
                error(['Could not create ' outdir ' Message:' mess]);
            end
            for channelidx = 1:numel(channels) 
                channel = channels(channelidx);
                infile = sprintf(fmt, indir, sprintf(fmt_infile, prefix, position, cycle, channel));
                outfile = sprintf(fmt, outdir, sprintf(fmt_outfile, types{channelidx}));
                copyfile(infile, outfile);
            end
        end
    end
end