function [ ] = write_config( output_dir )
    config = get_config();
    fn = sprintf('%s/%s', output_dir, config('config_output_filename'));
    mkdir_basename(fn);

    f = fopen(fn,'a');
    if f ~= -1
        printval(f, config);
        fclose(f);
    else
        error('Could not write config to output file [%s]', fn);
    end
end

function [] = printval(f, val)
    if isa(val, 'containers.Map')
        printmap(f, val);
    elseif isa(val, 'double')
        fprintf(f, '%f', val);
    elseif isa(val, 'char')
        fprintf(f, '''%s''', val);
    elseif isa(val, 'logical')
        if val
            fprintf(f, 'true');
        else
            fprintf(f, 'false');
        end
    elseif isa(val, 'cell')
        printcell(f,val);
    elseif isa(val, 'struct')
        printstruct(f,val);
    elseif ismatrix(val)
        printmatrix(f,val);
    else
        error('Cannot print type [%s]', class(val));
    end
end

function [] = printcell(f, val)
    fprintf(f, '{');
    % assume only 2d cells for config
    [r,~] = size(val);
    for j=1:r
        row = val(j,:);
        for k=1:numel(row)
            printval(f, val{j,k});
        end
        if r > 1
            fprintf(f, '\n');
        end
    end
    fprintf(f, '}\n');
end

function [] = printmatrix(f,val)
    % again, assume only 2d
    if isa(val(1,1), 'char')
        % assume string
        fprintf(f, '"%s"', val);
    else
        fprintf(f, '[');
        [r,~] = size(val);
        for j=1:r
            row = val(j,:);
            for k=1:numel(row)
                printval(f, val(j,k));
            end
            if r > 1
                fprintf(f, '\n');
            end
        end
        fprintf(f, ']\n');
    end
end

function [] = printstruct(f, val)
    fprintf(f, 'structure(');
    fields = fieldnames(val);
    for j=1:numel(fields) 
        fprintf(f, '"%s",', fields{j});
        printval(f, val.(fields{j}));
        fprintf(f, ',');
    end
    fprintf(f, ')\n');
end

function [] = printmap(f, val)
    fprintf(f, '{');
    keyset = keys(val);
    for i=1:numel(keyset)
        key = keyset{i};
        fprintf(f, '"%s": ', key);
        printval(f, val(key));
        fprintf(f, '\n');
    end
    fprintf(f, '}');
end