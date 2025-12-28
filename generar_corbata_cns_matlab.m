function corbata_cns = generate_corbata_cns()
% Generates rows of the Corbata Numeric System (CNS)
% Row 1 is [1]. Later rows follow the Corbata diamond pattern.

    % ----- 1. User input -----
    while true
        prompt = ['How many Corbata rows (starting from row 2) ' ...
                  'do you want to generate? '];
        num_rows = input(prompt);

        % check that input is a positive integer
        if isscalar(num_rows) && isnumeric(num_rows) && ...
                num_rows == floor(num_rows) && num_rows > 0
            break;
        else
            disp('Invalid input. Please enter a positive integer.');
        end
    end

    % total rows including the top row [1]
    max_row_index = num_rows + 1;

    % ----- 2. Storage -----
    corbata_cells = cell(max_row_index, 1);
    corbata_cells{1} = 1;            % Row 1: [1]

    % map used to detect repeated numbers
    used_numbers = containers.Map('KeyType','int32','ValueType','logical');
    used_numbers(1) = true;          % 1 already used

    % ----- 3. Build rows r = 2 ... max_row_index -----
    for r = 2:max_row_index
        row_vec = generate_corbata_row(r);

        % 3.1 parity check
        row_parity = mod(row_vec(1), 2);
        if any(mod(row_vec, 2) ~= row_parity)
            error('CNS error: row %d mixes even and odd numbers.', r);
        end

        % 3.2 repetition check
        for v = row_vec
            if isKey(used_numbers, v) && used_numbers(v)
                error('CNS error: number %d is repeated in row %d.', v, r);
            else
                used_numbers(v) = true;
            end
        end

        % 3.3 store row
        corbata_cells{r} = row_vec;
    end

    % ----- 4. Output on screen -----
    corbata_cns = corbata_cells(2:end);

    fprintf('\n--- Corbata Numeric System rows ---\n');
    for r = 2:max_row_index
        vec = corbata_cells{r};
        fprintf('Row %d (length %d): [', r, numel(vec));
        if numel(vec) > 1
            fprintf('%d, ', vec(1:end-1));
        end
        fprintf('%d]\n', vec(end));
    end

    % ----- 5. SAVE as columns in a CSV -----

    % Find the maximum row length to pad the columns
    max_len = max(cellfun(@numel, corbata_cells(2:end)));

    % Create a rectangular cell array for the CSV
    csv_matrix = cell(max_len, num_rows);

    for col = 1:num_rows
        row_vec = corbata_cells{col + 1};  % +1 to skip row 1
        L = numel(row_vec);

        % fill the column downward
        for r = 1:L
            csv_matrix{r, col} = row_vec(r);
        end

        % pad the remaining cells with empty strings
        for r = L+1:max_len
            csv_matrix{r, col} = '';
        end
    end

    % write the CSV
    csvname = 'corbata_output_columns.csv';
    fid = fopen(csvname, 'w');

    for r = 1:max_len
        for c = 1:num_rows
            val = csv_matrix{r,c};
            if c < num_rows
                if isempty(val)
                    fprintf(fid, ',');
                else
                    fprintf(fid, '%d,', val);
                end
            else
                if isempty(val)
                    fprintf(fid, '\n');
                else
                    fprintf(fid, '%d\n', val);
                end
            end
        end
    end

    fclose(fid);
    fprintf('Corbata saved as column-based CSV to "%s".\n', csvname);

end

% =========================
% Row generator
% =========================

function row_vec = generate_corbata_row(r)

    if r == 1
        row_vec = 1;
        return;
    end

    L = (r - 2)^2 + 2;

    if r == 2
        row_vec = L;
        return;
    end

    row_len = r - 1;
    row_vec = zeros(1, row_len);
    row_vec(1) = L;

    if mod(r, 2) == 0
        % even row r = 2k
        k = r/2;
        negcount = k - 1;
        idx = 1;

        % negative differences
        for p = negcount - 1 : -1 : 0
            idx = idx + 1;
            diff_val = -2 * (k + p - 1);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        % positive differences
        for t = 1 : negcount
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

    else
        % odd row r = 2k + 1
        k = (r - 1)/2;
        idx = 1;

        % negative differences
        for s = k - 2 : -1 : 0
            idx = idx + 1;
            diff_val = -2 * (k + s);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        % central +2
        idx = idx + 1;
        row_vec(idx) = row_vec(idx - 1) + 2;

        % positive differences
        for t = 2 : k
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end
    end
end
