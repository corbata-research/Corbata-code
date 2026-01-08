function corbata_compare_all_vs_primes()
% Compare all Corbata numbers vs only prime numbers in a plot.
% Figure 1: all numbers (grey) with primes overlaid (red).
% Figure 2: only primes (row index vs prime value).

    % ----- 1. User input -----
    while true
        prompt = ['How many Corbata rows (starting from row 2) ', ...
                  'do you want to generate and plot? '];
        num_rows = input(prompt);

        if isscalar(num_rows) && isnumeric(num_rows) && ...
                num_rows == floor(num_rows) && num_rows > 0
            break;
        else
            disp('Invalid input. Please enter a positive integer.');
        end
    end

    max_r = num_rows + 1;  % include top row [1]

    % Storage for the rows
    corbata_cells = cell(max_r, 1);
    corbata_cells{1} = 1;  % Row 1

    % Track used numbers to detect repeats (use double keys)
    used_numbers = containers.Map('KeyType','double','ValueType','logical');
    used_numbers(1) = true;

    % ----- 2. Build all Corbata rows -----
    for r = 2:max_r
        row_vec = generate_corbata_row_local(r);

        % Parity check: row cannot mix even and odd
        if any(mod(row_vec,2) ~= mod(row_vec(1),2))
            error('CNS error: row %d mixes even and odd numbers.', r);
        end

        % No repetition check
        for v = row_vec
            if isKey(used_numbers, v)
                error('CNS error: number %d is repeated in row %d.', v, r);
            end
            used_numbers(v) = true;
        end

        corbata_cells{r} = row_vec;
    end

    % ----- 3. Plot: all numbers and primes overlaid -----

    % Figure 1: all numbers (grey) and primes (red) in the same plot
    figure;
    hold on;

    % Loop again to plot row by row (no huge memory arrays)
    for r = 2:max_r
        row_vec = corbata_cells{r};

        % x-coordinates for this row (same r for all entries)
        x_all = r * ones(size(row_vec));

        % Plot ALL numbers in grey
        plot(x_all, row_vec, '.', ...
             'Color', [0.7 0.7 0.7], ...
             'MarkerSize', 12);

        % Now get primes only for this row
        prime_mask = isprime(row_vec);
        if any(prime_mask)
            primes_row = row_vec(prime_mask);
            x_primes   = r * ones(size(primes_row));

            % Plot primes in red over the grey cloud
            plot(x_primes, primes_row, '.', ...
                 'Color', [1 0 0], ...
                 'MarkerSize', 12);
        end
    end

    xlabel('Row index r');
    ylabel('Number value');
    title('All Corbata numbers (grey) vs primes (red)');
    grid on;
    hold off;

    % ----- 4. Second figure: only primes (clean view) -----
    all_rows_prime = [];
    all_vals_prime = [];

    for r = 2:max_r
        vec = corbata_cells{r};
        primes_in_row = vec(isprime(vec));
        if isempty(primes_in_row)
            continue;
        end
        all_rows_prime = [all_rows_prime; r * ones(numel(primes_in_row),1)];
        all_vals_prime = [all_vals_prime; primes_in_row(:)];
    end

    if ~isempty(all_rows_prime)
        figure;
        scatter(all_rows_prime, all_vals_prime, 12, 'filled');
        xlabel('Row index r');
        ylabel('Prime value');
        title('Prime values in Corbata rows (only primes)');
        grid on;
    else
        fprintf('No primes found in the selected range of rows.\n');
    end
end


% ================================================================
% Local row generator: Corbata rule
% ================================================================
function row_vec = generate_corbata_row_local(r)
% Produces row r of the Corbata (r >= 2).
% Row length: r-1
% Leftmost value: L = (r-2)^2 + 2

    r = double(r);

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

    if mod(r,2) == 0
        % EVEN ROW r = 2k
        k = r / 2;
        negcount = k - 1;
        idx = 1;

        % Left-side negative jumps
        for p = negcount-1:-1:0
            idx = idx + 1;
            diff_val = -2 * (k + p - 1);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        % Right-side positive jumps
        for t = 1:negcount
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

    else
        % ODD ROW r = 2k + 1
        k = (r - 1) / 2;
        idx = 1;

        % Negative jumps before center
        for s = k-2:-1:0
            idx = idx + 1;
            diff_val = -2 * (k + s);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        % Central +2 jump
        idx = idx + 1;
        row_vec(idx) = row_vec(idx - 1) + 2;

        % Positive jumps after center
        for t = 2:k
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end
    end
end
