function corbata_area_difference()
% Computes the "area" difference in the Corbata:
% A_total  = total number of Corbata points (all numbers)
% A_primes = total number of prime points
% DeltaA   = A_total - A_primes

    % ----- 1. User input -----
    while true
        prompt = ['How many Corbata rows (starting from row 2) ', ...
                  'do you want to include? '];
        num_rows = input(prompt);

        if isscalar(num_rows) && isnumeric(num_rows) && ...
                num_rows == floor(num_rows) && num_rows > 0
            break;
        else
            disp('Invalid input. Please enter a positive integer.');
        end
    end

    max_r = num_rows + 1;    % last row index (including top row [1])

    % Storage (optional, in case we later want plots)
    total_points_cum  = zeros(num_rows,1);  % cumulative all numbers
    prime_points_cum  = zeros(num_rows,1);  % cumulative primes

    % Row cache
    corbata_cells = cell(max_r,1);
    corbata_cells{1} = 1;   % row 1 = [1]

    % Repetition check map
    used_numbers = containers.Map('KeyType','double','ValueType','logical');
    used_numbers(1) = true;

    total_points = 0;
    prime_points = 0;

    % ----- 2. Build rows and count points -----
    for r = 2:max_r
        row_vec = generate_corbata_row_local(r);

        % parity check
        if any(mod(row_vec,2) ~= mod(row_vec(1),2))
            error('CNS error: row %d mixes even and odd numbers.', r);
        end

        % repetition check
        for v = row_vec
            if isKey(used_numbers, v)
                error('CNS error: number %d is repeated in row %d.', v, r);
            end
            used_numbers(v) = true;
        end

        corbata_cells{r} = row_vec;

        % count all numbers in this row
        total_points = total_points + numel(row_vec);

        % count primes in this row
        primes_in_row = row_vec(isprime(row_vec));
        prime_points  = prime_points + numel(primes_in_row);

        idx = r - 1;  % index from 1..num_rows
        total_points_cum(idx) = total_points;
        prime_points_cum(idx) = prime_points;
    end

    % ----- 3. Area difference and fractions -----
    DeltaA = total_points - prime_points;

    R = max_r;  % last row index
    frac_primes  = prime_points / total_points;      % observed
    frac_expected = 1 / log(R);                      % ~ 1/log(R)

    fprintf('\n=== Corbata area summary (rows 2 to %d) ===\n', max_r);
    fprintf('Total points (all numbers)     A_total  = %d\n', total_points);
    fprintf('Total prime points            A_primes = %d\n', prime_points);
    fprintf('Difference (composites area)  DeltaA   = %d\n', DeltaA);
    fprintf('\nObserved prime area fraction        = %.6g\n', frac_primes);
    fprintf('Approx expected fraction 1/log(R)   = %.6g  (R = %d)\n', ...
            frac_expected, R);
end


% ================================================================
% Local row generator (same Corbata rule as antes)
% ================================================================
function row_vec = generate_corbata_row_local(r)
% Returns row r (r >= 2) of the Corbata.
% Row length: r - 1
% Leftmost value: L = (r - 2)^2 + 2

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
        % EVEN row r = 2k
        k = r / 2;
        negcount = k - 1;
        idx = 1;

        % left negative jumps
        for p = negcount-1:-1:0
            idx = idx + 1;
            diff_val = -2 * (k + p - 1);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        % right positive jumps
        for t = 1:negcount
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

    else
        % ODD row r = 2k + 1
        k = (r - 1) / 2;
        idx = 1;

        % negative jumps before center
        for s = k-2:-1:0
            idx = idx + 1;
            diff_val = -2 * (k + s);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        % central +2 jump
        idx = idx + 1;
        row_vec(idx) = row_vec(idx - 1) + 2;

        % positive jumps after center
        for t = 2:k
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end
    end
end
