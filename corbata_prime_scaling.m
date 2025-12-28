function corbata_prime_scaling()

    % -------------------------------
    % 1. Ask user for number of rows
    % -------------------------------
    num_rows = input('How many Corbata rows (starting from row 2) do you want to generate? ');

    max_r = num_rows + 1;  % because row 1 = [1]

    % Storage
    prime_count = zeros(max_r, 1);
    valid_r     = [];      % rows that contain at least one prime

    fprintf('\nGenerating rows and counting primes...\n');

    % -------------------------------
    % 2. Generate Corbata rows
    % -------------------------------
    for r = 2:max_r
        row_vec = generate_corbata_row(r);

        % extract primes
        primes_in_row = row_vec(isprime(row_vec));
        prime_count(r) = numel(primes_in_row);

        if ~isempty(primes_in_row)
            valid_r(end+1) = r;
        end
    end

    % only consider rows that actually had primes
    valid_prime_counts = prime_count(valid_r);

    % -------------------------------
    % 3. Compute scaling factor c(r)
    % -------------------------------
    c_vals = zeros(size(valid_r));

    for k = 1:length(valid_r)
        r = valid_r(k);
        expected_log = (r - 1) / log(r);
        c_vals(k) = valid_prime_counts(k) / expected_log;
    end

    % -------------------------------
    % 4. Plot c(r)
    % -------------------------------
    figure;
    plot(valid_r, c_vals, '.', 'MarkerSize', 8, 'Color', [0 0.5 1]); hold on;
    yline(0.5, 'r', 'LineWidth', 1.5);
    xlabel('Row index r');
    ylabel('c(r) = observed primes / ((r-1)/log(r))');
    title('Scaling factor c(r) for prime density in the Corbata');

    legend('Observed c(r)', 'c = 0.5 (hypothesis)');
    grid on;

    fprintf('\nDone. Plot generated.\n');
end


% =====================================================================
% Corbata row generator (same version as before)
% =====================================================================
function row_vec = generate_corbata_row(r)

    if r == 1
        row_vec = 1;
        return;
    end

    % leftmost value:
    L = (r - 2)^2 + 2;

    if r == 2
        row_vec = L;
        return;
    end

    row_len = r - 1;
    row_vec = zeros(1, row_len);
    row_vec(1) = L;

    if mod(r, 2) == 0
        % EVEN rows
        k = r / 2;
        negcount = k - 1;
        idx = 1;

        % negative steps
        for p = negcount - 1 : -1 : 0
            idx = idx + 1;
            diff_val = -2 * (k + p - 1);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        % positive steps
        for t = 1 : negcount
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

    else
        % ODD rows
        k = (r - 1) / 2;
        idx = 1;

        % negative steps
        for s = k - 2 : -1 : 0
            idx = idx + 1;
            diff_val = -2 * (k + s);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        % central +2
        idx = idx + 1;
        row_vec(idx) = row_vec(idx - 1) + 2;

        % positive steps
        for t = 2 : k
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end
    end
end
