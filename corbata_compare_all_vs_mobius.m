function corbata_compare_all_vs_mobius()
% Compare all Corbata numbers vs squarefree numbers (Mobius mu(n) ≠ 0).
% Figure 1: all numbers in grey with squarefree values overlaid in red.
% Figure 2: only squarefree values (Mobius mu = ±1).
% Additionally, the code prints:
%   - the count of squarefree numbers per row
%   - the density per row
%   - comparison with the theoretical densities 4/pi^2 and 8/pi^2.

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

    % Storage for all Corbata rows
    corbata_cells = cell(max_r, 1);
    corbata_cells{1} = 1;  % Row 1 contains only the value 1

    % Track used numbers to ensure no repeats
    used_numbers = containers.Map('KeyType','double','ValueType','logical');
    used_numbers(1) = true;

    % Storage for counts and lengths per row (r = 2..max_r)
    sqfree_counts  = zeros(max_r,1);
    row_lengths    = zeros(max_r,1);
    sqfree_density = zeros(max_r,1);

    % ----- 2. Build all Corbata rows -----
    for r = 2:max_r
        row_vec = generate_corbata_row_local(r);

        % Parity check: row cannot mix even and odd values
        if any(mod(row_vec,2) ~= mod(row_vec(1),2))
            error('CNS error: row %d mixes even and odd numbers.', r);
        end

        % Repetition check
        for v = row_vec
            if isKey(used_numbers, v)
                error('CNS error: number %d is repeated in row %d.', v, r);
            end
            used_numbers(v) = true;
        end

        corbata_cells{r} = row_vec;
        row_lengths(r)   = numel(row_vec);
    end

    % ----- 3. Figure 1: all numbers vs Mobius (squarefree) -----
    figure;
    hold on;

    for r = 2:max_r
        row_vec = corbata_cells{r};

        % Same x coordinate for all row entries
        x_all = r * ones(size(row_vec));

        % Plot all values in grey
        plot(x_all, row_vec, '.', 'Color', [0.7 0.7 0.7]);

        % Compute Möbius mu(n) for the row
        mu_vals = mobius_mu(row_vec);

        % Squarefree numbers correspond to mu = ±1 (mu ≠ 0)
        mask_sqfree = (mu_vals ~= 0);

        % Store counts and density
        sqfree_counts(r)  = sum(mask_sqfree);
        if row_lengths(r) > 0
            sqfree_density(r) = sqfree_counts(r) / row_lengths(r);
        end

        if any(mask_sqfree)
            vals_sqfree = row_vec(mask_sqfree);
            x_sqfree    = r * ones(size(vals_sqfree));

            % Plot squarefree numbers in red
            plot(x_sqfree, vals_sqfree, '.', 'Color', [1 0 0], 'MarkerSize', 6);
        end
    end

    xlabel('Row index r');
    ylabel('Number value');
    title('All Corbata numbers (grey) vs squarefree numbers (Mobius μ ≠ 0) (red)');
    grid on;
    hold off;

    % ----- 4. Figure 2: only squarefree values -----
    all_rows_sqfree = [];
    all_vals_sqfree = [];

    for r = 2:max_r
        vec = corbata_cells{r};
        mu_vals = mobius_mu(vec);
        mask_sqfree = (mu_vals ~= 0);

        if ~any(mask_sqfree)
            continue;
        end

        vals_row = vec(mask_sqfree);
        all_rows_sqfree = [all_rows_sqfree; r * ones(numel(vals_row),1)];
        all_vals_sqfree = [all_vals_sqfree; vals_row(:)];
    end

    if ~isempty(all_rows_sqfree)
        figure;
        scatter(all_rows_sqfree, all_vals_sqfree, 6, 'filled');
        xlabel('Row index r');
        ylabel('Squarefree value');
        title('Squarefree numbers in Corbata rows (Mobius μ ≠ 0)');
        grid on;
    else
        fprintf('No squarefree numbers found in the selected range of rows.\n');
    end

    % ----- 5. Print density comparison per row -----
    dens_even_theory = 4 / (pi^2);  % theoretical density among even integers
    dens_odd_theory  = 8 / (pi^2);  % theoretical density among odd integers

    fprintf('\nSquarefree density per row (Mobius μ ≠ 0):\n');
    fprintf('Row   length   count   density   theory    difference\n');

    for r = 2:max_r
        len_r   = row_lengths(r);
        cnt_r   = sqfree_counts(r);
        dens_r  = sqfree_density(r);

        if mod(r,2) == 0
            theory_r = dens_even_theory;
        else
            theory_r = dens_odd_theory;
        end

        diff_r = dens_r - theory_r;

        fprintf(' %3d   %5d   %5d   %7.4f   %7.4f   %+8.4f\n', ...
                r, len_r, cnt_r, dens_r, theory_r, diff_r);
    end

    % ----- 6. Global averages for even and odd rows -----
    rows = (2:max_r).';

    even_mask = (mod(rows,2) == 0);
    odd_mask  = (mod(rows,2) == 1);

    total_sqfree_even = sum(sqfree_counts(rows(even_mask)));
    total_len_even    = sum(row_lengths(rows(even_mask)));
    avg_even_density  = total_sqfree_even / total_len_even;

    total_sqfree_odd = sum(sqfree_counts(rows(odd_mask)));
    total_len_odd    = sum(row_lengths(rows(odd_mask)));
    avg_odd_density  = total_sqfree_odd / total_len_odd;

    fprintf('\nGlobal density summary:\n');
    fprintf('Even rows: observed %.6f, theory 4/pi^2 = %.6f, diff = %+8.6f\n', ...
            avg_even_density, dens_even_theory, avg_even_density - dens_even_theory);
    fprintf('Odd rows : observed %.6f, theory 8/pi^2 = %.6f, diff = %+8.6f\n', ...
            avg_odd_density,  dens_odd_theory,  avg_odd_density  - dens_odd_theory);

end



% ===================================================================
% Local row generator for the Corbata Numerical System
% ===================================================================
function row_vec = generate_corbata_row_local(r)
% Produces row r of the Corbata pattern.
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
        % EVEN ROW r = 2k
        k = r / 2;
        negcount = k - 1;
        idx = 1;

        % Left side negative jumps
        for p = negcount-1:-1:0
            idx = idx + 1;
            diff_val = -2 * (k + p - 1);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        % Right side positive jumps
        for t = 1:negcount
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

    else
        % ODD ROW r = 2k + 1
        k = (r - 1) / 2;
        idx = 1;

        % Negative jumps before the center
        for s = k-2:-1:0
            idx = idx + 1;
            diff_val = -2 * (k + s);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        % Central +2 jump
        idx = idx + 1;
        row_vec(idx) = row_vec(idx - 1) + 2;

        % Positive jumps after the center
        for t = 2:k
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end
    end
end



% ===================================================================
% Möbius function μ(n)
% Returns:
%   μ(n) = 1   if n is squarefree with an even number of prime factors
%   μ(n) = -1  if n is squarefree with an odd number of prime factors
%   μ(n) = 0   if n has any repeated prime factor (not squarefree)
% ===================================================================
function mu_vals = mobius_mu(n_vec)
    n_vec = double(n_vec(:));         % ensure column vector
    mu_vals = zeros(size(n_vec));     % initialize output

    for i = 1:numel(n_vec)
        n = n_vec(i);

        if n == 1
            mu_vals(i) = 1;
            continue;
        end

        f = factor(n);                % prime factorization

        % If any prime factor is repeated → μ(n) = 0
        if numel(f) ~= numel(unique(f))
            mu_vals(i) = 0;
        else
            % Squarefree case: μ(n) = (-1)^(number of distinct primes)
            k = numel(unique(f));
            if mod(k,2) == 0
                mu_vals(i) = 1;
            else
                mu_vals(i) = -1;
            end
        end
    end

    % return the same shape as input
    mu_vals = reshape(mu_vals, size(n_vec));
end
