function corbata_prime_equilibrium_tool()
% Corbata prime equilibrium tool
% Generates Corbata rows and tests prime density equilibrium using
% horizontal normalized coordinate bands x in [-1, 1].
%
% What it does
% 1) Generates Corbata rows using your deterministic construction
% 2) Computes prime density profiles by x normalized bands for two windows of rows
% 3) Computes distances between the two profiles
% 4) Optional sweep where the second window moves upward to test convergence

    clc;
    fprintf('\nCorbata prime equilibrium tool\n');

    % User inputs
    num_rows = ask_positive_int('How many Corbata rows do you want to generate (starting at row 2)? ');
    num_bins = ask_positive_int('How many x bands (bins) do you want (example 40 or 60)? ');

    fprintf('\nChoose two row windows to compare.\n');
    fprintf('Important: windows are in real row index r (same r used in your plots).\n');
    fprintf('Your generated rows cover r = 2 up to r = %d.\n', num_rows + 1);

    wA1 = ask_positive_int('Window A start row rA_start: ');
    wA2 = ask_positive_int('Window A end row rA_end: ');
    wB1 = ask_positive_int('Window B start row rB_start: ');
    wB2 = ask_positive_int('Window B end row rB_end: ');

    do_sweep = ask_yes_no('Do you want a convergence sweep for Window B (yes or no)? ');

    % Generate rows
    fprintf('\nGenerating Corbata rows...\n');
    corbata_rows = generate_corbata_rows(num_rows);

    % Basic checks on window bounds
    max_real_r = num_rows + 1;
    check_window([wA1 wA2], max_real_r);
    check_window([wB1 wB2], max_real_r);

    % Compute profiles for the two windows
    fprintf('\nComputing prime density profiles...\n');

    edges = linspace(-1, 1, num_bins + 1);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));

    profA = profile_for_window(corbata_rows, wA1, wA2, edges);
    profB = profile_for_window(corbata_rows, wB1, wB2, edges);

    diffVec = profA.density - profB.density;

    L1 = mean(abs(diffVec), 'omitnan');
    L2 = sqrt(mean(diffVec.^2, 'omitnan'));

    fprintf('\nResults for Window A [%d, %d] vs Window B [%d, %d]\n', wA1, wA2, wB1, wB2);
    fprintf('Mean absolute difference L1 = %.6g\n', L1);
    fprintf('Root mean square difference L2 = %.6g\n', L2);

    % Plots
    figure;
    plot(centers, profA.density, 'LineWidth', 1.5);
    hold on;
    plot(centers, profB.density, 'LineWidth', 1.5);
    grid on;
    xlabel('x normalized position in row');
    ylabel('prime density per band');
    title(sprintf('Prime density profiles, L1 %.4g, L2 %.4g', L1, L2));
    legend(sprintf('Window A r %d to %d', wA1, wA2), sprintf('Window B r %d to %d', wB1, wB2));

    figure;
    plot(centers, diffVec, 'LineWidth', 1.5);
    grid on;
    xlabel('x normalized position in row');
    ylabel('density difference A minus B');
    title('Difference of prime density profiles');

    % Optional sweep to test convergence
    if do_sweep
        fprintf('\nConvergence sweep setup.\n');
        fprintf('We keep Window A fixed and slide Window B upward.\n');
        b_len = wB2 - wB1;
        fprintf('Current Window B length is %d rows.\n', b_len);

        sweep_start = ask_positive_int('Sweep start row for Window B: ');
        sweep_end   = ask_positive_int('Sweep end row for Window B: ');
        sweep_step  = ask_positive_int('Sweep step in rows: ');

        check_window([sweep_start sweep_start + b_len], max_real_r);
        check_window([sweep_end   sweep_end   + b_len], max_real_r);

        centersB = [];
        L1s = [];
        L2s = [];

        fprintf('\nRunning sweep...\n');
        for b1 = sweep_start:sweep_step:sweep_end
            b2 = b1 + b_len;
            if b2 > max_real_r
                break;
            end

            profB_s = profile_for_window(corbata_rows, b1, b2, edges);
            d = profA.density - profB_s.density;

            centersB(end+1) = 0.5 * (b1 + b2);
            L1s(end+1) = mean(abs(d), 'omitnan');
            L2s(end+1) = sqrt(mean(d.^2, 'omitnan'));

            fprintf('Window B r %d to %d, L1 %.6g, L2 %.6g\n', b1, b2, L1s(end), L2s(end));
        end

        figure;
        plot(centersB, L1s, 'LineWidth', 1.5);
        grid on;
        xlabel('center of Window B');
        ylabel('L1 mean absolute difference');
        title('Convergence test, does L1 decrease as Window B moves upward');

        figure;
        plot(centersB, L2s, 'LineWidth', 1.5);
        grid on;
        xlabel('center of Window B');
        ylabel('L2 root mean square difference');
        title('Convergence test, does L2 decrease as Window B moves upward');
    end

    fprintf('\nDone.\n');

end

function rows = generate_corbata_rows(num_rows)
% Returns a cell array where rows{r} is the row vector for real row index r
% The cell array includes rows from r = 1 to r = num_rows + 1
% Row 1 is [1]
% Row 2.. are generated by the Corbata rule

    max_real_r = num_rows + 1;

    rows = cell(max_real_r, 1);
    rows{1} = 1;

    used_numbers = containers.Map('KeyType','int32','ValueType','logical');
    used_numbers(int32(1)) = true;

    for r = 2:max_real_r
        row_vec = generate_corbata_row(r);

        row_parity = mod(row_vec(1), 2);
        if any(mod(row_vec, 2) ~= row_parity)
            error('CNS error: row %d mixes even and odd numbers.', r);
        end

        for v = row_vec
            vv = int32(v);
            if isKey(used_numbers, vv) && used_numbers(vv)
                error('CNS error: number %d is repeated in row %d.', v, r);
            else
                used_numbers(vv) = true;
            end
        end

        rows{r} = row_vec;
    end
end

function prof = profile_for_window(rows, rStart, rEnd, edges)
% Computes prime density profile by x normalized bands for rows rStart..rEnd
% rows is indexed by real r, so rows{r} is row r
% x normalized is in [-1, 1] within each row

    numBins = numel(edges) - 1;
    primeCounts = zeros(1, numBins);
    totalCounts = zeros(1, numBins);

    for r = rStart:rEnd
        rowVec = rows{r};
        n = numel(rowVec);

        center = (n + 1) / 2;
        halfwidth = (n - 1) / 2;

        if halfwidth == 0
            x = 0;
        else
            x = ((1:n) - center) / halfwidth;
        end

        binIdx = discretize(x, edges);
        isP = isprime(rowVec);

        for b = 1:numBins
            mask = (binIdx == b);
            totalCounts(b) = totalCounts(b) + sum(mask);
            primeCounts(b) = primeCounts(b) + sum(isP(mask));
        end
    end

    density = primeCounts ./ totalCounts;

    prof.primeCounts = primeCounts;
    prof.totalCounts = totalCounts;
    prof.density = density;
end

function row_vec = generate_corbata_row(r)
% Your row generator logic, unchanged

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
        k = r/2;
        negcount = k - 1;
        idx = 1;

        for p = negcount - 1 : -1 : 0
            idx = idx + 1;
            diff_val = -2 * (k + p - 1);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        for t = 1 : negcount
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

    else
        k = (r - 1)/2;
        idx = 1;

        for s = k - 2 : -1 : 0
            idx = idx + 1;
            diff_val = -2 * (k + s);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end

        idx = idx + 1;
        row_vec(idx) = row_vec(idx - 1) + 2;

        for t = 2 : k
            idx = idx + 1;
            diff_val = 2 * (k + t);
            row_vec(idx) = row_vec(idx - 1) + diff_val;
        end
    end
end

function x = ask_positive_int(prompt)
    while true
        x = input(prompt);
        if isscalar(x) && isnumeric(x) && x == floor(x) && x > 0
            return;
        end
        fprintf('Invalid input. Enter a positive integer.\n');
    end
end

function ok = ask_yes_no(prompt)
    while true
        s = input(prompt, 's');
        s = lower(strtrim(s));
        if strcmp(s, 'yes') || strcmp(s, 'y')
            ok = true;
            return;
        end
        if strcmp(s, 'no') || strcmp(s, 'n')
            ok = false;
            return;
        end
        fprintf('Please type yes or no.\n');
    end
end

function check_window(win, max_real_r)
    if win(1) < 1 || win(2) < 1 || win(2) < win(1)
        error('Invalid window. Check start and end.');
    end
    if win(2) > max_real_r
        error('Window end exceeds generated rows. Generate more rows or reduce the window.');
    end
end
