function EROs = practicalEROs_check(A, B)
    % practicalEROs_check: Find EROs to transform A into B if possible.
    % If impossible, throws an error with the reason.

    tol = 1e-10; % numerical tolerance
    [m,n] = size(A);

    if ~isequal(size(B), [m n])
        error('Impossible: A and B must be the same size.');
    end

    % Rank must match
    if rank(A) ~= rank(B)
        error('Impossible: rank(A) ~= rank(B). Row ops preserve rank.');
    end

    % Zero column invariant
    for j = 1:n
        if all(abs(A(:,j)) < tol) && any(abs(B(:,j)) > tol)
            error('Impossible: Column %d is zero in A but nonzero in B.', j);
        end
    end

    % Start transformation
    C = A;
    EROs = {}; % store operations

    for col = 1:min(m,n)
        % Find pivot
        pivot_row = find(abs(C(col:end,col)) > tol,1) + col - 1;
        if isempty(pivot_row)
            continue;
        end

        % Swap rows if needed
        if pivot_row ~= col
            C([col pivot_row],:) = C([pivot_row col],:);
            EROs{end+1} = sprintf('C([%d %d],:) = C([%d %d],:); %% swap R%d <-> R%d', ...
                                   col, pivot_row, pivot_row, col, col, pivot_row);
        end

        % Scale pivot row to match B
        if abs(C(col,col)) > tol && abs(B(col,col)) > tol
            factor = B(col,col) / C(col,col);
            if abs(factor - 1) > tol
                C(col,:) = C(col,:) * factor;
                EROs{end+1} = sprintf('C(%d,:) = C(%d,:) * %.5g; %% R%d -> R%d * %.5g', ...
                                       col, col, factor, col, col, factor);
            end
        end

        % Eliminate other entries in column
        for row = 1:m
            if row ~= col && abs(C(row,col) - B(row,col)) > tol
                factor2 = (C(row,col) - B(row,col)) / C(col,col);
                C(row,:) = C(row,:) - factor2 * C(col,:);
                EROs{end+1} = sprintf('C(%d,:) = C(%d,:) - %.5g*C(%d,:); %% R%d -> R%d - %.5g*R%d', ...
                                       row, row, factor2, col, row, row, factor2, col);
            end
        end
    end

    % Final check
    if norm(C - B, 'fro') > 1e-8
        error('Impossible: after EROs, result does not match B. Likely not reachable.');
    end

    % Display results
    disp('%% Sequence of EROs to go from A to B:');
    for k = 1:length(EROs)
        disp(EROs{k});
    end
end