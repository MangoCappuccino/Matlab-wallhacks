function A = ref(A)
    % Temporarily disable warnings related to symbolic evaluation
    warning('off');

    % Ensure the matrix is symbolic if needed
    A = sym(A); 

    [m, n] = size(A); % Get the size of the matrix
    lead = 1; % Leading entry column

    for r = 1:m
        if n < lead
            return; % If there are no more columns to process, exit
        end

        % Find a row with a non-zero leading entry
        i = r;
        while isAlways(A(i, lead) == 0)
            i = i + 1;
            if i > m
                i = r;
                lead = lead + 1; % Move to the next column
                if n < lead
                    return; % If no more leads, exit
                end
            end
        end

        % Swap rows if needed
        if i ~= r
            A([i, r], :) = A([r, i], :);
        end

        % Eliminate entries below the pivot (leading entry)
        for i = r+1:m
            lv = A(i, lead) / A(r, lead); % Calculate the elimination factor
            if ~isAlways(lv == 0)
                A(i, :) = simplify(A(i, :) - lv * A(r, :)); % Row operation
            end
        end

        lead = lead + 1; % Move to the next lead
    end

    % Re-enable warnings after the function finishes
    warning('on');
end
