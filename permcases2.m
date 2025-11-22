function perm_cases(matrix, varargin)
    % varargin = { {'a', [0]}, {'b', [1 0 -1 2]} }

    varNames = cellfun(@(x) x{1}, varargin, 'UniformOutput', false);
    varValues = cellfun(@(x) x{2}, varargin, 'UniformOutput', false);
    symsList = sym(varNames);

    % Combine symbolic + numeric values for each variable
    combinedValues = cell(size(varValues));
    for i = 1:numel(varValues)
        combinedValues{i} = [{symsList(i)}, num2cell(varValues{i})];
    end

    % Create all combinations
    [C{1:numel(combinedValues)}] = ndgrid(combinedValues{:});
    combinationsList = cellfun(@(x) x(:), C, 'UniformOutput', false);
    combinationsList = [combinationsList{:}];

    % Store classification results
    results = {};

    % Loop through each combination
    for i = 1:size(combinationsList,1)
        caseMatrix = matrix;
        caseDesc = struct();

        % Build description + substitute numeric parameters
        for j = 1:numel(symsList)
            val = combinationsList{i,j};
            if isa(val,'sym')      % symbolic
                caseDesc.(varNames{j}) = 'symbolic';
            else                   % numeric
                caseMatrix = subs(caseMatrix, symsList(j), val);
                caseDesc.(varNames{j}) = val;
            end
        end

        % Compute RREF
        r = rref(caseMatrix);

        % Detect inconsistency / uniqueness / infinite solutions
        [m,n] = size(r);
        rankA = rank(r(:,1:n-1));   % rank of coefficient matrix
        rankAug = rank(r);          % rank of augmented matrix
        numVars = n-1;

        if rankAug > rankA
            category = "INCONSISTENT";
        elseif rankA == numVars
            category = "CONSISTENT & UNIQUE";
        else
            category = "CONSISTENT & INFINITE";
        end

        % Print the case
        fprintf("\n==== Case: ");
        disp(caseDesc);
        disp(r);

        % Save for summary
        results{end+1} = struct("params", caseDesc, "type", category);
    end

    % Final summary table
    fprintf("\n================ SUMMARY ================\n");

    for k = 1:numel(results)
        p = results{k}.params;
        t = results{k}.type;

        fprintf("Case %d: ", k);
        fields = fieldnames(p);

        for f = 1:numel(fields)
            v = p.(fields{f});
            if ischar(v)
                fprintf("%s = symbolic, ", fields{f});
            else
                fprintf("%s = %g, ", fields{f}, v);
            end
        end

        fprintf(" ==> %s\n", t);
    end

    fprintf("=========================================\n");
end