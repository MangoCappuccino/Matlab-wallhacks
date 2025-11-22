function perm_cases(matrix, varargin)
    % varargin = { {'a', [0]}, {'b', [1 0 -1 2]} }
    
    varNames = cellfun(@(x) x{1}, varargin, 'UniformOutput', false);   % {'a','b'}
    varValues = cellfun(@(x) x{2}, varargin, 'UniformOutput', false);  % {[0],[1 0 -1 2]}
    
    symsList = sym(varNames); % create symbolic variables
    
    % For each variable, include symbolic version plus numeric values
    combinedValues = cell(size(varValues));
    for i = 1:numel(varValues)
        combinedValues{i} = [{symsList(i)}, num2cell(varValues{i})];
    end
    
    % Generate all combinations (symbolic + numeric)
    [C{1:numel(combinedValues)}] = ndgrid(combinedValues{:});
    combinationsList = cellfun(@(x) x(:), C, 'UniformOutput', false);
    combinationsList = [combinationsList{:}];
    
    % Loop over all combinations
    for i = 1:size(combinationsList,1)
        caseMatrix = matrix;
        caseParts = {};
        for j = 1:numel(symsList)
            val = combinationsList{i,j};
            if isa(val,'sym') % symbolic, keep as variable
                caseParts{end+1} = char(val);
            else % numeric
                caseMatrix = subs(caseMatrix, symsList(j), val);
                caseParts{end+1} = sprintf('%s = %g', varNames{j}, val);
            end
        end
        caseString = strjoin(caseParts, ', ');
        resultMatrix = rref(caseMatrix);
        fprintf('Case: %s\n', caseString);
        disp(resultMatrix);
    end
end
