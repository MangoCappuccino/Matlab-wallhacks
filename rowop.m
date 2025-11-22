function A = rowop(A, i, op, k, j)
% ROWOP performs an elementary row operation on matrix A.
%
% Supported operations:
%   Row scaling:
%       A = rowop(A, i, '*', k)        % Ri = k * Ri
%
%   Row replacement:
%       A = rowop(A, i, '+', k, j)     % Ri = Ri + k * Rj
%       A = rowop(A, i, '-', k, j)     % Ri = Ri - k * Rj
%
%   Row swap:
%       A = rowop(A, i, '<->', [], j)  % Ri <-> Rj
%
% Inputs:
%   A : matrix
%   i : row index being modified (or first row for swap)
%   op: operation ('*', '+', '-', '<->')
%   k : scalar multiplier (ignored for swap)
%   j : second row index (only needed for replacement/swap)

    switch op
        case '*'
            % Row scaling
            A(i,:) = k * A(i,:);
            
        case '+'
            % Row replacement with addition
            A(i,:) = A(i,:) + k * A(j,:);
            
        case '-'
            % Row replacement with subtraction
            A(i,:) = A(i,:) - k * A(j,:);
            
        case '<->'
            % Row swap
            temp = A(i,:);
            A(i,:) = A(j,:);
            A(j,:) = temp;
            
        otherwise
            error('Invalid operator. Use "*", "+", "-", or "<->".');
    end
end