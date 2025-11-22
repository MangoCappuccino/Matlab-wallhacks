function C = cofactorMatrix(A)
% cofactorMatrix(A)
% Returns the cofactor matrix of a square matrix A
%
% Works for both numeric and symbolic matrices.
%
% Example:
% syms a b c d
% A = [a b; c d];
% C = cofactorMatrix(A)

    % Check that A is square
    [m, n] = size(A);
    if m ~= n
        error('Matrix must be square.');
    end

    % Initialize cofactor matrix
    C = sym(zeros(m));  % Use sym to handle symbolic variables automatically

    % Loop over each element
    for i = 1:m
        for j = 1:n
            % Minor matrix (delete i-th row and j-th column)
            M = A;
            M(i, :) = [];
            M(:, j) = [];

            % Compute cofactor: (-1)^(i+j) * det(M)
            C(i,j) = (-1)^(i+j) * det(M);
        end
    end
end