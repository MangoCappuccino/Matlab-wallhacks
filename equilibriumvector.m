function pi = equilibriumvector(A)
% EQUILIBRIUMVECTOR  Return the stationary distribution of a stochastic matrix.
%
%   pi = equilibriumvector(A)
%
%   If A is row-stochastic (rows sum to 1), pi satisfies  pi' = pi' * A.
%   If A is column-stochastic (columns sum to 1), pi satisfies  A * pi = pi.
%
%   Output:
%       pi  - stationary vector (column, sums to 1)
%            returns [] if A is not stochastic

tol = 1e-10;

% Check if square
[n,m] = size(A);
if n ~= m
    pi = [];
    return
end

% Determine type of stochasticity
rowSums = sum(A,2);
colSums = sum(A,1)';

isRow = all(abs(rowSums - 1) < tol);
isCol = all(abs(colSums - 1) < tol);

if ~isRow && ~isCol
    pi = [];
    return
end

% Compute the stationary vector
if isRow
    % Row-stochastic: find left eigenvector of A with eigenvalue 1
    [V,D] = eig(A');      % A' * v = v
    d = diag(D);
    [~, idx] = min(abs(d - 1));
    v = V(:, idx);
else
    % Column-stochastic: find right eigenvector of A with eigenvalue 1
    [V,D] = eig(A);       % A * v = v
    d = diag(D);
    [~, idx] = min(abs(d - 1));
    v = V(:, idx);
end

% Clean and normalize
v = real(v);              % remove tiny imaginary noise
v = max(v,0);             % clamp negatives
pi = v / sum(v);          % normalize to sum to 1

end
