function [isStoch, stochType, stationary, info] = is_stochastic_and_stationary(A, tol)
% IS_STOCHASTIC_AND_STATIONARY  Check whether A is stochastic and return equilibrium.
%
% Usage:
%   [isStoch, stochType, stationary, info] = is_stochastic_and_stationary(A)
%   [isStoch, stochType, stationary, info] = is_stochastic_and_stationary(A, tol)
%
% Inputs:
%   A   - square matrix
%   tol - (optional) numerical tolerance (default 1e-8)
%
% Outputs:
%   isStoch    - logical: true if A is (row- or column-) stochastic within tol
%   stochType  - 'row' or 'column' (which kind of stochastic), [] if not stochastic
%   stationary - equilibrium vector (column vector summing to 1) if isStoch true, []
%   info       - diagnostic string (warnings, why failed, or empty)
%
% Method:
%   - If rows sum to 1 (within tol) -> treat as row-stochastic (Markov rows).
%     stationary satisfies pi' = pi' * A -> pi is left eigenvector of A with eigenvalue 1.
%     We compute right-eigenvector of A' for eigenvalue 1, normalize, fix small noise.
%   - Else if columns sum to 1 -> treat as column-stochastic and use right eigenvector of A.
%   - Uses eigen decomposition first; if numerics are bad or multiplicity issues arise,
%     falls back to a power iteration (left or right) to find the steady state.
%
% Robustness:
%   - tolerates tiny negative noise (clamps tiny negative to 0)
%   - converts tiny complex rounding to real
%   - returns messages in 'info'

if nargin < 2 || isempty(tol)
    tol = 1e-8;
end

info = '';
stationary = [];
stochType = [];
isStoch = false;

% Basic checks
if ~ismatrix(A) || size(A,1) ~= size(A,2)
    info = 'A must be a square matrix.';
    return
end
n = size(A,1);

% Allow tiny negative rounding; flag larger negatives as invalid
minEntry = min(A(:));
if minEntry < -1e-12
    info = 'Matrix contains significantly negative entries (not stochastic).';
    return
end

% compute sums
rowSums = sum(A,2);
colSums = sum(A,1)';

isRow = all(abs(rowSums - 1) <= tol);
isCol = all(abs(colSums - 1) <= tol);

if ~isRow && ~isCol
    info = 'Rows or columns do not sum to 1 within tolerance.';
    return
end

% Decide type (prefer row if both within tol)
if isRow
    stochType = 'row';
    % left eigenvector pi' satisfies pi' = pi'*A
    % Equivalent: A'*v = 1*v where v = pi (right eigenvector of A')
    try
        [V,D] = eig(A.');
        d = diag(D);
        % find index closest to 1
        [~, idx] = min(abs(d - 1));
        v = V(:, idx);
    catch
        v = [];
    end
    % If eigen returned, try to clean; otherwise fallback
    stationary = process_eigvector(v);
    if isempty(stationary)
        % fallback: power iteration on left: iterate pi = pi * A
        stationary = power_method_left(A, n, tol);
        info = strcat(info, ' Fell back to power iteration for left eigenvector.');
    end
else
    stochType = 'column';
    try
        [V,D] = eig(A);
        d = diag(D);
        [~, idx] = min(abs(d - 1));
        v = V(:, idx);
    catch
        v = [];
    end
    stationary = process_eigvector(v);
    if isempty(stationary)
        stationary = power_method_right(A, n, tol);
        info = strcat(info, ' Fell back to power iteration for right eigenvector.');
    end
end

% Final cleanup: remove tiny negs, re-normalize
if ~isempty(stationary)
    stationary = real(stationary);
    stationary(stationary < 0 & stationary > -1e-10) = 0; % clamp tiny negative noise
    if all(stationary == 0)
        info = strcat(info, ' Stationary vector collapsed to zero.');
        isStoch = false;
        stationary = [];
        stochType = [];
        return
    end
    stationary = stationary / sum(stationary);
    isStoch = true;
end

end

%% helper: process eigenvector v -> normalized stationary (or [] if unacceptable)
function s = process_eigvector(v)
s = [];
if isempty(v), return; end
% discard if zeros or huge imaginary part
if all(abs(v) < 1e-12)
    return
end
% if small imaginary, take real part
if max(abs(imag(v))) < 1e-9
    v = real(v);
else
    % if significant complex, fail
    return
end
% if mostly negative entries then try absolute (numeric issues)
if any(v < -1e-9)
    % if negatives are tiny, take abs; otherwise fail
    if max(abs(v(v<0))) > 1e-6
        return
    else
        v = abs(v);
    end
end
% normalize
if sum(v) == 0
    return
end
s = v / sum(v);
end

%% helper: left power method (iterate row vector pi = pi * A)
function pi = power_method_left(A, n, tol)
pi = ones(1,n)/n;
maxit = 20000;
for k=1:maxit
    new = pi * A;
    if norm(new-pi,1) < tol
        pi = new';
        return
    end
    pi = new;
end
pi = pi';
end

%% helper: right power method (iterate column vector x_{k+1} = A * x_k)
function x = power_method_right(A, n, tol)
x = ones(n,1)/n;
maxit = 20000;
for k=1:maxit
    new = A * x;
    if norm(new-x,1) < tol
        x = new;
        return
    end
    x = new;
end
end