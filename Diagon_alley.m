function Diagon_alley_fixed(A, tolerance)
% DIAGON_ALLEY_FIXED  Compute eigenvalues/eigenspaces and build P_GSP/D_GSP.
%   Diagon_alley_fixed(A)
%   Diagon_alley_fixed(A, tolerance)
%
% Automatically orthonormalizes if A is symmetric (orthogonally diagonalizable).
% Works with symbolic input (auto-converts integer matrices) and numeric input.
% Returns P_GSP and D_GSP in base workspace and prints a summary.
%
% NOTE: This is a corrected version; it fixes the Gram-Schmidt projection,
% symbolic zero checks, and type handling.

    if nargin < 2 || isempty(tolerance)
        tolerance = 1e-8;
    end

    % Auto-convert integer/rational matrices to symbolic for exact results.
    use_sym = false;
    if exist('sym', 'file') && ~isa(A, 'sym') && all(abs(A(:) - round(A(:))) < 1e-10)
        A = sym(A);
        use_sym = true;
        disp('→ Auto-converted matrix to symbolic for exact results.');
    elseif isa(A, 'sym')
        use_sym = true;
    end

    n = size(A, 1);
    if size(A, 1) ~= size(A, 2)
        error('Matrix A must be square.');
    end

    fprintf('\n=== Eigenvalue Analysis ===\n');

    % Symmetry test
    if use_sym
        is_symmetric = isequal(A, A.');
    else
        is_symmetric = issymmetric(A);
    end

    % Compute eigenvalues (symbolic or numeric)
    if use_sym
        eigvals_sym = eig(A); % symbolic eigenvalues
        eigvals_all = double(vpa(eigvals_sym, 12));
    else
        eigvals_all = eig(A);
        eigvals_sym = []; % placeholder
    end

    % Filter to keep only real eigenvalues (within tol)
    eigvals_all = eigvals_all(abs(imag(eigvals_all)) < tolerance);
    eigvals_all(abs(eigvals_all) < tolerance) = 0;
    unique_vals = unique(round(eigvals_all(:) * 1e12) / 1e12);
    unique_vals = sort(unique_vals);

    total_geom_mult = 0;
    summary = struct();
    if use_sym
        P_cols = sym([]);
        D_diag = sym([]);
    else
        P_cols = [];
        D_diag = [];
    end

    for k = 1:length(unique_vals)
        lambda_val = unique_vals(k);

        % pick a symbolic eigenvalue if available that matches numerically
        if use_sym
            idx = find(abs(double(eigvals_sym - lambda_val)) < max(1e-10, 100*tolerance));
            if ~isempty(idx)
                lambda = eigvals_sym(idx(1));
            else
                lambda = sym(lambda_val);
            end
        else
            lambda = lambda_val;
        end

        fprintf('\nEigenvalue #%d: %s\n', k, char(vpa(lambda, 6)));
        alg_mult = sum(abs(eigvals_all - lambda_val) < max(1e-8, 10*tolerance));

        % ---- Eigenspace ----
        if use_sym
            E = null(A - lambda * eye(n)); % symbolic nullspace
            if ~isempty(E)
                E = simplify(E, 'Steps', 20);
                try
                    E = rationalize(E);
                catch
                    % ignore rationalize failure
                end
            end
        else
            E = null(double(A) - lambda * eye(n), tolerance);
            if isempty(E)
                % fallback: try eigenvectors from eig
                [V, D] = eig(double(A));
                idxv = find(abs(diag(D) - lambda_val) < max(1e-6, 100*tolerance));
                if ~isempty(idxv)
                    E = V(:, idxv);
                end
            end
        end

        % Filter out complex components (keep columns that are real within tol)
        if ~isempty(E)
            if use_sym
                % convert columns that are purely real
                real_cols = true(1, size(E, 2));
                for c = 1:size(E, 2)
                    col = E(:, c);
                    if ~isequal(simplify(imag(col)), sym(zeros(size(col))))
                        real_cols(c) = false;
                    end
                end
                E = E(:, real_cols);
            else
                real_cols = all(abs(imag(E)) < tolerance, 1);
                E = real(E(:, real_cols));
            end
        end

        geom_mult = size(E, 2);
        total_geom_mult = total_geom_mult + geom_mult;

        % ---- Orthonormalization ----
        if is_symmetric && ~isempty(E)
            E = gs_orthonormal(E, use_sym, tolerance);
        elseif ~is_symmetric && ~isempty(E)
            E = normalize_columns(E, use_sym);
        end

        summary(k).eigenvalue = lambda;
        summary(k).alg_mult = alg_mult;
        summary(k).geom_mult = geom_mult;
        summary(k).eigenspace = E;

        fprintf('  Algebraic multiplicity: %d\n', alg_mult);
        fprintf('  Geometric multiplicity: %d\n', geom_mult);
        fprintf('  Eigenspace basis (columns):\n');
        disp(E);

        assignin('base', sprintf('eigenspace_basis_%d', k), E);

        if geom_mult > 0
            if use_sym
                P_cols = [P_cols, E];
                D_diag = [D_diag; repmat(lambda, size(E, 2), 1)];
            else
                P_cols = [P_cols, E];
                D_diag = [D_diag; repmat(lambda, size(E, 2), 1)];
            end
        end
    end

    % ---- Summary ----
    is_diag = (total_geom_mult == n);
    assignin('base', 'eigen_summary_GSP', summary);
    assignin('base', 'is_diagonalizable_GSP', is_diag);

    fprintf('\n=== Summary ===\n');
    if is_diag
        fprintf('✓ Matrix IS diagonalizable.\n');
        P = P_cols;
        if use_sym
            D = diag(D_diag);
        else
            D = diag(D_diag);
        end
        assignin('base', 'P_GSP', P);
        assignin('base', 'D_GSP', D);
        if is_symmetric
            fprintf('→ Symmetric matrix: P_GSP is orthonormal.\n');
        else
            fprintf('→ P_GSP and D_GSP saved in original form.\n');
        end
        fprintf('   Verify: simplify(P_GSP^-1 * A * P_GSP) = D_GSP\n');
    else
        fprintf('✗ Matrix is NOT diagonalizable.\n');
    end
    fprintf('===========================\n\n');
end

%% ===== Helper: Gram-Schmidt then normalize =====
function Qn = gs_orthonormal(V, use_sym, tolerance)
    if use_sym
        m = size(V, 1);
        Q = sym(zeros(m, 0));
        for j = 1:size(V, 2)
            v = V(:, j);
            for i = 1:size(Q, 2)
                qi = Q(:, i);
                denom = simplify(qi' * qi);
                if isequal(simplify(denom), sym(0))
                    continue;
                end
                proj = simplify((qi' * v) / denom) * qi;
                v = simplify(v - proj);
            end
            if ~isequal(simplify(v' * v), sym(0))
                Q(:, end+1) = simplify(v);
            end
        end
        Qn = normalize_columns(Q, true);
    else
        m = size(V, 1);
        Q = zeros(m, 0);
        for j = 1:size(V, 2)
            v = V(:, j);
            for i = 1:size(Q, 2)
                qi = Q(:, i);
                denom = (qi' * qi);
                if abs(denom) < eps, continue; end
                proj = ((qi' * v) / denom) * qi;
                v = v - proj;
            end
            if norm(v) > tolerance
                Q(:, end+1) = v;
            end
        end
        Qn = normalize_columns(Q, false);
    end
end

%% ===== Helper: normalize columns safely =====
function Qn = normalize_columns(Q, use_sym)
    [m, p] = size(Q);
    if use_sym
        Qn = sym(zeros(m, p));
        for i = 1:p
            v = Q(:, i);
            nrm = simplify(sqrt(v' * v));
            if ~isequal(simplify(nrm), sym(0))
                Qn(:, i) = simplify(v / nrm);
            end
        end
    else
        Qn = zeros(m, p);
        for i = 1:p
            v = Q(:, i);
            nrm = norm(v);
            if nrm > eps
                Qn(:, i) = v / nrm;
            end
        end
    end
end