function [U, Sigma, V] = symbolicSVD(A_input)
    % symbolicSVD  Calculates the SVD of a matrix using the symbolic, step-by-step method.
    %
    %   [U, Sigma, V] = symbolicSVD(A)
    %
    %   This function follows the exact method from your notes:
    %   1. Finds V and Sigma from the eigendecomposition of A'*A.
    %   2. Sorts the eigenvalues in descending order.
    %   3. Finds U using the formula u_i = (1/sigma_i) * A * v_i.
    %   4. Completes U by finding an orthonormal basis for the null space
    %      of the columns of U found so far (your "Method 1").
    %
    %   Input:
    %     A_input : The matrix to decompose. It will be converted to
    %               symbolic to produce "nice" fractional answers.
    %
    %   Outputs:
    %     U     : The m x m orthogonal U matrix.
    %     Sigma : The m x n diagonal matrix of singular values (S).
    %     V     : The n x n orthogonal V matrix.
    
    % Ensure input is symbolic for exact fractional math
    A = sym(A_input);
    
    [m, n] = size(A);

    % --- Step 1: Find V and Sigma ---
    % (Based on your "Problem 1(a)" notes)
    
    % A'*A = V * (S'*S) * V'
    % The columns of V are the eigenvectors of A'*A
    ATA = A'*A;
    
    % Find eigenvectors (P) and eigenvalues (D)
    % P's columns are eigenvectors
    % D is a diagonal matrix of eigenvalues
    [P, D] = eig(ATA);
    
    % Sort eigenvalues in descending order (from your "Epilogue" hint)
    [sorted_eigenvalues, sort_indices] = sort(diag(D), 'descend');
    
    % Get singular values (sqrt of eigenvalues)
    singular_values = sqrt(sorted_eigenvalues);
    
    % Sort eigenvectors (columns of V) to match the sorted eigenvalues
    V = P(:, sort_indices);
    
    % --- FIX 1: Normalize eigenvectors in V ---
    % The 'eig' function on symbolic matrices does not guarantee normalized vectors.
    % We must normalize them to ensure V is orthogonal.
    for i = 1:n
        col_vec = V(:, i);
        if norm(col_vec) ~= 0 % Avoid division by zero
            V(:, i) = col_vec / norm(col_vec);
        end
    end
    % --- End Fix 1 ---
    
    % Build the Sigma matrix (must be size m x n)
    Sigma = sym(zeros(m, n));
    num_singular_values = length(singular_values);
    
    % Place the diagonal values
    Sigma(1:num_singular_values, 1:num_singular_values) = diag(singular_values);

    % --- Step 2: Find U ---
    % (Based on your "Problem 1(b)" notes)
    
    U_partial = sym([]); % Initialize the first part of U
    
    % Find the rank (number of non-zero singular values)
    % Use a small tolerance for any potential symbolic float artifacts
    rank_r = sum(double(singular_values) > 1e-10);
    
    % Calculate the first 'r' columns of U
    % u_i = (1/sigma_i) * A * v_i
    % Since V is now normalized, this will produce normalized u_i columns.
    for i = 1:rank_r
        sigma_i = singular_values(i);
        v_i = V(:, i);
        
        u_i = (1/sigma_i) * A * v_i;
        U_partial = [U_partial, u_i]; % Add the new vector as a column
    end
    
    % If A is "tall" (m > r), we need to find the remaining m-r columns of U
    % These columns form an orthonormal basis for the null space of U_partial'
    % (This is your "Method 1: null(U')")
    if m > rank_r
        % The 'null' function on symbolic matrices returns a basis,
        % but it is not guaranteed to be orthonormal.
        U_complement_basis = null(U_partial');
        
        % --- FIX 2: Normalize the null space basis vectors ---
        U_complement = sym([]); % Initialize empty
        num_complement_vecs = size(U_complement_basis, 2);
        
        for i = 1:num_complement_vecs
            vec = U_complement_basis(:, i);
            if norm(vec) ~= 0
                U_complement = [U_complement, vec / norm(vec)];
            end
        end
        % --- End Fix 2 ---
        
        % Combine the two parts to form the full U matrix
        U = [U_partial, U_complement];
    else
        % If m <= r, U is already complete
        U = U_partial;
    end
    
    % --- Step 3: Simplify and Return ---
    % simplify() makes the symbolic output much cleaner
    U = simplify(U);
    Sigma = simplify(Sigma);
    V = simplify(V);
    
end