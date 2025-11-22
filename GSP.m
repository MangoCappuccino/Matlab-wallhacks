function [U, Q, dependent_indices] = GSP(V, tolerance)
% GSP Performs Gram-Schmidt orthogonalization with linear dependence detection
%
% Inputs:
%   V - Matrix where each column is a vector
%   tolerance - (optional) Threshold for dependence detection (default 1e-10)
%
% Outputs:
%   U - Orthogonal (unnormalized) vectors
%   Q - Orthonormal vectors
%   dependent_indices - Indices of linearly dependent columns in V
%
% Saved automatically in workspace as:
%   orthogonal_set_GSP
%   orthonormal_set_GSP
%   dependent_indices_GSP
%
% Example:
%   V = [1 0 1; 1 1 0; 0 1 1];
%   GSP(V);

    if nargin < 2
        tolerance = 1e-10;
    end

    [~, n] = size(V);
    U = [];
    Q = [];
    dependent_indices = [];

    fprintf('\n=== Gram-Schmidt Process ===\n\n');

    for j = 1:n
        fprintf('Processing vector v%d:\n', j);
        u_j = V(:, j);

        % Subtract projections onto previous orthogonal vectors
        for i = 1:size(U, 2)
            proj = (dot(V(:, j), U(:, i)) / dot(U(:, i), U(:, i))) * U(:, i);
            u_j = u_j - proj;
        end

        % Check for linear dependence
        if norm(u_j) < tolerance
            fprintf('  → v%d is LINEARLY DEPENDENT (ignored)\n\n', j);
            dependent_indices = [dependent_indices, j];
            continue;
        end

        % Add orthogonal and orthonormal vectors
        U = [U, u_j];
        q_j = u_j / norm(u_j);
        Q = [Q, q_j];

        fprintf('  → v%d added to orthogonal and orthonormal sets (‖u%d‖ = %.4f)\n\n', j, j, norm(u_j));
    end

    % Display results
    fprintf('=== Summary ===\n');
    fprintf('Independent vectors: %d\n', size(Q, 2));
    if isempty(dependent_indices)
        fprintf('✓ All vectors are linearly independent!\n\n');
    else
        fprintf('✗ Dependent vector indices: [ ');
        fprintf('%d ', dependent_indices);
        fprintf(']\n\n');
    end

    fprintf('=== Orthogonal Set (U) ===\n');
    disp(U);
    fprintf('=== Orthonormal Set (Q) ===\n');
    disp(Q);
    fprintf('=== Verification (Q^T * Q) ===\n');
    disp(Q' * Q);

    % Assign results to base workspace
    assignin('base', 'orthogonal_set_GSP', U);
    assignin('base', 'orthonormal_set_GSP', Q);
    assignin('base', 'dependent_indices_GSP', dependent_indices);

    fprintf('\nResults saved to workspace:\n');
    fprintf('  orthogonal_set_GSP\n');
    fprintf('  orthonormal_set_GSP\n');
    fprintf('  dependent_indices_GSP\n');
    fprintf('=============================================\n\n');
end