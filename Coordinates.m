function [coords, basisType, projection, lsError] = coordInBasis(B, v)
    % B: matrix whose columns are basis vectors
    % v: vector to express in basis B
    % coords: coordinates of v (or its projection) in B
    % basisType: 'Orthonormal', 'Orthogonal', or 'Not orthogonal'
    % projection: projection of v onto span(B)
    % lsError: norm(v - projection)
    
    tol = 1e-10;  % tolerance for checks
    
    % Check orthogonality
    G = B' * B;
    diagOnes = all(abs(diag(G) - 1) < tol);
    offDiagZeros = all(all(abs(G - diag(diag(G))) < tol));
    
    if diagOnes && offDiagZeros
        basisType = 'Orthonormal';
        coordsFull = B' * v;
    elseif offDiagZeros
        basisType = 'Orthogonal';
        coordsFull = (diag(1./diag(G)) * B') * v;
    else
        basisType = 'Not orthogonal';
        coordsFull = (B' * B) \ (B' * v); % least squares solution
    end
    
    % Compute projection
    projection = B * coordsFull;
    
    % Compute least squares error
    lsError = norm(v - projection);
    
    % Determine if v is in subspace
    if lsError < tol
        coords = coordsFull;
        fprintf('v is in the subspace. Coordinates computed.\n');
    else
        coords = coordsFull;
        fprintf(['v is NOT in the subspace.\n', ...
                 'Returning coordinates of the projection, projection vector, and least squares error.\n']);
    end
    
    % Display basis type
    fprintf('Basis type: %s\n', basisType);
end