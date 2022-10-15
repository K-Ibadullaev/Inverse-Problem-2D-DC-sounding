function [A, b] = assemble_system(type, x, y, sigma, I, global_idx)
    % Assemble linear system for FD continuity equation discretization.
    %
    % SYNTAX
    %   [A, b] = assemble_system(type, x, y, sigma, I, global_idx)
    %
    % INPUT PARAMETER
    %   type       ... Character, denoting the FD operator type
    %                  ['Laplace', 'BWT', 'DM2']
    %   x          ... Vector of mesh nodes in x
    %   y          ... Vector of mesh nodes in y
    %   sigma      ... Vector of cell parameter, i.e. conductivities
    %   global_idx ... Vector of (linear) node indices of source locations
    %   I          ... Vector of source strengths
    %
    % OUTPUT PARAMETER
    %   b ... Vector of system rhs
    %
    % Mathias Scheunert, Sascha Weit 2022
    
    nx = length(x);
    ny = length(y);
    n = nx*ny;
    
    % Get coupling coefficients.
    [C0, C1, C2, C3, C4] = assemble_coefficients(x, y, sigma, type);

    % Assemble system matrix.
    A_ = spdiags([C3(:), C1(:), C0(:), C2(:), C4(:)], ...
                  nx+[-nx, -1, 0, 1, nx], n, n+2*nx); % initialize as recangular matrix
    A = A_(:, nx+1:end-nx);                           % remove redundant columns
    
    % Assemble rhs.
    b = assemble_rhs(x, y, global_idx, I);
end