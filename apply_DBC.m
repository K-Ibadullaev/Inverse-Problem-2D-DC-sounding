function [A, b] = apply_DBC(A, b, sigma, x, y, global_idx, g, type)
    % Apply Dirichlet boundary conditions on linear system.
    %
    % SYNTAX
    %   [A, b] = apply_DBC(A, b, sigma, x, y, global_idx, g, type)
    %
    % INPUT PARAMETER
    %   A          ... Matrix representing the FD operator.
    %   b          ... RHS vector.
    %   sigma      ... Vector of cell parameter, i.e. conductivities
    %   x          ... Vector of mesh nodes in x.
    %   y          ... Vector of mesh nodes in y.
    %   global_idx ... Indices where value has to be applied.
    %   g          ... Dirichlet value at idx.
    %   type       ... Char denoting discretization type.
    %
    % OUTPUT PARAMETER
    %   A ... Updated operator matrix.
    %   b ... Updated rhs vector.
    %
    % REMARKS
    %   Always ensure to apply Dirichlet BC after all others!
    %
    % Mathias Scheunert, Sascha Weit 2022

    % Sanity checks.
    nx = numel(x);
    ny = numel(y);
    n = nx*ny;
    % Check data type and length.
    assert(isvector(g) && isvector(global_idx) && ...
           length(g) == length(global_idx));
    % Identify all boundary nodes (indices) and sort w.r.t. 'side'.
    all_global_idx = reshape(1:(n), nx, ny);
    bnd_global_idx = {all_global_idx(1, :).', ...   % at x_min all y
                      all_global_idx(end, :).', ... % at x_max all y
                      all_global_idx(:, 1), ...     % at y_min all x
                      all_global_idx(:, end)};      % at y_max all x
    % Check that given global_idx are boundary indices.
    assert(all(ismember(global_idx, vertcat(bnd_global_idx{:}))));

    % Get main diagonal (faster to re-assemble than to do A(sub2ind())!)
    [C0, ~, ~, ~, ~] = assemble_coefficients(x, y, sigma, type);
    C0 = C0(:);
    
    % Remove truncated FD expression and insert identity.
    A(global_idx, :) = 0;
    C0(global_idx) = 1;
    A = spdiags(C0, 0, A);

    % Insert Dirichlet value in system rhs vector.
    b(global_idx) = g;
end
