function [A, b] = apply_NBC(A, b, sigma, x, y, global_idx, g, type)
    % Apply Neumann boundary conditions on linear system.
    %
    % SYNTAX
    %   [A, b] = apply_NBC(A, b, sigma, x, y, global_idx, g, type)
    %
    % INPUT PARAMETER
    %   A          ... Matrix representing the FD operator.
    %   b          ... RHS vector.
    %   sigma      ... Vector of cell parameter, i.e. conductivities
    %   x          ... Vector of mesh nodes in x.
    %   y          ... Vector of mesh nodes in y.
    %   global_idx ... Indices where value has to be applied.
    %   g          ... Neumann value at idx.
    %   type       ... Char denoting discretization type.
    %
    % OUTPUT PARAMETER
    %   A ... Updated operator matrix.
    %   b ... Updated rhs vector.
    %
    % REMARKS
    %   Always ensure to apply Neumann BC before all others!
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
    
    % Match each given (global) boundary index with the corresponding side
    % of the model domain.
	bnd_idx = cellfun(@(i) ismember(global_idx, i), bnd_global_idx, ...
                      'UniformOutput', false);
    
    % Enlarge spacing vectors.
    dx = diff(x);
    dy = diff(y);
    d_ = {[dx(1); dx(:); dx(end)];
          [dy(1); dy(:); dy(end)]};

    % Get system coefficients.
    [~, C1, C2, C3, C4] = assemble_coefficients(x, y, sigma, type);

    % Loop over bnd parts.
    for s = 1:length(bnd_idx)
        if s <= 2   % w.r.t x
            grad_idx = 1;
            shift = 1;
            coeff_LHS = C1(:) + C2(:);
            coeff_RHS = {C1(:), C2(:)};
        else        % w.r.t y
            grad_idx = 2;
            shift = nx;
            coeff_LHS = C3(:) + C4(:);
            coeff_RHS = {C3(:), C4(:)};
        end
        cur_bnd_idx = bnd_idx{s};
        cur_global_idx = global_idx(cur_bnd_idx);       
        if mod(s, 2) == 1 % x/y min
            % Adjust rhs related to bnd grid node and direction.
            b(cur_global_idx) = b(cur_global_idx) + 2*g(cur_bnd_idx) .* ...
                (d_{grad_idx}(1)) .* coeff_RHS{grad_idx}(cur_global_idx);
            % Adjust neighbour entries of system related to bnd grid node  
            % and direction.
            A(sub2ind([n, n], cur_global_idx, cur_global_idx+shift)) = ...
                coeff_LHS(cur_global_idx); 
        else              % x/y max
            b(cur_global_idx) = b(cur_global_idx) - 2*g(cur_bnd_idx) .* ...
                (d_{grad_idx}(end)) .* coeff_RHS{grad_idx}(cur_global_idx);
            A(sub2ind([n, n], cur_global_idx, cur_global_idx-shift)) = ...
                coeff_LHS(cur_global_idx); 
        end
    end
end