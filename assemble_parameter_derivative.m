function O = assemble_parameter_derivative(x, y, type)
    % Set up sparse 1st or 2nd order derivative acting on cells in 2D mesh.
    %
    % SYNTAX
    %   O = assemble_parameter_derivative(x, y, type)
    %
    % INPUT PARAMETER
    %   x    ... Vector [nx, 1] of mesh nodes in x direction.
    %   y    ... Vector [ny, 1] of mesh nodes in y direction.
    %   type ... Char denoting desired finite difference scheme.
    %
    % OUTPUT PARAMETER
    %   O ... Matrix [nx-1 x ny-1] representing the 1st derivative for 
    %         cells in mesh(x, y).
    %
    % REMARKS
    %   Using mesh nodes to define 'inner mesh' living on cell midpoints.
    %   Boundary conditions:
    %   - hom. Neumann  ... ensure that parameter at boundary continues

    % Get cell midpoints.
    x =(x(1:end-1) + x(2:end)) / 2;
    y =(y(1:end-1) + y(2:end)) / 2;
    nx = length(x);
    ny = length(y);

    % Get boundaries.
    [idx_x, idx_y] = ndgrid(1:nx, 1:ny);
    idx_x_min = idx_x == 1;
    idx_x_max = idx_x == nx;
    idx_y_min = idx_y == 1;
    idx_y_max = idx_y == ny;
    bnd_idx = find(idx_x_min | idx_x_max | idx_y_min | idx_y_max);

    % Set up operator.
    O = assem_2D_operator(type, x, y);

    % Apply hom. Neumann BC.
    O = add_hom_Neumann_BC(type, O, x, y, bnd_idx);
end

%% Helper.

function O = assem_1D_operator(type, x, dx)
    % Sparse central difference operator in 1D.
    %
    % INPUT PARAMETER
    %   type ... Char denoting desired finite difference scheme.
    %   x    ... Vector of mesh nodes.
    %   dx   ... Vector of step sizes.
    %
    % OUTPUT PARAMETER
    %   O ... Operator matrix.
    
    % Initialize.
    n = length(x);
    dx_ = [dx(1), dx, dx(end)]; % mirror mesh spacing at end nodes
    
    % Assemble.
    weights_lines = 1 ./ ((dx_(1:end-1) .* dx_(2:end)) .* ...
                          (dx_(1:end-1) + dx_(2:end)));
    switch type
        case 'Gradient'
            % 1st order central differences.
            weights_cols = [-dx_(2:end).^2; ...
                             dx_(2:end).^2 - dx_(1:end-1).^2; ...
                             dx_(1:end-1).^2].';            
                        
        case 'Laplace'
            % 2nd-order central differences.
            weights_cols = 2 * [  dx_(2:end); ...
                                -(dx_(1:end-1) + dx_(2:end)); ...
                                  dx_(1:end-1)].';
            
        otherwise
            error('Unkown operator type.');
    end
    % Note:use trick to correctly place columns of d_cols on diagonals of 
    % D, as spdiags([col1, col2, ...], [row_idx1, row_idx2, ...], ,..)
    % for square matrices!
    O = weights_lines(:) .* spdiags(fliplr(weights_cols), ...
                                    [-1, 0, 1], n, n).';
end

function O = assem_2D_operator(type, x, y)
    % Sparse central difference operator at mesh nodes z with step size dz.
    %
    % INPUT PARAMETER
    %   type ... Char denoting desired finite difference scheme.
    %   x    ... Vector of mesh nodes in x.
    %   y    ... Vector of step sizes in y.
    %
    % OUTPUT PARAMETER
    %   O ... Operator matrix.
    
    % Initialize 1D operators.
    Ox_ = assem_1D_operator(type, x, diff(x));
    Oy_ = assem_1D_operator(type, y, diff(y));

    % Map (local) 1D operators to global.
    switch type
        case 'Gradient'
            O = [kron(eye(length(y)), Ox_); kron(Oy_, eye(length(x)))];
        case 'Laplace'
            O = kron(eye(length(y)), Ox_) + kron(Oy_, eye(length(x)));
        otherwise
            error('Unsupported typer of derivative.');
    end
end

function O = add_hom_Neumann_BC(type, O, x, y, global_idx)
    % Incorporate homogeneous Neumann values by adjusting O.
    %
    % INPUT PARAMETER
    %   O          ... Matrix representing the FD operator.
    %   x          ... Vector of mesh nodes in x.
    %   y          ... Vector of mesh nodes in y.
    %   global_idx ... Index where value has to be applied.
    %
    % OUTPUT PARAMETER
    %   O ... Updated operator matrix.

    % Sanity checks.
    [n, m] = size(O);
    nx = numel(x);
    ny = numel(y);
    all_global_idx = reshape(1:(nx*ny), nx, ny);
    bnd_global_idx = {all_global_idx(1, :).', ...
                      all_global_idx(end, :).', ...
                      all_global_idx(:, 1), ...
                      all_global_idx(:, end)};
    assert(all(ismember(global_idx, vertcat(bnd_global_idx{:}))));
	bnd_idx = cellfun(@(i) ismember(global_idx, i), bnd_global_idx, ...
                      'UniformOutput', false);
    switch type
        case 'Gradient'
            assert(n == 2*m);
            % Apply BC on both parts.
            O = mat2cell(O, [m, m], m);
        case 'Laplace'
            assert(n == m)
            O = {O};
        otherwise
                error('Derivative not supported.');
    end
    
    % Enlarge spacing vectors.
    dx = diff(x);
    dy = diff(y);
    d_ = {[dx(1), dx, dx(end)];
          [dy(1), dy, dy(end)]};

    % Set rhs value depending on operator.
    % Loop over bnd parts.
    for s = 1:length(bnd_idx)
        if s <= 2   % w.r.t x
            grad_idx = 1;
            shift = 1;
        else        % w.r.t y
            grad_idx = 2;
            shift = nx;
        end
        cur_global_idx = global_idx(bnd_idx{s});
        switch type
            case 'Gradient'                          
                % Adjust entry related bnd grid node.
                if mod(s, 2) == 1    % bnd at x/y min
                    O{grad_idx}(sub2ind([m, m], cur_global_idx, cur_global_idx)) = ...
                        -1 / (2*d_{grad_idx}(1));
                else                 % bnd at x/y max
                    O{grad_idx}(sub2ind([m, m], cur_global_idx, cur_global_idx)) =  ...
                        1 / (2*d_{grad_idx}(end));
                end

            case 'Laplace'
                % Set rhs entry adjust neighbour entry related to bnd 
                % grid node.
                if mod(s, 2) == 1
                    O{1}(sub2ind([m, m], cur_global_idx, cur_global_idx+shift)) = ...
                        2*O{1}(sub2ind([m, m], cur_global_idx, cur_global_idx+shift));
                else
                    O{1}(sub2ind([m, m], cur_global_idx, cur_global_idx-shift)) = ...
                        2*O{1}(sub2ind([m, m], cur_global_idx, cur_global_idx-shift));
                end
        end
    end

    % Map Neumann values.
    O = cell2mat(O);
end
