function T = assemble_tensor(type, x, y, bnd_D, bnd_N)
    % Assembles tensor by concatenation of matrix assemblys.
    %
    % SYNTAX
    %   T = assemble_tensor(type, x, y, bnd_D, bnd_N)
    %
    % INPUT PARAMETER
    %   type           ... Character, denoting the FD operator type
    %                      ['Laplace', 'BWT', 'DM2']
    %   x              ... Vector of mesh nodes in x
    %   y              ... Vector of mesh nodes in y
    %   bnd_N, bnd_D   ... Vectors of vertex indices at which
    %                      Neumann/Dirichlet boundaries are defined
    %
    % OUTPUT PARAMETER
    %   T ... Matrix [n_cell*n_dof x n_dof] representing the derivative 
    %         tensor d_A/d_m [n_dof x n_dof x n_cell] where all tensor  
    %         slices are vertically concatenated
    %
    % REMARKS
    %   The routine exploits the linear dependence of the system matrix A
    %   from model parameter m:
    %   - the derivative w.r.t. m_i sets all entries of A which depend on
    %     m_j (i~=j) to zero
    %   - the derivative w.r.t. m_i sets all entries containing k*m_i to k
    %   Hence, assembling the derivative w.r.t. m_i can be done by
    %   assembling the system matrix with the ith unit vector as model
    %   parameter
    %
    % Mathias Scheunert, Sascha Weit 2022

    % Loop over parameter.
    nc = (length(x)-1) * (length(y)-1);
    T = cell(nc, 1);
    assem_handle = @(k) assem_DC_system(k, type, x, y, bnd_D, bnd_N);
    fprintf('Assemble T\n');
    fprintf(sprintf('Loop over %d parameter: ', nc));
    count = 0;
    for i = 1:nc
        % Assemble system for unit parameter vector.
        T{i} = assem_handle(get_unit_vector(i));
        count = count + 1; % add some verbosity
        if count >= round(nc/20)-1
            fprintf('.');
            count = 0;
        end
    end
    T = cell2mat(T);
    fprintf(' done.\n');

    %% Helper
    
    function e = get_unit_vector(i)
        % Construct unit vector. 
        e = zeros(nc, 1);
        e(i) = 1;
    end

    function A = assem_DC_system(k, type, x, y, bnd_D, bnd_N)
        % Assembling of the DC fwd problem lhs.
        A = assemble_system(type, x, y, k, [], []);
        A = apply_tensor_NBC(A, k, x, y, bnd_N, type);
        A = apply_tensor_DBC(A, bnd_D);
    end

    function A = apply_tensor_NBC(A, k, x, y, global_idx, type)
        % Apply Neumann boundary conditions on linear system lhs.
        % Identify all boundary nodes (indices) and sort w.r.t. 'side'.
        nx = numel(x);
        ny = numel(y);
        n = nx*ny;
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
    
        % Get system coefficients.
        [~, C1, C2, C3, C4] = assemble_coefficients(x, y, k, type);
    
        % Loop over bnd parts.
        for s = 1:length(bnd_idx)
            if s <= 2   % w.r.t x
                shift = 1;
                coeff_LHS = C1(:) + C2(:);
            else        % w.r.t y
                shift = nx;
                coeff_LHS = C3(:) + C4(:);
            end
            cur_global_idx = global_idx(bnd_idx{s});       
            if mod(s, 2) == 1 % x/y min
                % Adjust neighbour entries of system related to bnd grid node  
                % and direction.
                A(sub2ind([n, n], cur_global_idx, cur_global_idx+shift)) = ...
                    coeff_LHS(cur_global_idx); 
            else              % x/y max
                A(sub2ind([n, n], cur_global_idx, cur_global_idx-shift)) = ...
                    coeff_LHS(cur_global_idx); 
            end
        end
    end

    function A = apply_tensor_DBC(A, global_idx)      
        % Replace truncated FD expression by zeros.
        A(global_idx, :) = 0;
    end
end
