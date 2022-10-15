function d = fwp_dc_problem(type, x, y, X, Y, x_a, x_b, x_m, x_n, y_abmn, I_A, bnd_idx_N, bnd_idx_D, k, param)
         % calculates forward problem solution for parameter param (conductivity)

         % input
         % type ... discretization scheme
         % x, y ... x and y coordinates
         % X, Y ... grid
         % x_a ... positions of current electrode A
         % x_b ... positions of current electrode B
         % x_m ... position of measurement electrode m
         % x_n ... position of measurement electrode n
         % I_A ... current strength in A
         % bnd_idx_N ... index of boundary values for Neumann
         % bnd_idx_D ... index of boundary values for Dirichlet
         % k ... configuration factors k depending on geometry
         % param ... parameters of grid, provided as vector

         % output d: observation data d

         %% measurement setup
        
         % unique positions in x_a, x_b
         x_a_u = unique(x_a);
         x_b_u = unique(x_b);
        
         % number of current dipoles
         nm = numel(x_a_u);
        
         I_B = -I_A; % in A
        
         % Define boundary values (homogeneous BC at both boundary parts).
         bnd_g_N = zeros(length(bnd_idx_N), 1);
         bnd_g_D = zeros(length(bnd_idx_D), 1);

         %% calculate A, b for linear system of equations and solve it for potential phi
         %  do the observation at electrodes M, N
         %  calculate rho_a using observation operator q
         %  assemble the Jacobian matrix J and sensitivity matrix S
        
         for i = 1:nm
        
            % indices of source electrode
            A_idx = find(X == x_a_u(i) & Y == y_abmn(i));
            B_idx = find(X == x_b_u(i) & Y == y_abmn(i));
        
            % Assemble system and apply boundary conditions.
            % only the rhs changes, so A can be left out for the second source
            [A{i}, b_a] = assemble_system(type, x, y, param, I_A, A_idx);
            [~, b_b] = assemble_system(type, x, y, param, I_B, B_idx);
        
            % use linearity 
            b = b_a + b_b;
        
            [A{i}, b] = apply_NBC(A{i}, b, param, x, y, bnd_idx_N, bnd_g_N, type);
            [A{i}, b] = apply_DBC(A{i}, b, param, x, y, bnd_idx_D, bnd_g_D, type); 
        
            % solve system
            phi{i} = full(A{i}\b);
        
         end
        
         %% Observation 
        
         for i = 1:numel(x_m)
        
            idx = find(x_a(i) == x_a_u);
        
            M_idx = find(X == x_m(i) & Y == y_abmn(i));
            N_idx = find(X == x_n(i) & Y == y_abmn(i));
        
            q = zeros(1, length(phi{idx}));
            q(M_idx) = k(i) * 1;
            q(N_idx) = -k(i) * 1;
        
            % extract apparent resistivities
            d(i,1) =  q * phi{idx}; 
        
         end
end