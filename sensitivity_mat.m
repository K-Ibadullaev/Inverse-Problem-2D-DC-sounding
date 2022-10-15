function S = sensitivity_mat(type, x, y, X, Y, x_a, x_b, x_m, x_n, y_abmn, I_A, bnd_idx_N, bnd_idx_D, k, T, param)
         
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
         % T ... tensor
         % param ... parameters of grid, provided as vector

         % output S ... sensitivity matrix

         % solve direct problem 
        % solve direct problem 
         [~,A, phi] = fwp_dc_problem_enlarged(type, x, y, X, Y, x_a, x_b, x_m, x_n, y_abmn, I_A, bnd_idx_N, bnd_idx_D, k, param);

         % unique positions in x_a, x_b
         x_a_u = unique(x_a);
         x_b_u = unique(x_b);
         
         nm = numel(x_a_u);

         % assemble Jacobian

         for i = 1:nm
             L = ttv(T,phi{i},2);
             J{i} = (-A{i}\L);
         end


         for i = 1:numel(x_a)
             
             idx = find(x_a(i) == x_a_u);

             M_idx = find(X == x_m(i) & Y == y_abmn(i));
             N_idx = find(X == x_n(i) & Y == y_abmn(i));

             q = zeros(1, length(phi{idx}));
             q(M_idx) = k(i) * 1;
             q(N_idx) = -k(i) * 1;

             S(i,:) = q*J{idx};

         end
         
end