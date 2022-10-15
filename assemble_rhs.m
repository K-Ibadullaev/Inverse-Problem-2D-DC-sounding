function b = assemble_rhs(x, y, global_idx, I)
    % Build rhs vector.
    %
    % SYNTAX
    %   b = assemble_rhs(x, y, global_idx, I)
    %
    % INPUT PARAMETER
    %   type       ... Character, denoting the FD operator type
    %   x          ... Vector of mesh nodes in x
    %   y          ... Vector of mesh nodes in y
    %   global_idx ... Vector of (linear) node indices of source locations
    %   I          ... Vector of source strengths
    %
    % OUTPUT PARAMETER
    %   b ... Vector of system rhs
    %
    % Mathias Scheunert, Sascha Weit 2022
        
    % Initialize.
    n = length(x) * length(y);
    b = sparse([], [], [], n, 1, 0);
    if nargin > 3
        % Add dummy spacings.
        f = diff(x(:)).'; 
        g = diff(y(:)).';
        [Fl, Gl] = ndgrid([0, f], [0, g]);
        [Fr, Gr] = ndgrid([f, 0], [g, 0]);
        
        % Predefine source volume (use shifted index due to dummies).
        V = (Fl(global_idx) + Fr(global_idx)) .* ...
            (Gl(global_idx) + Gr(global_idx)) ./ 4;
        
        % Add source term.
        b(global_idx) = -I(:) ./ V(:);
    end
end
