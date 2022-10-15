function A = ttv(T, v, dim)
    % Multiply tensor with vector along dim.
    %
    % SYNTAX
    %   A = ttv(T, v, dim)
    %
    % INPUT
    %   T   ... Matrix [n_cell*n_dof x n_dof] representing a tensor 
    %           [n_dof x n_dof x n_cell] 
    %   v   ... Vector [n_cell/n_dof x 1]
    %   dim ... Scalar denoting the dimension at which v should be
    %           multiplied at T
    %
    % OUTPUT
    %   A ... Resulting matrix
    %
    % Mathias Scheunert 2022

    % Get tensor sizes.
    [i, j] = size(T);
    dim_T = [j, j, i/j];
    assert(round(dim_T(3)) == dim_T(3));
    assert(length(v) == dim_T(dim));

    % Do multiplication.
    if dim == 1
        error('not supported');
    elseif dim == 2
        A = reshape(T*v, dim_T(1), dim_T(3)); 
    elseif dim == 3
        T_ = reshape(T.', [], dim_T(3));
        A = reshape(T_*v, dim_T(2), dim_T(1)).';
    else
        error('dim exceed tensor size.');
    end
end