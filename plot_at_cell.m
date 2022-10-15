function plot_at_cell(K, X, Y)
    % Plots K defined at cells between nodes X,Y for tensor-product mesh.
    %
    % SYNTAX
    %   plot_at_cell(K, X, Y)
    %
    % INPUT PARAMETER
    %   K ... Matrix [nx-1, ny-1] of cell values.
    %   X ... Matrix [nx, ny] of mesh nodes in x.
    %   Y ... Matrix [nx, ny] of mesh nodes in y.
    %
    % Mathias Scheunert 2022
    
    assert(all(size(K) == size(X)-1) && all(size(K) == size(Y)-1));
    
    % Add dummy entries.
    K_ = NaN(size(K) + 1);
    K_(1:end-1, 1:end-1) = K;
    
    % Plot mesh with K as polygon value.
    plot_handle = pcolor(X, Y, K_);
    plot_handle.EdgeAlpha = 0.05;
    xlabel('x in [m]');
    ylabel('y in [m]');
    axis('equal');
    colorbar();
end