function plot_at_node(F, X, Y, type)
    % Plots F defined at nodes X,Y for tensor-product mesh.
    %
    % SYNTAX
    %   plot_at_node(F, X, Y, type)
    %
    % INPUT PARAMETER
    %   F ... Matrix [nx, ny] of node values.
    %   X ... Matrix [nx, ny] of mesh nodes in x.
    %   Y ... Matrix [nx, ny] of mesh nodes in y.
    %
    % OPTIONAL PARAMETER
    %   type ... Char, denoting a plot type.
    %
    % REMARKS
    %   If no interpolation is used, given values F (defined at nodes X and 
    %   Y) are plotted at fictional cells that are defined around the
    %   nodes.
    %   The mesh nodes will be visualized as small points.
    %
    % Mathias Scheunert 2022
    
    assert(all(size(F) == size(X)) && all(size(F) == size(Y)));
    if nargin > 3
        switch type
            case 'interp'
                % Plot interpolation between nodes.
                % No enlargement of the mesh is required here, as (corner) 
                % nodes of the plot will have its true value.
                plot_handle = pcolor(X, Y, F);
                plot_handle.FaceColor = 'interp';
                
            otherwise
                error('Unsupported plot type.');
        end
    else
        % Add dummy entries.
        F_ = NaN(size(F) + 1);
        F_(1:end-1, 1:end-1) = F;

        % Redefine mesh. 
        % Rectangular cells around the mesh nodes are constructed (using 
        % the given mesh spacings) such that the resulting cells will have 
        % the value of the nodes located at its centers.
        dx = diff(X(:, 1)).';
        dy = diff(Y(1, :));
        dx = [dx(1), dx, dx(end)]; % mirror boundary cells
        dy = [dy(1), dy, dy(end)];
        dx_ = [dx; dx] ./ 2;
        dy_ = [dy; dy] ./ 2;
        dx_ = dx_(:);
        dy_ = dy_(:);
        dx_ = reshape(dx_(2:end-1), 2, []);
        dy_ = reshape(dy_(2:end-1), 2, []);
        x_ = [(X(1)-dx_(1)), (X(1)-dx_(1)) + cumsum(sum(dx_, 1))];
        y_ = [(Y(1)-dy_(1)), (Y(1)-dy_(1)) + cumsum(sum(dy_, 1))];
        [X_, Y_] = ndgrid(x_, y_);

        % Plot redefined mesh with F as polygon value.
        plot_handle = pcolor(X_, Y_, F_);
    end
    plot_handle.EdgeColor = 'none';
    xlabel('x in [m]');
    ylabel('y in [m]');
    axis('equal');
    colorbar();
    hold on
        plot(X, Y, '.k', 'MarkerSize', 1);   % add nodes
    hold off
end