function [sig_gradx, sig_grady] = get_sig_grad(sig, x, y)
    % Get (volume weighted) derivatives of sig at mesh nodes x, y.
    %
    % SYNTAX
    %   [sig_gradx, sig_grady] = get_sig_grad(sig, x, y)
    % 
    % INPUT PARAMETER
    %   sig ... Vector/matrix of parameter defined at cells
    %   x   ... Vector of mesh nodes in x
    %   y   ... Vector of mesh nodes in y
    %
    % OUTPUT PARAMETER
    %   sig_gradx ... Vector of parameter derivative in x direction
    %   sig_grady ... Vector of parameter derivative in y direction
    %
    % Mathias Scheunert, Sascha Weit 2022

    % Get matrix representation of mesh spacings.
    f = diff(x); 
    g = diff(y);
    % Mirror mesh spacing: in x                  in y           direction
    [F, G] = ndgrid([f(1); f(:); f(end)], [g(1); g(:); g(end)]);
    
    % May transform.
    if isvector(sig)
        sig = reshape(sig, length(x)-1, length(y)-1);
    end

    % Mirror outermost (cell)parameters.
    sig_enlarged = zeros(size(sig)+2);
    sig_enlarged(2:end-1,2:end-1) = sig;
    sig_enlarged(1, :) = sig_enlarged(2, :);
    sig_enlarged(end, :) = sig_enlarged(end-1, :);
    sig_enlarged(:, 1) = sig_enlarged(:, 2); % corners are included
    sig_enlarged(:, end) = sig_enlarged(:, end-1);

    % Build auxiliary S.
    S = sig_enlarged .* F .* G;

    % Get gradients:
    % in x (w.r.t. i)
    sig_gradx = (S(2:end, 2:end) + S(2:end, 1:end-1)) ./ ...
                    (F(2:end, 2:end) .* (G(1:end-1, 1:end-1) + G(2:end, 2:end))) - (S(1:end-1, 2:end) + S(1:end-1, 1:end-1)) ./ ...
                    (F(1:end-1, 1:end-1) .* (G(1:end-1, 1:end-1) + G(2:end, 2:end)));
    sig_gradx = 2*sig_gradx ./ (F(2:end, 2:end) + F(1:end-1, 1:end-1));
    % y direction (w.r.t. j)
    sig_grady = (S(1:end-1, 2:end) + S(2:end, 2:end)) ./ ...
                    (G(2:end, 2:end) .* (F(1:end-1, 1:end-1) + F(2:end, 2:end))) - (S(1:end-1, 1:end-1) + S(2:end, 1:end-1)) ./ ...
                    (G(1:end-1, 1:end-1) .* (F(1:end-1, 1:end-1) + F(2:end, 2:end)));
    sig_grady = 2*sig_grady ./ (G(2:end, 2:end) + G(1:end-1, 1:end-1));
end
