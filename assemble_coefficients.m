function [C0, C1, C2, C3, C4] = assemble_coefficients(x, y, k, type)
    % Assemble coupling coefficients for 2D FD discretizations.
    %
    % SYNTAX
    %   [C0, C1, C2, C3, C4] = assemble_coefficients(x, y, k, type)
    %
    % INPUT PARAMETER
    %   x    ... Vector of mesh nodes in x.
    %   y    ... Vector of mesh nodes in y.
    %   k    ... Vector of cell parameter.
    %   type ... Char denoting discretization type.
    %
    % OUTPUT PARAMETER
    %   C0 ... C4 ... Coupling coefficient vectors (matrix representation).
    %
    % Sascha Weit 2022

    % Interpolate parameter and calculate its derivatives.
    sig_point = get_sig_interpol(k, x, y);
    [sig_gradx, sig_grady] = get_sig_grad(k, x, y);

    % Enlarge parameter (incorporate mirror cells).
    sig_enlarged = zeros(size(sig_point)+2);
    sig_enlarged(2:end-1,2:end-1) = sig_point;
    sig_enlarged(1,:) = sig_enlarged(2,:);
    sig_enlarged(end,:) = sig_enlarged(end-1,:);
    sig_enlarged(:,1) = sig_enlarged(:,2);
    sig_enlarged(:,end) = sig_enlarged(:,end-1);

    % Get distances.
    f = diff(x(:));
    g = diff(y(:));
    f = [f(1); f(:); f(end)];
    g = [g(1); g(:); g(end)];
    [fi, gi] = ndgrid(f(2:end), g(2:end));
    [fni, gni] = ndgrid(f(1:end-1), g(1:end-1));
    
    % Get operator coefficients, 
    %   C2 for x+ i.e. f(1:end-1) and C1 for x- i.e. f(2:end), 
    %   C3 for y+ i.e. g(1:end-1) and C4 for y- i.e. g(2:end).
    switch type
        case 'BWT' % Brewitt-Taylor and Weaver
            C1 = (2*sig_point - sig_gradx .* fi) ./ (fni.*(fni+fi));
            C2 = (2*sig_point + sig_gradx .* fni) ./ (fi.*(fni+fi));
            C3 = (2*sig_point - sig_grady .* gi) ./ (gni.*(gni+gi));
            C4 = (2*sig_point + sig_grady .* gni) ./ (gi.*(gni+gi));
        case 'L'  % Laplace
            C1 = (2*sig_point) ./ (fni.*(fni+fi));
            C2 = (2*sig_point) ./ (fi.*(fni+fi));
            C3 = (2*sig_point) ./ (gni.*(gni+gi));
            C4 = (2*sig_point) ./ (gi.*(gni+gi));
        case 'DM2' % Dey and Morrison 2
            C1 = (sig_enlarged(1:end-2, 2:end-1) + sig_enlarged(2:end-1, 2:end-1)) ./ ...
                    (fni .* (fni + fi));
            C2 = (sig_enlarged(3:end, 2:end-1) + sig_enlarged(2:end-1, 2:end-1)) ./ ...
                    (fi .* (fni + fi));
            C3 = (sig_enlarged(2:end-1, 1:end-2) + sig_enlarged(2:end-1, 2:end-1)) ./ ...
                    (gni .* (gni + gi));
            C4 = (sig_enlarged(2:end-1, 3:end) + sig_enlarged(2:end-1, 2:end-1)) ./ ...
                    (gi .* (gni + gi));
        case 'Gx'   % Gradient in x-direction
            C1 = (-fi) ./ (fni.*(fni+fi));
            C2 = (fni) ./ (fi.*(fni+fi));
            C3 = 0*C1;
            C4 = 0*C1;
        case 'Gy'   % Gradient in y-direction
            C3 = (-gi) ./ (gni.*(gni+gi));
            C4 = (gni) ./ (gi.*(gni+gi));
            C1 = 0*C3;
            C2 = 0*C3;
    end
    C0 = -(C1 + C2 + C3 + C4);
end