function taylor_test_default(fwp_handle, s, J, Js2Jm, s2m, m2s)
    % Apply taylor test to validate given Jacobian matrix J.
    %
    % SYNTAX
    %   taylor_test_default(fwp_handle, s, J, Js2Jm, s2m, m2s)
    %
    % INPUT PARAMETER
    %   fwp_handle ... Handle to forward problem solution, taking model
    %                  parameter vector s [n_cell x 1] and providing data
    %                  d [n_obs x 1].
    %   s          ... Model parameter vector [n_cell x 1]
    %   d          ... d(s) data vector [n_obs x 1], depending on s
    %   J          ... J(s) sensitivity matrix [n_obs x n_cell], depending 
    %                  on s
    %
    % OPTIONAL INPUT PARAMETER
    %   Js2Jm      ... Handle to parameter transformation for sensitivity 
    %                  matrix, taking sensitivity matrix [n_obs x n_cell] 
    %                  w.r.t. non-transformed model s and the vector 
    %                  [n_cell x 1] s and providing the transformed 
    %                  sensitivity matrix [n_obs x n_cell] 
    %   s2m        ... Handle to transform model parameter s, taking vector
    %                  s [n_cell x 1] and provide transformed vector m 
    %                  [n_cell x 1]
    %   m2s        ... Handle to back-transform model parameter m, taking 
    %                  vector m [n_cell x 1] and provide original vector s 
    %                  [n_cell x 1]
    %
    % Mathias Scheunert 2022

    % Check input arguments.
    if nargin < 4
        % Define default transformation handles.
        Js2Jm = @(J, s) J;
        s2m = @(s) s;
        m2s = @(m) m;
    else
        assert(nargin == 6);
    end

    % Transform.
    m = s2m(s);
    J = Js2Jm(J, s);

    % Initialize.
    fprintf('Taylor-test, calculate: ');
    dm = randn(size(m));
    dm = dm(:)./norm(dm);
    h_e_min = -7;
    h_e_max = -1;
    h = fliplr(logspace(h_e_max, h_e_min, (h_e_max-h_e_min+1))).';
    d_0 = size(h);
    d_1 = size(h);
    d = fwp_handle(m2s(m));
    
    % Calculate norms.
    for i = 1:length(h)
        mod_curr = m(:) + h(i).*dm;
        d_curr = fwp_handle(m2s(mod_curr)); % beware back-trafo of m!
        d_0(i) = norm(d_curr - d);
        d_1(i) = norm(d_curr - d - J * h(i)*dm);
    end
    fprintf(' done.\n');
    
    % remove NaN
    nan_true = isnan(d_1) & isnan(d_0); 
    h_plot = h(~nan_true);
    d_0_plot = d_0(~nan_true);
    d_1_plot = d_1(~nan_true);
    
    % Plot behaviour.
    loglog(h_plot,h_plot./h_plot(end));
    hold on
    loglog(h_plot, (h_plot.^2)/(h_plot(end)^2), '-r');
    loglog(h_plot, d_0_plot./d_0_plot(end),     'xk');
    loglog(h_plot, d_1_plot./d_1_plot(end),     'ok');
    set(gca, 'Xdir', 'reverse');
    legend('O(h)', 'O(h^2)', '0-order', '1-order');
    xlabel('h');
    ylabel('Residual norm(h)')
    title('Taylor-Test');
    hold off
end