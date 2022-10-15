%% Short summary
% The most precise method in the given setting is Levenberg-Marquardt. 
% More sophisticated damped LS method has shown a bit worse performance, though
% C='grad' demonstrated quite accurate results in only 23 iterations,
% which is 2 times less iterations then for Levenberg
% Unfortunately, for C='Laplace' it yields matrix is close to singular
% and condition number is extremely small

% Also worth to note changed alpha and alpha_min for Line search
% The golden section ratio for alpha and higher alpha_min led to 
% faster convergence, though it's necessary that they don't lead to
% skipping the of global minimum

% Please, run the script by section.

%%
clc;
clear;

%% load data
load data.mat
%% set up forward problem


% spacing
dx = 1.5;

% inhomogeneous grid
x_h = (-30:dx:30);
x_h1 = logspace(log10(50), log10(2000), 30);
x_h2 = -flip(x_h1);
x = [x_h2 x_h x_h1];
 
y_h = -30:dx:0;
y_h2 = -logspace(log10(2000),log10(50), 15);
y = [y_h2 y_h];



[X, Y] = ndgrid(x,y);

%% model parameter

xc = x(1:end-1) + diff(x)./2;
yc = y(1:end-1) + diff(y)./2;
[Xc, Yc] = ndgrid(xc, yc);

% parameter - conductivity
P = zeros(size(Xc))+ 1.0361e-04;%0.0001;%1.0361e-04


param = P(:);

%% measurement 

% electrode positions
x_e = data.ele(:,1);
y_e = data.ele(:,2);

% positions of A, B, M, N 
x_a = x_e(data.measurement(:,1));
x_b = x_e(data.measurement(:,2));
x_m = x_e(data.measurement(:,3));
x_n = x_e(data.measurement(:,4));

% unique positions for x_a, x_b
x_a_u = unique(x_a);
x_b_u = unique(x_b);

% number of dipoles
dp_nm = numel(x_a_u);

y_abmn = y_e(data.measurement(:,1)); % same for all electrodes


I_A = 1; % in A
I_B = -1; % in B

% measured data

rho_measured = data.measurement(:,6);

%% indices for boundary conditions

% Get linear node indices of all four boundaries.
idx_x_min = X == min(x);
idx_x_max = X == max(x);
idx_y_min = Y == min(y);
idx_y_max = Y == max(y); % air-earth surface
bnd_idx_Neum = find(idx_y_max);
bnd_idx_Dch = find(idx_x_min | idx_x_max | idx_y_min);

% Define boundary values (homogeneous BC at both boundary parts).
bnd_g_Neum = zeros(length(bnd_idx_Neum), 1);
bnd_g_Dch = zeros(length(bnd_idx_Dch), 1);

type = 'BWT'; % set  discretization type by Brewit-Tayor & Weaver

%% assemble tensor 
T = assemble_tensor(type,x,y,bnd_idx_Dch, bnd_idx_Neum);

%% calculate A, b for linear system of equations and solve it for potential phi
%  do the observation at electrodes M, N
%  calculate rho_a using measurement operator q
%  assemble the Jacobian matrix J and sensitivity matrix S
% profile on
for i = 1:dp_nm

    % indices of source electrode
    A_idx = find(X == x_a_u(i) & Y == y_abmn(i));
    B_idx = find(X == x_b_u(i) & Y == y_abmn(i));

    % Assemble system and apply boundary conditions.
    % only rho changes
    [A{i}, b_a] = assemble_system(type, x, y, P(:), I_A, A_idx);
    [~, b_b] = assemble_system(type, x, y, P(:), I_B, B_idx);

    % use linearity 
    b = b_a + b_b;

    [A{i}, b] = apply_NBC(A{i}, b, P(:), x, y, bnd_idx_Neum, bnd_g_Neum, type);
    [A{i}, b] = apply_DBC(A{i}, b, P(:), x, y, bnd_idx_Dch, bnd_g_Dch , type); 

    % solve system
    phi{i} = full(A{i}\b);

end
% profile viewer
%% Observation 


for i = 1:dp_nm
    L = ttv(T,phi{i},2);
    J{i} = -(A{i}\L);  
end

for i = 1:numel(x_m)

    idx = find(x_a(i) == x_a_u);

    M_idx = find(X == x_m(i) & Y == y_abmn(i));
    N_idx = find(X == x_n(i) & Y == y_abmn(i));

    k(i) = data.measurement(i,5); % configuration factor

    q = zeros(1, length(phi{idx}));
    q(M_idx) = k(i) * 1;
    q(N_idx) = -k(i) * 1;
%     save measurement operator for all configurations
    Q(:,i) = q;
    % extract apparent resistivities
    rho_a(i) = q * phi{idx}; 

    % find reference point for plot
    x_ref(i) = (x_a(i) + x_n(i))/2;  

    %%%% assemble sensitivity matrix %%%%

    S{i} = q * J{idx};
end
%% plot

elstp = figure;

for i = 1:dp_nm
    plot_at_node(reshape(phi{i}, length(x), length(y)), X, Y);
    title(['fwp solution of configuration ' num2str(i)]);
    a = colorbar;
    a.Label.String = 'electrical potential [V]';
    ax = gca;
    ax.CLim = [-5000 5000];
    xlabel('x [m]');
    ylabel('y [m]');
    xlim([-20, 20]);
    ylim([-10, 0]);
    pause(0.5);
    set(gcf,'color','w');
%     saveas(elstp,sprintf('norms/el_setup_%d.png',i));
end
%% pseudoplot

ffsct=figure;
t = tiledlayout(1,2);
% 
nexttile
plot_pseudosection(data.ele, data.measurement(:,1:4), data.measurement(:,6))
title('Measured apparent $\rho_{a}$ [$\Omega m$]','Interpreter', 'latex')

nexttile
plot_pseudosection(data.ele, data.measurement(:,1:4), rho_a.')
title('First forward calculated $\rho_{a}$  [$\Omega m$]','Interpreter', 'latex')
set(gcf,'color','w');
% saveas(ffsct,sprintf('norms/sections.png'));



%% taylor test

% resort of sensitivity matrix

for i = 1:numel(x_a)%171

    
    S_all(i,:) = (S{i});

end

% set up handle for dc problem
fwp_dc_problem_handle = @(param) fwp_dc_problem(type, x, y, X, Y, x_a, x_b, x_m, x_n, y_abmn, I_A, bnd_idx_Neum, bnd_idx_Dch, k, param);

% taylor test

% taylor_test_default(fwp_dc_problem_handle, param, S_all);
Js2Jm = @(J,s) J * diag(s);
s2m = @(s) log(s);
m2s = @(m) exp(m);
ftt = figure;
taylor_test_default(fwp_dc_problem_handle, param, S_all,Js2Jm,s2m,m2s);
% saveas(ftt,sprintf('norms/first_taylor.png'));


%% Lcurve trial


% Set the list of lambdas 
lambdas_list = 100000:10000:10^6;


% Set the reference model and approach
%% Levenberg
WeightType = 'Levenberg'; % weighting strategy for Lcurve and further computations
% ref_model = P(:); %homogeneous
ref_model = zeros(size(param)); % zero

% % try L-curve approach
% % use optimal lambda as initial value for inversion 
[~,~,~,~,lambda_opt_list] = Lcurve(lambdas_list,S_all,rho_a.', WeightType,x,y,ref_model );

%% Tikhonov, C - 'Gradient' 
% works with ref_model = P(:) or ref_model = zeros(size(param))
% and parameters = P(:). Should be 23 iterations

% WeightType = 'Gradient'; % weighting strategy for Lcurve and further computations
% ref_model = P(:); % homogeneous
% ref_model = zeros(size(param)); $ zero model

% L-curve approach
% use optimal lambda as initial value for inversion 
% doesn't help in this case
% [~,~,~,~,lambda_opt_list_GR] = Lcurve(lambdas_list,S_all,rho_a.', WeightType,x,y,ref_model );
%% Tikhonov, C - 'Laplace'
% WeightType = 'Laplace'; % weighting strategy for Lcurve and further computations
% ref_model = P(:); %causes Nans
% ref_model = zeros(size(param));% causes NaNs in S
% try L-curve approach
% use optimal lambda as initial value for inversion 
% [~,~,~,~,lambda_opt_list_LA] = Lcurve(lambdas_list,S_all,rho_a.', WeightType,x,y,ref_model );


%% Discrepancy approach

% try  discrepancy approach
error_level = 1e-5; 
d_syn = rho_a.' + ((data.rel_err * mean(rho_measured))/100) * (0.99+0.02*rand) ;
decrease_factor = 0.5;% reduces lambda
lambda_max = 10^12; % maximal lambda
[lambda_err,err] = lambdas_error_level(S_all,d_syn,error_level,x,y,WeightType,ref_model,decrease_factor,lambda_max);

% plot results
figure()
% plot oprimal lambda with the lowest error level
semilogy(lambda_err,err,'-o',lambda_err(find(min(err) == err,1)),min(err(find(min(err) == err,1))),'*r')
xlabel('$\lambda$','Interpreter','latex')
ylabel('Error','Interpreter','latex')
title("Error vs Lambda")
legend('error','optimal value')
% define optimal lambda according the discrepancy approach
indl = find(min(err) == err,1);
lambda_opt_discrepancy = lambda_err(indl);
%% Set optimal lambda value
l_type ='Constant';
% l_type = 'Lcurve' ;
% l_type = 'Discrepancy';
switch l_type
    case 'Constant'
        lambda_opt = 120000;% best option  for all 
    case 'Lcurve'
        lambda_opt = lambda_opt_list;  
    case 'Discrepancy'
        lambda_opt = lambda_opt_discrepancy;
    otherwise
        lambda_opt = 120000;
end



%% Gauss-Newton process
clear d
clear d_norm
clear param_num

% max iterations
k_max = 50; %50;%15

% set lambda
lambda0 = lambda_opt;

% set parameters
param_init = param; % initial model parameters
d_init = rho_a.'; % initial observation data

% set up matrices which will be further filled later
param_num(:,1) = log(param_init);
d(:,1) = d_init;

% set handle for sensitivity matrix
sens_handle = @(param) sensitivity_mat(type, x, y, X, Y, x_a, x_b, x_m, x_n, ...
              y_abmn, I_A, bnd_idx_Neum, bnd_idx_Dch, k, T, param);

% for line search
a_min = 0.07;%1*10^(-3);


b = 0.4;% for cooling approach

for i = 1:k_max
    lambda = lambda0 * 0.5;
%   lambda = lambda0 * exp(-b*(i-1));
    delta_d(:,i) = rho_measured-d(:,i);
    
 
    
    if i == 1
        S_inv{i} = S_all;
    else
        S_inv{i} = sens_handle(exp(param_num(:,i)));
    end

    S_num = S_inv{i}*diag(exp(param_num(:,i)));% transformed S

    
    [C,dp_num(:,i+1)] = WeightedSolveReg(x,y, S_num,(rho_measured-d(:,i)),lambda,WeightType,ref_model);% regression
    
    alpha  = lineSearch_DC(fwp_dc_problem_handle, exp(param_num(:,i)), dp_num(:,i+1), rho_measured, a_min);%Line Search

    param_num(:,i+1) = param_num(:,i) + alpha * dp_num(:,i+1); % model update

    d(:,i+1) = fwp_dc_problem_handle(exp(param_num(:,i+1))); % calculate new data  

%   d_norm(:,i) = norm(rho_measured-d(:,i+1)).^2;
    d_norm(:,i) = norm(rho_measured-d(:,i+1));
    
    % termination condition 1
    % if the ratio of current norms and measured data less than relative error, abort the GN
     if (d_norm(:,i)/norm(rho_measured,2)) < data.rel_err/100 
        sprintf('Termination condition is reached at %1.d iteration',i)
        k_max = i;
        break
     end
    % termination condition 2: slower 
    % if norms grow, abort the GN
%      if (d_norm(:,i)) > d_norm(:,i-1) 
%         sprintf('Termination condition is reached at %1.d iteration',i)
%         k_max = i;
%         break
%     end
    
    
end

param_res = exp(param_num);
%% taylor test
 gt = figure;

for i=1:k_max
    taylor_test_default(fwp_dc_problem_handle, exp(param_num(:,i)), S_inv{i}, Js2Jm,s2m,m2s);
    pause(0.5)
    
    % uncomment and edit name to save desired plot
    % saveas(gt,sprintf('taylor_test_iters/dampedgrad_test_%d.png',i));
end
%% Norm behaviour
figure()
plot(d_norm,'o')
ylabel('$\|$ $\Delta$d $\|$','Interpreter','latex');
xlabel('Iterations');
title('Norms trend')

%% plot pseudosection
% figure(8);
hh = figure;
for i= 1:k_max

plot_pseudosection(data.ele, data.measurement(:,1:4), d(:,i+1))
pause(0.5)
title(['inverse data $\rho_{a}$ [$\Omega m$]' num2str(i)],'Interpreter', 'latex')
set(gcf,'color','w');
% saveas(hh,sprintf('inverse_iters/bad_alpha_%d.png',i));
% saveas(hh,sprintf('inverse_iters/damped_grad%d.png',i));
end

%% Remark
% The following sections are intended for the implementation of algorithms
% having apriori data, i.e reference model.
% For example, one can take any output model vector of any iteration 
% to explore the performance of the regularized regression in general case.


%% with reference model from
%Doesn't work
ref_model = param_res(:,round(k_max/10));
% Tikhonov, C - 'Gradient'
% WeightType = 'Gradient'; % weighting strategy for Lcurve and further computations

% Tikhonov, C - 'Laplace'
 WeightType = 'Laplace'; % weighting strategy for Lcurve and further computations

%% Gauss-Newton
k_max_C = 40 ;%50;%15

% lambda0_C = 12e4;% set higher value
lambda0_C = lambda_opt; 
param_init_C = param_res(:,round(k_max/10)); % initial model parameters
d_init_C = rho_a.'; % initial observation data




% set up matrices which will be further filled later
param_num_C(:,1) = log(param_init_C);
d_C(:,1) = d_init_C;


% for line search
a_min_C = 0.07;%1*10^(-3);

% for cooling approach
b_C = 0.2;

for i = 1:k_max_C
    lambda_C = lambda0_C * 0.5;
%     lambda = lambda0_C * exp(-b*(i-1));
    delta_d_C(:,i) = rho_measured-d_C(:,i);
    
 
    
    if i == 1
        S_inv_C{i} = S_all;
    else
        S_inv_C{i} = sens_handle(exp(param_num_C(:,i)));
    end

    S_num_C = S_inv_C{i}*diag(exp(param_num_C(:,i)));

    
    [C,dp_num_C(:,i+1)] = WeightedSolveReg(x,y, S_num_C,(rho_measured-d_C(:,i)),lambda_C,WeightType,ref_model);
    alpha_C  = lineSearch_DC(fwp_dc_problem_handle, exp(param_num_C(:,i)), dp_num_C(:,i+1), rho_measured, a_min_C);

    param_num_C(:,i+1) = param_num_C(:,i) + alpha_C * dp_num_C(:,i+1); % model update

    d_C(:,i+1) = fwp_dc_problem_handle(exp(param_num_C(:,i+1))); % calculate new data  

%     d_norm(:,i) = norm(rho_measured-d(:,i+1)).^2;
    d_norm_C(:,i) = norm(rho_measured-d_C(:,i+1));

    % termination condition

    if (d_norm_C(:,i)/rho_measured) < data.rel_err/100 
        sprintf('Termination condition is reached at %1.d iteration',i)
        k_max_C=i;
        break
    end
    
end

param_res_C = exp(param_num_C);



%% Norm behaviour
figure()
plot(d_norm_C,'o')
ylabel('$\|$ $\Delta$d $\|$','Interpreter','latex');
xlabel('Iterations');
title('Norms trend')
%% plot pseudosection
% figure(8);
hc = figure;
for i= 1:k_max_C

plot_pseudosection(data.ele, data.measurement(:,1:4), d_C(:,i+1))
pause(0.5)
title(['inverse data $\rho_{a}$ [$\Omega m$]' num2str(i)],'Interpreter', 'latex')
set(gcf,'color','w');
%  saveas(hc,sprintf('inverse_iters/damped_grad%d.png',i));
end
%% taylor test
 gtc = figure;

for i=1:k_max_C
    taylor_test_default(fwp_dc_problem_handle, exp(param_num(:,i)), S_inv{i}, Js2Jm,s2m,m2s);
    pause(0.5)
%     saveas(gtc,sprintf('taylor_test_iters/grad_test_%d.png',i));
end



