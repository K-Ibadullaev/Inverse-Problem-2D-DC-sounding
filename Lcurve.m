function [x_lambdas,norm_dAx,norm_Cx,curv_lambda,lambda_opt] = Lcurve(lambdas_list,S,d, WeightType,x,y,ref_model)
for i=1:length(lambdas_list)
%     x_lambdas(:,i)= mySolveReg(S,d,lambdas_list(i), C,ref_model);
    [C,x_lambdas(:, i)]= WeightedSolveReg(x, y, S, d, lambdas_list(i), WeightType, ref_model);
    
%     norm_dAx(i) = norm(d - S * x_lambdas(:,i), 2);
%     norm_Cx(i) = norm(C *( x_lambdas(:,i)- ref_model), 2);

    norm_dAx(i) = norm(d - S * x_lambdas(:,i), 2).^2;
    norm_Cx(i) = norm(C *( x_lambdas(:,i)- ref_model), 2).^2;
end
figure(length(lambdas_list))
subplot(2,1,1)
% plot(norm_dAx,norm_Cx,'o')
semilogy(norm_dAx,norm_Cx,'o')
xlabel([0, 3])


ylabel('$\|$ C(x($\lambda$) - h))$\|$','Interpreter','latex');
title('L - Curve')
subplot(2,1,2)

 curv_lambda = curvature(norm_dAx,norm_Cx,2);
 indl=find(max(curv_lambda) == curv_lambda);
 lambda_opt = max(lambdas_list(indl));
 plot(norm_dAx,curv_lambda,'-o')
% semilogy(norm_dAx,curv_lambda,'-o')
 
 
 ylabel('curv($\lambda$)','Interpreter','latex');
 xlabel('$\|$ Ax($\lambda$) - d $\|$','Interpreter','latex');
 title('Maximum curvature')
 
 
 
end
