function [lb,err_l] = lambdas_error_level(S,d_syn,error_level,x,y,WeightType,ref_model,decrease_factor,lambda_max)
         % returns error and the optimal lambda within the scope of noise
         % S...sensitivity matrix
         % d_syn...synthetic data, including rel. error
         % lambda_0...initial value for lambda
         % error_level...error level of measured data
         % x,y...mesh nodes in x/y direction, vector
         % type...weighted with identity matrix, first deriviative or
         % second deriviative   
         % decrease_factor - reduces lambda
         % lambda_max - threshold for lambda's
         k = 100;
         lambda_0 = lambda_max;
         
        
         i = 1;   
         while (k>0)
               [~,param_lambda] = WeightedSolveReg(x,y, S,d_syn,lambda_0,WeightType,ref_model);
%                d_cmp = S*param_lambda;
%                err_update = std(d_cmp);

              err_update = sum(abs(S * param_lambda - d_syn)./d_syn)./length(d_syn);
              err_l(i) = err_update;
              lb(i) = lambda_0;
              lambda_0 = round(lambda_0 * decrease_factor);
              k=k-1;
              i=i+1;
              % termination condition
              if (err_update<error_level)
                  i =i-1;
                  fprintf('Stopped at %d iteration',i)
                  break
              end
              
              
         end
    lb = flip(lb);% sort the list in the ascending order
    lb = lb(:,1:i);
          return


end