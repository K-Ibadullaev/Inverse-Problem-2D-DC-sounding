function alpha  = lineSearch_DC(fwp_handle, param, dp_num, d, a_min)

         % fwp_handle ... handle for direct problem
         % param ... model parameters (vector)
         % dp_num ... model update
         % d ... measured data 
         % a_min ... lowest alpha

%          a = 1;
            a = (1+5^0.5)/2;
         phi_0 = (norm(d - fwp_handle(param)));
         
         p_update = param + a * dp_num;

         phi_check = (norm(d - fwp_handle(p_update)));
             
         while (phi_check > phi_0) && (a > a_min)

                     a = a/2;
                     
                     p_update = param + a * dp_num;

                     phi_check = (norm(d - fwp_handle(p_update)));
         end

         alpha = a;

end