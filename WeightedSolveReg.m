function [C,x ] = WeightedSolveReg(x,y, S,d,lambda,WeightType,ref_model)
% ref_model is reference model parameters. For Levenberg method it's zero
% vector
% S - sensitivity S
% d is a data d; for nonlinear it's delta d
% x, y mesh nodes for assembling weight 
% WeightType - select an approach for C:
% 'Laplace', 'Gradient' or 'Levenberg'

% assemble weight matrix

switch WeightType
    case 'Laplace'
      C = assemble_parameter_derivative(x,y,'Laplace'); 

    case 'Gradient'
      C = assemble_parameter_derivative(x,y,'Gradient');

    case 'Levenberg'
      C = eye(size(S.'*S));
      
    otherwise
        disp('Wrong entry for approach')
      
end

x =(S.'*S + lambda * C.'*C)\(S.'*d+lambda * C.'*C*ref_model(:)) ;
end

