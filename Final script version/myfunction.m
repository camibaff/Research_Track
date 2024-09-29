function CC_ = myfunction(par,x,NPAR,delay,tau,reg)
% Function to be fitted by the LM_function
% Input:
%       par : array function handle @(p,x)myfunction
%       x : array function handle @(p,x)myfunction
%       NPAR : total number of parameters
%       delay : synaptic transmission delay
%       tau : typical time scale of synaptic impact
%       reg : regularization constant - is an optional input to adjust the fitting
% Output:
%       CC_ : estimated parameters of the Cross-Correlogram

% Function: c(t) = exp(a(t) + Jij*f(t) + Jji*f(-t))



% for ii = 1:reg-1
%     par((ii+1):reg:(NPAR-reg)) = par(ii:reg:NPAR-(reg+1));
% end % ii

CC_ = par(1:NPAR-2) + par(NPAR-1)*func_f(x,delay,tau(1)) + par(NPAR)*func_f(-x,delay,tau(2));
CC_ = exp(CC_);

end
