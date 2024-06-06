function CC_ = myfunction(par,x,NPAR, delay, tau)
% Function used by the LM.

CC_ = par(1:NPAR-2) + par(NPAR-1)*func_f(x,delay,tau(1)) + par(NPAR)*func_f(-x,delay,tau(2));
CC_ = exp(CC_);
end