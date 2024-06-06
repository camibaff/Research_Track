function Gk = calc_Gk(par,beta,tau,c,n_sp,t_sp,delay_synapse)
% Calculate Gk

Gk = zeros(1,NPAR-2); %0 for i=1:NPAR-2
for i=1:NPAR-2
    x_k = i * DELTA - WIN;
    tmp = 0;
    if x_k <= -delay_synapse && abs(par(NPAR) * func_f(-x_k,delay_synapse,tau(2))) > 1.0e-6
        tmp = expint(par(NPAR) * func_f(-x_k,delay_synapse,tau(2)));
        tmp = tmp - expint(par(NPAR) * func_f(-x_k+DELTA,delay_synapse,tau(2)));

        Gk(i) = tmp * exp(par(i) * tau(1));
    elseif x_k > delay_synapse && abs(par(NPAR-1) * func_f(x_k,delay_synapse,tau(1))) > 1.0e-6
        tmp = expint(par(NPAR-1) * func_f(x_k-DELTA,delay_synapse,tau(1)));
        tmp = tmp - expint(par(NPAR-1) * func_f(x_k,delay_synapse,tau(1)));

        Gk(i) = tmp * exp(par(i) * tau(1));
    else
        Gk(i) = DELTA * exp(par(i));
    end
end

end