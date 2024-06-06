function [log_post,log_likelihood] = calc_log_posterior(par,beta,tau,c,n_sp,t_sp,delay_synapse,Gk)
% Calculate log posterior probability

log_likelihood = 0;

for i=1:NPAR
    log_likelihood = log_likelihood + par(i) * c(i);
end

for i=1:NPAR-2
    log_likelihood = log_likelihood - Gk(i);
end

tmp = 0;
for i=1:NPAR-3
    tmp = tmp + (par(i+1)-par(i))^2;
end
tmp = (beta/(2*DELTA)) * tmp;

log_post = log_likelihood - tmp;

end