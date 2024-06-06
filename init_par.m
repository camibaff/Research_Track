function par = init_par(rate,NPAR)
% Initialize parameter

par = ones(NPAR,1);
par = log(rate)*par;
par(NPAR-1) = 0.1;
par(NPAR) = 0.1;

end