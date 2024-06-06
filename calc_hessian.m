function hessian = calc_hessian(par,beta,tau,delay_synapse,Gk,NPAR)
% Calculate Hessian of log posterior probability --> CONTROLLARE DI NUOVO

hessian = zeros(NPAR,NPAR);
% d^2P/da_kda_l, d^2P/da_kdj
for i=1:NPAR-2 

    for j=1:NPAR 

    %d^2P/da_kdj
    if j==NPAR-1
        x_k = i*DELTA - WIN;

        if x_k > delay_synapse
            
            % if abs(J) < 1.0e-3 , approximate J=0
            if abs(par(NPAR-1)) < 1.0e-3
                hessian(i,j) = tau(1) * exp(par(i)) * func_f(x_k-DELTA,delay_synapse,tau(1)) * (1-exp(-DELTA/tau(1)));
            else
                hessian(i,j) = -tau(1) * exp(par(i))/par(NPAR-1) * exp(par(NPAR-1) * func_f(x_k-DELTA,delay_synapse,tau(1))) - exp(par(NPAR-1) * func_f(x_k,delay_synapse,tau(1)));
            end

        end

    elseif j==NPAR
            x_k = i*DELTA - WIN;

            if x_k <= -delay_synapse

                %if abs(J) < 1.0e-3, approximate J=0
                if abs(par(NPAR)) < 1.0e-3
                    hessian(i,j) = tau(2) * exp(par(i)) * func_f(-x_k,delay_synapse,tau(2)) * (1-exp(-DELTA/tau(2)));
                else
                    hessian(i,j) = -tau(2) * exp(par(i))/par(NPAR) * (exp(par(NPAR) * func_f(-x_k,delay_synapse,tau(2))) - exp(par(NPAR) * func_f(-x_k+DELTA,delay_synapse,tau(2))));
                end
            end

    %d^2p/da_kda_l
    else
        if i==j               % 1° K_delta: l'indice i deve essere uguale al primo? 2° K_delta: indice i deve essere = al terzultimo? 
            hessian(i,j) = -Gk(i) + beta/DELTA * (K_delta(i,1) + K_delta(i,NPAR-2) - 2);
        else                            % (i-2,j-1) e (i,j-1) ???
            hessian(i,j) = beta/DELTA * (K_delta(i-1,j) + K_delta(i+1,j));
        end
    end
    end
end

%d^2P/dJ^2
for i=NPAR-2:NPAR

    for j=1:NPAR

        if j>NPAR-1

            if i==j && j==NPAR-1
                tmp = 0;

                for k=1:NPAR-1
                    x_k = k *DELTA - WIN; %nel codice c'era k+1
                    
                    if x_k > delay_synapse

                        %if abs(J) < 1.0e-3, approximate J=0
                        if abs(par(NPAR-1)) < 1.0e-3
                            tmp = (tau(1)/2) * (func_f(x_k-DELTA,delay_synapse,tau(1))^2) * (1-exp(-2*DELTA/tau(1)));
                            hessian(i,j) = hessian(i,j)-tmp;
                        else
                            tmp = (par(NPAR-1) * func_f(x_k-DELTA,delay_synapse,tau(1))-1) * exp(par(NPAR-1) * func_f(x_k-DELTA,delay_synapse,tau(1)));
                            tmp = tmp - (par(NPAR-1) * func_f(x_k,delay_synapse,tau(1)-1)) * exp(par(NPAR-1) * func_f(x_k,delay_synapse,tau(1)));
                            hessian(i,j) = hessian(i,j) - (tau(1) * exp(par(k)))/(par(NPAR-1)^2) * tmp;
                        end
                    end
                end

            elseif i==j && j==NPAR
                    tmp = 0;

                    for k=1:NPAR-1
                        x_k = k * DELTA - WIN; %nel codice c'era k+1

                        if x_k <= -delay_synapse

                            % if abs(J) < 1.0e-3, aprroximate J=0
                            if abs(par(NPAR,1)) < 1.0e-3
                                tmp = (tau(2)/2) * (func_f(-x_k,delay_synapse,tau(2))^2) * (1-exp(-2*DELTA/tau(2)));
                                hessian(i,j) = hessian(i,j) - tmp;
                            else
                                tmp = (par(NPAR) * func_f(-x_k,delay_synapse,tau(2))-1) * exp(par(NPAR) * func_f(-x_k,delay_synapse,tau(2)));
                                tmp = tmp - (par(NPAR) * func_f(-x_k+DELTA,delay_synapse,tau(2))-1) * exp(par(NPAR) * func_f(-x_k+DELTA,delay_synapse,tau(2)));
                                hessian(i,j) = hessian(i,j) - (tau(2) * exp(par(k)))/(par(NPAR)^2) * tmp;
                            end
                        end
                    end
            end

        else
            hessian(i,j) = hessian (j,i);
        end
    end
end


end