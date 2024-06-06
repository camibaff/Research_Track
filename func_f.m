function time_profile = func_f(sec,delay,tau)
% The time profile of the synaptic interaction

time_profile = zeros(1,length(sec));
for i=1:length(sec)
    if sec(i) >= delay
    time_profile(i) = exp(-(sec(i)-delay)/tau);
    else
    %vettore di zeri
    time_profile(i) = 0;
    end
end

end