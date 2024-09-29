function time_profile = func_f(sec,delay,tau)
% The time profile of the synaptic interaction
% Input:
%       sec: time from the spikes of the reference neuron
%       delay: synaptic transmission delay
%       tau: typical time scale of synaptic impact
% Output:
%       time_profile: time profile of the monosynaptic interaction

% The monosynaptic interaction of the timescale tau is modelled with a fixed time profile.

time_profile = zeros(1,length(sec));
for i=1:length(sec)

    if sec(i) >= delay
    time_profile(i) = exp(-(sec(i)-delay)/tau);
    
    else    
    time_profile(i) = 0;
    
    end
end

end