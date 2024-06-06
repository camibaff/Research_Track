function delta_value = K_delta(i,j)
% Kronecker delta
if i == j
    delta_value = 1;
else 
    delta_value = 0;
end

end