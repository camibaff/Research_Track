function solution_LM = plt_GLMCC_singleCC (CC_i_j,i,j,xdata,NPAR,WIN,bin_num,delay,tau,reg)
% Fitting Cross-Correlogram to GLMCC for the i-th cell compared to j-th cell
% Input:
%       CC_i_j : cell array (1 x 4) containing the Cross-Correlogram and the informations of the comparison between the i-th neuron and the j-th neuron
%       i : identification number of the source neuron (i.e. cell_i )
%       j : identification number of the target neuron (i.e. cell_j )
%       xdata : spike time differences
%       NPAR : total number of parameters
%       WIN : window size of the Cross-Correlogram
%       bin_num : bin width of the histogram
%       delay : synaptic transmission delay
%       tau : typical time scale of synaptic impact
%       reg : regularization constant for the fitting
% Output:
%       solution_LM : coefficients that best fit the nonlinear function myfunction to the Cross-Correlogram

% Estimate the parameters
hist_i_j = CC_i_j;      % CC_i_j{2} : list of histogram values
solution_LM = LM_function (NPAR,WIN,bin_num,hist_i_j{2},delay,tau,reg);

% Fitting Cross-correlogram with GLM
f=figure;
bar(xdata,hist_i_j{2},'FaceAlpha',0.4,'EdgeAlpha',0.1)
[t,s] = title('Cross-Correlogram fitting GLMCC',['From cell',num2str(i),' to cell',num2str(j)]);

hold on

plot(xdata,exp(solution_LM(1:NPAR-2)),'LineWidth',1.5)  % a(t)
plot(xdata(1:end-1),exp(solution_LM(NPAR-1)*func_f(xdata(1:end-1),delay,tau(1))),'LineWidth',1.5)   % J12
plot(xdata(1:end-1),exp(solution_LM(NPAR)*func_f(-xdata(1:end-1),delay,tau(2))),'LineWidth',1.5)    % J21

plot(xdata,myfunction(solution_LM,xdata,NPAR, delay, tau,reg),'LineWidth',1.5)  % final result

legend(f.Children.Children(end-1:-1:1),{'a(t)','J12','J21','total'})

end

