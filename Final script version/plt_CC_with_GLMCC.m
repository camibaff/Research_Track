function [cc_list_tot,solution_LM_tot] = plt_CC_with_GLMCC (source,target,M,T,WIN,bin_num,NPAR,xdata,delay,tau,reg)
% Delivery function
% This function allows both to draw the Cross-Correlograms and to fit them with the GLMCC method, starting from 2 lists of reference numbers, for the source neurons and the target neurons. 
% This allows to directly select the required cells for the comparison.
% Input:
%       source : vector containing the identification numbers of the source neurons
%       target : vector containing the identification numbers of the target neurons
%       M : number of comparisons
%       T : duration of recording
%       WIN : window size of the Cross-Correlogram
%       bin_num : bin width of the histogram
%       NPAR : total number of parameters
%       xdata : spike time differences
%       delay : synaptic transmission delay
%       tau : typical time scale of synaptic impact
%       reg : regularization constant for the fitting
% Output:
%       cc_list_tot : cell array (1 x N-1) containing the Cross-Correlograms and the informations of each comparison between two neurons
%       solution_LM_tot : cell array (1 x N-1) containing the coefficients that best fit the Cross-Correlograms - each cell solution_LM_tot{ii} contains the estimated parameters for the ii-th Cross-Correlogram

figure
sgtitle('Delivery: Cross-correlograms fitting GLMCC') 

for ii = 1:M
    s = sprintf('cell%d.txt',source(ii));
    t = sprintf('cell%d.txt',target(ii));
    
    cc_list = linear_crossCorrelogram(s,t,T,WIN,bin_num);
    cc_list_tot{ii} = cc_list;
    
    solution_LM = LM_function (NPAR,WIN,bin_num,cc_list{2},delay,tau,reg);
    solution_LM_tot{ii} = solution_LM;
    
    ax = subplot(4,ceil(M/4),ii);
    hold(ax,"on");
    
    bar(ax,xdata,cc_list{2},'FaceAlpha',0.4,'EdgeAlpha',0.1)
    
    plot(ax,xdata,exp(solution_LM(1:NPAR-2)),'LineWidth',1.5)  % a(t)
    plot(ax,xdata(1:end-1),exp(solution_LM(NPAR-1)*func_f(xdata(1:end-1),delay,tau(1))),'LineWidth',1.5)    % J12
    plot(ax,xdata(1:end-1),exp(solution_LM(NPAR)*func_f(-xdata(1:end-1),delay,tau(2))),'LineWidth',1.5)     % J21
    
    plot(ax,xdata,myfunction(solution_LM,xdata,NPAR, delay, tau,reg),'LineWidth',1.5)   % final result
    title(ax,['from cell',num2str(source(ii)),' to cell',num2str(target(ii))])
end

end
