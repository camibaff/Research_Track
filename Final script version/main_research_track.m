%% Reconstructing neuronal circuitry applying GLMCC method

%% 1. Introduction:
%The aim of this work is to test the method developed by Kobayashi and his team to reconstruct neuronal circuitry from parallel spike trains.
% Starting from parallel spike trains, it allows to infer the connectivity among neurons by applying a Generalized Linear Model (GLM) to spike Cross-Correlations (CC). 
% Due to the simultaneous use of these two techniques, the method is called GLMCC method.
% Their approach does not reduce dimensionality, but uses all the data to carry out mesoscopic neuroanatomy, that is, to reveal the fine neuronal circuitry in which neural circuit computation is carried out.
% From such an high channel count recordings, they are able to estimate neuronal connectivity by quantifying the degree to which firing from a given neuron is influenced by the firing of neurons from which the index neuron is receiving input.

%% 2. Method implementation
%From their Python version of the GLMCC, we implemented in Matlab a different version of the method, following the main steps of the original one.
% In particular, we reduced the steps and the number of function necessary to obtain the final result, taking advantage of some functions already implemented in Matlab.
% In this phase we have used two recordings from the simulation data to test the method: "cell4sim" and "cell9sim".

clear all
close all
clc

% 2.1 Constraints definition
WIN = 50;       % window size for cross-correlogram [ms]
DELTA = 1;      % bin width for histogram [ms]--> step size for excitatory neurons
MAX = 200000;
T = 5400;       % duration of recording [s]
tau = [4,4];    % typical time scale of synaptic impact
gamma = 0.0005;
beta = 2/gamma;
delay = 1;      % transmission delay
reg = 3;        % regulation constant to plot the fit model
bin_width = DELTA;
bin_num = 2* WIN / bin_width + 1;   % number of bins
NPAR = bin_num + 2;                 % 103 total parameters
xdata = linspace(-WIN,WIN,bin_num);

% 2.2 Plot Cross-Correlogram
% To estimate neuronal connectivity between each pair of neurons the first step is to construct the Cross-Correlogram.
% Cross-correlation (CC): obtained by collecting spike times of a neuron measured relative to every spike of a reference neuron.
% In order to do that we have created the function linear_crossCorrelogram to obtain the list of differential spike times and the list of histogram values. This function also returns the number of spikes of both cells.
cc_list = linear_crossCorrelogram('cell9sim.txt','cell4sim.txt',T,WIN,bin_num);
figure
bar(xdata,cc_list{2})
xlabel('Spike Time Differences [ms]');
ylabel('Count');
title('Cross-Correlogram');

% 2.3. Estimate the GLMCC parameters
% The second step is to compute the parameters of the GLMCC.
% The GLMCC equation can be derived from a two-body Generalized Linear Model (GLM) describing individual neurons interacting with each other, given by the firing rates of two neurons.
% Putting some constraints on the fluctuations and on the temporal profile of the monosynaptic impact, they obtained the GLMCC equation.
% The maximum a posteriori (MAP) inference for the set of parameters of the GLMCC is performed using the Levenberg-Marquardt (LM) method, which interpolates between the Newton method and the gradient descent.
% In order to do that we have created the function LM_function, which applies the Matlab function lsqcurvefit, that is a nonlinear least-squares solver to minimize the function.
solution_LM = LM_function (NPAR,WIN,bin_num,cc_list{2},delay,tau,reg);

% 2.4 Fitting Cross-correlogram with GLMCC
% In the end the last step is to fit the CC with the estimated parameters.
% The figure represents both, all the contributions of the GLMCC equation, and the final result of the parameter estimation. 
f=figure;
bar(xdata,cc_list{2},'FaceAlpha',0.4,'EdgeAlpha',0.1)
[t,s] = title('Cross-Correlogram fitting GLMCC','From cell9sim to cell4sim');

hold on

plot(xdata,exp(solution_LM(1:NPAR-2)),'LineWidth',1.5)  % a(t)
plot(xdata(1:end-1),exp(solution_LM(NPAR-1)*func_f(xdata(1:end-1),delay,tau(1))),'LineWidth',1.5)   % J12
plot(xdata(1:end-1),exp(solution_LM(NPAR)*func_f(-xdata(1:end-1),delay,tau(2))),'LineWidth',1.5)    % J21

plot(xdata,myfunction(solution_LM,xdata,NPAR, delay, tau,reg),'LineWidth',1.5)  % final result
legend(f.Children.Children(end-1:-1:1),{'a(t)','J12','J21','total'})
hold off

% Representation without a(t):
f1=figure;
bar(xdata,cc_list{2},'FaceAlpha',0.4,'LineStyle','none')
[t1,s1] = title('Cross-Correlogram fitting GLMCC','From cell9sim to cell4sim');

hold on

plot(xdata(1:end-1),exp(solution_LM(NPAR-1)*func_f(xdata(1:end-1),delay,tau(1))),'LineWidth',1.5)   % J12
plot(xdata(1:end-1),exp(solution_LM(NPAR)*func_f(-xdata(1:end-1),delay,tau(2))),'LineWidth',1.5)    % J21

plot(xdata,myfunction(solution_LM,xdata,NPAR, delay, tau,reg),'LineWidth',1.5)  % final result

legend(f1.Children.Children(end-1:-1:1),{'J12','J21','total'})
hold off

%% 3. Testing experimental data
% In this section we have tested the implemented method with experimental data, for a total of 20 spike trains recordings available.
% First, we have plotted all the Cross-Correlograms, concerning the comparison between each couple of neurons.
% Then, looking at the obtained plots, we have chosen the most significant Cross-Correlograms to apply the GLMCC method. The results are shown in the 3.2 subsection.
addpath('experimental_data')
experimental_data = natsortfiles(dir('experimental_data'));
experimental_data = experimental_data(4:end);
N = length(experimental_data);

for i=1:N

    filename{i} = experimental_data(i).name;
    
end

% 3.1 Cross-correlograms
% Each figure shows the Cross-Correlograms relative to the comparison of a specific reference neuron with all the other cells.
for j=1:N

CC_j = plt_crosscorr(filename,N,T,WIN,bin_num,xdata,j);
eval(['CC_' num2str(j-1) '=CC_j']);

end

%% 3.2 Fitting Cross-Correlogram with GLMCC
% In each subsection, the GLMCC method is applied to a specific Cross-Correlogram selected between the 19 ones contained in the i-th cell array CC_ (computed in the previous section).
% NOTE (because of the naming of the cells starting from 0): i.e. the CC_3{14} refers to cell3 and cell13.

% cell1 compared to cell15
solution_LM_1_15 = plt_GLMCC_singleCC (CC_1{16},1,15,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell1 compared to cell16
solution_LM_1_16 = plt_GLMCC_singleCC (CC_1{17},1,16,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell3 compared to cell2
solution_LM_3_2 = plt_GLMCC_singleCC (CC_3{3},3,2,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell3 compared to cell13
solution_LM_3_13 = plt_GLMCC_singleCC (CC_3{14},3,13,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell4 compared to cell9
solution_LM_4_9 = plt_GLMCC_singleCC (CC_4{10},4,9,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell4 compared to cell19
solution_LM_4_19 = plt_GLMCC_singleCC (CC_4{20},4,19,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell6 compared to cell11
solution_LM_6_11 = plt_GLMCC_singleCC (CC_6{12},6,11,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell6 compared to cell16
solution_LM_6_16 = plt_GLMCC_singleCC (CC_6{17},6,16,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell8 compared to cell18
solution_LM_8_18 = plt_GLMCC_singleCC (CC_8{19},8,18,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell9 compared to cell13
solution_LM_9_13 = plt_GLMCC_singleCC (CC_9{14},9,13,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell9 compared to cell19
solution_LM_9_19 = plt_GLMCC_singleCC (CC_9{20},9,19,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell11 compared to cell6
solution_LM_11_6 = plt_GLMCC_singleCC (CC_11{7},11,6,xdata,NPAR,WIN,bin_num,delay,tau,reg);

% cell11 compared to cell18
solution_LM_11_18 = plt_GLMCC_singleCC (CC_11{19},11,18,xdata,NPAR,WIN,bin_num,delay,tau,reg);

%% 4. Delivery
% In this last section we tried to reconstruct the result shown in the paper.
% In order to do that, starting from the identification numbers of source neurons and target neurons, we have created another function that includes both the steps of the method.
source = [1,4,4,7,7,7,8,8,8,9,10,10,10,10,12,13,13,13,13,14,14];
target = [14,11,14,8,14,15,11,14,15,2,1,8,11,14,15,3,12,14,15,11,15];
M = length(source);

[cc_list_tot,solution_LM_tot] = plt_CC_with_GLMCC (source,target,M,T,WIN,bin_num,NPAR,xdata,delay,tau,reg);

% The differences between the Cross-Correlograms shown in the paper and those obtained here, are due to the fact that there is a lack of data among those provided by the authors, necessary to produce this type of figures. 
% They only shared a part of the analysis, so we are not able to reconstruct the entire procedures they described. 
% To demonstrate this fact, even trying to retrieve the images from Python code, the outputs are not aligned with the attached results.

%% 5. Conclusions
% This implementation of the method is more conservative than the original one, indeed the estimated values obtained from the curve-fitting process follows exactly the evolution of the Cross-Correlogram, instead of give a smoothed baseline for the large-scale fluctuations.
% However, it allows to detect the synaptic connections that can be employed in the reconstruction of the neuronal connectivity (i.e. with the Hinton diagram, in which connections are proportional to the PSPs).
% Moreover, the advantage is that it does not require to manually compute any posterior probability like in the original code, resulting in a faster implementation.

% Indeed in the Python version, the set of parameters that characterize c(t) in the GLMCC equation, are determined starting from the Bayes theorem.
% In particular, they obtained the posterior distribution of the set of parameters, given the spike data t_k. 
% To detect short-term synaptic impacts hidden in the large-scale
% fluctuations in the Cross-correlogram, they have adapted a(t) to the slow part of the fluctuations, providing a prior distribution that penalizes a large gradient of a(t).
% Finally, to determine the parameters with the maximum a posteriori estimate, they maximized the logarithm of the posterior distribution.

% In the end, looking at their final result, a possible future direction regarding our script may be to develop an adaptation of the regularization provided by "gamma", in order to highlight short-term synaptic impacts more, while reducing the large-scale fluctuations.

% See the Live Script for formulas. 