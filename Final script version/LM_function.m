function solution_LM = LM_function (NPAR,WIN,bin_num,hist_array,delay,tau,reg)
% LM function
% Input:
%       NPAR: total number of parameters
%       WIN: window size of the Cross-Correlogram
%       bin_num: bin width of the histogram
%       hist_array: list of histogram values of the Cross-Correlogram
%       delay: synaptic transmission delay
%       tau: typical time scale of synaptic impact
%       reg: regularization constant for the fitting
% Output: 
%       solution_LM: coefficients that best fit the nonlinear function myfunction to the Cross-Correlogram

% The MAP inference for theta={J12,J21,a(t)} is performed using the LM method. 
% Using the Matlab function lsqcurvefit to fit the function myfunction.

x0 = zeros(1,NPAR);     % initial parameters
x0(end-1:end) = 0.1;    % set the last two parameters at 0.1

xdata = linspace(-WIN,WIN,bin_num);     % input data for the model
ydata = hist_array;                     % histogram values to be fitted

options = optimset('Algorithm', 'levenberg-marquardt');
solution_LM = lsqcurvefit(@(p,x)myfunction(p,x,NPAR,delay,tau,reg),x0,xdata,ydata,[],[],options);

end