function cc_list = linear_crossCorrelogram(filename1,filename2,T,WS,bin_num)
% Make Cross-Correlogram (2nd version)
% This version has been optimized with respect to the first one: use of find function instead of implementing the index_linear_search function and use of bsxfun function.
% Input:
%       filename1 : first file containing the spike times (reference neuron)
%       filename2 : second file containing the spike times
%       T : duration of recording [s] - float or int
%       WS : window size for Cross-Correlogram [ms]
%       bin_num : bin width for histogram [ms]
% Output: cell array cc_list containing:
%       cc_list{1} : list of differential spike times
%       cc_list{2} : list of histogram values
%       cc_list{3} : number of cell1's spikes
%       cc_list{4} : number of cell2's spikes

% Open cell data
cell1 = load(filename1);
cell2 = load(filename2);

% Convert spike times to float and filter based on T
cell1 = cell1(cell1 >= 0 & cell1 < T*1000);
cell2 = cell2(cell2 >= 0 & cell2 < T*1000);

fprintf("n_pre (%s): %d  ;  ",erase(filename1,".txt"),length(cell1));
fprintf("n_post (%s): %d\n",erase(filename2,".txt"),length(cell2));

% Make Correlogramm: c_ij (spike time differences within window)
% Compute spike time differences within window WS.
c = [];

for i=1:length(cell2)

    min_val = cell2(i) - WS;
    max_val = cell2(i) + WS;

    % Search for the index of the smallest value in the list greater than the target
    % or the length of the list if no such value exists:
    min_j = find(cell1 > min_val, 1, 'first');
    max_j = find(cell1 <= max_val, 1, 'last');    

    % Select only the differences less than WS
    differences = bsxfun(@minus, cell1(min_j:max_j)', cell2(i));
    valid_differences = differences(abs(differences) < WS);

    c = [c, valid_differences];

end

% Make histogram
% Create the histogram of the differences.
bin_num = bin_num - (mod(bin_num,2) - 1);

hist_array = histcounts(c,linspace(-WS,WS,bin_num+1));

% Remove the central value of the histogram.
hist_array(ceil(bin_num/2)) = 0;

% Result:
cc_list = {c,hist_array,length(cell1),length(cell2)};

end