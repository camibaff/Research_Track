function CC_ = plt_crosscorr(filename,N,T,WIN,bin_num,xdata,num_ref)
% Plot the Cross-Correlograms between a reference neuron and N-target neurons
% Input:
%       filename : cell array containing cells' file names (i.e. cell0, cell1, cell2, ...)
%       N : number of cells
%       T : duration of recording
%       WIN : window size for Cross-Correlogram
%       bin_num : number of bins
%       xdata : spike time differences
%       num_ref : index to select the reference neuron in the cell array of N filenames (i.e. cell4 is denoted as filename{5})
% Output:
%       CC_ : cell array (1 x N-1) containing the Cross-Correlograms and the informations of each comparison between two neurons
% NOTE: 
% Each cell CC_{i} : contains a cell array which is the output of the function linear_crossCorrelogram applied between two neurons.
% i.e. Between the reference neuron and the i-th cell (CC_{i}{1}: list of differential spike times, CC_{i}{2}: list of histogram values, CC_{i}{3}: number of spikes of the reference neuron, CC_{i}{4}: number of spikes of the ith neuron).
figure
sgtitle(['Cross-Correlograms for reference neuron: cell',num2str(num_ref-1)]) 

for i=1:N
    
    if i~=num_ref

    cc_list = linear_crossCorrelogram(filename{num_ref},filename{i},T,WIN,bin_num);
    CC_{i} = cc_list;
    
    subplot(4,5,i)
    bar(xdata,cc_list{2})
    
    title(['From cell',num2str(num_ref-1),' to cell',num2str(i-1)]);    
    
    end
end

end
