function result = linear_crossCorrelogram(filename1,filename2,T,WS,BW)
% Make Cross Correlogram

cell1 = load(filename1);
cell2 = load(filename2);
% Convert spike times to float and filter based on T
cell1 = cell1(cell1 >= 0 & cell1 < T*1000);
cell2 = cell2(cell2 >= 0 & cell2 < T*1000);

fprintf("n_pre: %d\n",length(cell1));
fprintf("n_post: %d\n",length(cell2));
%% Make Correlogramm: c_ij (spike time differences within window)
% Compute spike time differences within window WIN.
c = [];
min_index = 1;
max_index = 1;

for i=1:length(cell2)
    min = cell2(i) - WS;
    max = cell2(i) + WS;

    min_j = index_linear_search(cell1, min, min_index);
    min_index = min_j;

    max_j = index_linear_search(cell1, max, max_index);
    max_index = max_j;

    
    for j=min_j:max_j-1
        if abs(cell1(j) - cell2(i)) < WS
            
            c(i) = cell1(j) - cell2(i);
        end
    end
end
    
%% Make histogram
% Create the histogram of the differences
bin_num = 2*WS / BW; % number of bin

hist_array = histcounts(c,linspace(-WS,WS,bin_num+1));

result = {c,hist_array,length(cell1),length(cell2)};

end