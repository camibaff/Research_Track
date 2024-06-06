function result = index_linear_search(list,target,index)
% Search for index with the smallest value in the list which is bigger than target

if index == -1
    index = 1;
    result = index; % the search starts from the beginning

    while result <= length(list) && list(result) <= target
        result = result+1;
    end

else
    result = index; % the search starts from the given index
    while result <= length(list) && list(result) <= target
        result = result+1;
    end

end

end