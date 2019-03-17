% This function sorts the trials of the dexterit table. every second trial 
% ist the movement backwords and must be removed.

function [output_data] = f_sort_trials(input_data)

sorted_trials = [];

for ii = 1: size(input_data,2)
    vector_trials = [];
    for i = 2:size(input_data,1)
        block = str2num(cell2mat(input_data(i,ii)));  
        del_idx = 2:2:length(block);
        block(del_idx) = [];
        vector_trials(end+1:end+length(block)) = block;
    end
    sorted_trials(ii,1:length(vector_trials)) = vector_trials;    
end

output_data = sorted_trials';

end