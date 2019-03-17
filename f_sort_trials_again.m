function [output_data] = f_sort_trials_again(data_struct,vector_trials)

trialName = [];
nameStrings = [];
idx_trials = [];
sort_idx = [];

for i = 1:length(data_struct)
    trialName = data_struct{i}.file_name;
    nameStrings = strfind(trialName,'_');
    idx_trials = vector_trials == str2num(trialName(nameStrings(1)+1:nameStrings(1)+2));
    idx_trials = find(idx_trials);
    sort_idx(i) = idx_trials(str2num(trialName(nameStrings(2)+1:nameStrings(2)+2)));
end

output_data = sort_idx;

end