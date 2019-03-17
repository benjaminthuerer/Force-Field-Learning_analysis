

min_vel = min(0, min(min([meanTrial_Vel_SB, meanTrial_Vel_SR,meanTrial_Vel_WB,meanTrial_Vel_WR])));
max_vel = max(0, max(max([meanTrial_Vel_SB, meanTrial_Vel_SR,meanTrial_Vel_WB,meanTrial_Vel_WR])));

groups_trial = {meanTrial_SB, meanTrial_SR, meanTrial_WB, meanTrial_WR};
groups_sem = {semTrial_SB, semTrial_SR, semTrial_WB, semTrial_WR};
groups_vel = {meanTrial_Vel_SB, meanTrial_Vel_SR, meanTrial_Vel_WB, meanTrial_Vel_WR};
groups_acc = {meanTrial_Acc_SB, meanTrial_Acc_SR, meanTrial_Acc_WB, meanTrial_Acc_WR};
groups_change = {meanTrial_change_SB, meanTrial_change_SR, meanTrial_change_WB, meanTrial_change_WR};
groups_change_sem = {semTrial_change_SB, semTrial_change_SR, semTrial_change_WB, semTrial_change_WR};
groups_color = {'r', 'g', 'b', 'k'};

%% mean PD plot
figure;
title(['PD error with PD change as color ' szenarios(s)  ' trial ' num2str(min(trials_indices)) '-' num2str(max(trials_indices)) ])
hold on
for iGroup = 1:numel(groups_trial)

    plot(0,0,groups_color{iGroup},'LineWidth',2)

end
legend('SB','SR','WB','WR')
for iGroup = 1:numel(groups_trial)
%for iTrial = 1:numel(folder_struct_selected)

    x = 1:length(groups_trial{iGroup})';
    x = x(1:meanTrlDuration);
    y = groups_trial{iGroup}';
    y = y(1:meanTrlDuration);
    z = zeros(size(x));
    
    %col = linspace(0,1,numel(x));  % This is the color
    %col = [0 groups_vel{iGroup}'];  % This is the color for velocity
    %col = [0 0 groups_acc{iGroup}'];  % This is the color for accelatation
    col = [0 groups_change{iGroup}']; col = col(1:meanTrlDuration);  % This is the color for error(PD) change
    %col = [0 abs(groups_change{iGroup}')]; col = col(1:meanTrlDuration);  % This is the color for error(PD) change rectified values



    surface([x;x],[y;y],[z;z],[col;col], ...
        'facecol','no', ...
        'edgecol','interp', ...
        'linew',8);
    
    plot(x,y,groups_color{iGroup},'LineWidth',2)
    
    
end
hold off

colormap(jet(256))
colorbar

%% errorbar plot
figure;
title(['PD error (+-SEM ' szenarios(s)  ' trial ' num2str(min(trials_indices)) '-' num2str(max(trials_indices)) ])
hold on
for iGroup = 1:numel(groups_trial)
    
    x = 1:length(groups_trial{iGroup})';
    x = x(1:meanTrlDuration);
    y = groups_trial{iGroup}';
    y = y(1:meanTrlDuration);

    sem = groups_sem{iGroup}'; sem = sem(1:meanTrlDuration);  % This is the color for error(PD) change   

    errorbar(x,y,sem,'Color',groups_color{iGroup},'LineWidth',2)

    
end
legend('SB','SR','WB','WR')

%% PD-change (diff of PD progress)
figure;
title(['PD change ' szenarios(s)  ' trial ' num2str(min(trials_indices)) '-' num2str(min(trials_indices)) ])
hold on
for iGroup = 1:numel(groups_trial)
%for iTrial = 1:numel(folder_struct_selected)

    x = 1:length(groups_trial{iGroup})';
    x = x(1:meanTrlDuration);
    
    %col = linspace(0,1,numel(x));  % This is the color
    %col = [0 groups_vel{iGroup}'];  % This is the color for velocity
    %col = [0 0 groups_acc{iGroup}'];  % This is the color for accelatation
    col = [0 groups_change{iGroup}']; col = col(1:meanTrlDuration);  % This is the color for error(PD) change
    %col = [0 abs(groups_change{iGroup}')]; col = col(1:meanTrlDuration);  % This is the color for error(PD) change rectified values

    
    plot(x,col,groups_color{iGroup},'LineWidth',2)
    
    
end
hold off
legend('SB','SR','WB','WR')

%% subtracted PD by the second sample
figure;
title(['PD change subtract offset ' szenarios(s)  ' trial ' num2str(min(trials_indices)) '-' num2str(max(trials_indices)) ])
hold on
for iGroup = 1:numel(groups_trial)
%for iTrial = 1:numel(folder_struct_selected)

    x = 1:length(groups_trial{iGroup})';
    x = x(1:meanTrlDuration);
    
    
    %col = linspace(0,1,numel(x));  % This is the color
    %col = [0 groups_vel{iGroup}'];  % This is the color for velocity
    %col = [0 0 groups_acc{iGroup}'];  % This is the color for accelatation
    col = [0 groups_change{iGroup}']; col = col(1:meanTrlDuration);  % This is the color for error(PD) change
    %col = [0 abs(groups_change{iGroup}')]; col = col(1:meanTrlDuration);  % This is the color for error(PD) change rectified values

    col = col - col(2);
    
    %col = [0 diff(col)];
  
    %col = abs(col);
    
    plot(x,col,groups_color{iGroup},'LineWidth',2)
    
end
hold off
legend('SB','SR','WB','WR')

%% subtracted PD by the second sample +- standarddeviation
figure;
title(['PD change subtract offset (+-SEM) ' szenarios(s)  ' trial ' num2str(min(trials_indices)) '-' num2str(max(trials_indices)) ])
hold on
for iGroup = 1:numel(groups_trial)
%for iTrial = 1:numel(folder_struct_selected)

    x = 1:length(groups_trial{iGroup})';
    x = x(1:meanTrlDuration);
    
    
    %col = linspace(0,1,numel(x));  % This is the color
    %col = [0 groups_vel{iGroup}'];  % This is the color for velocity
    %col = [0 0 groups_acc{iGroup}'];  % This is the color for accelatation
    col = [0 groups_change{iGroup}']; col = col(1:meanTrlDuration);  % This is the color for error(PD) change
    %col = [0 abs(groups_change{iGroup}')]; col = col(1:meanTrlDuration);  % This is the color for error(PD) change rectified values
    col_sem = [0 groups_change_sem{iGroup}']; col_sem = col_sem(1:meanTrlDuration);
    col = col - col(2);
    
    %col = [0 diff(col)];
  
    %col = abs(col);
    
    errorbar(x,col,col_sem,'Color',groups_color{iGroup},'LineWidth',2)
    
end
hold off
legend('SB','SR','WB','WR')

%% absolute of the subtracted PD change by the second sample
figure;
title(['PD change subtract offset rectified values ' szenarios(s)  ' trial ' num2str(min(trials_indices)) '-' num2str(max(trials_indices)) ])
hold on
for iGroup = 1:numel(groups_trial)
%for iTrial = 1:numel(folder_struct_selected)

    x = 1:length(groups_trial{iGroup})';
    x = x(1:meanTrlDuration);
    
    
    %col = linspace(0,1,numel(x));  % This is the color
    %col = [0 groups_vel{iGroup}'];  % This is the color for velocity
    %col = [0 0 groups_acc{iGroup}'];  % This is the color for accelatation
    col = [0 groups_change{iGroup}']; col = col(1:meanTrlDuration);  % This is the color for error(PD) change
    %col = [0 abs(groups_change{iGroup}')]; col = col(1:meanTrlDuration);  % This is the color for error(PD) change rectified values

    col = col - col(2);
    
   % col = [0 diff(col)];
  
    col = abs(col);
    
    plot(x,col,groups_color{iGroup},'LineWidth',2)
    
end
hold off
legend('SB','SR','WB','WR')

%% diff of the subtracted PD change (second derivation of PD)
figure;
title(['PD change subtract offset deviation rectified values (inicated first time of adaptation arround 200 ms) ' szenarios(s)  ' trial ' num2str(min(trials_indices)) '-' num2str(max(trials_indices)) ])
hold on
for iGroup = 1:numel(groups_trial)
%for iTrial = 1:numel(folder_struct_selected)

    x = 1:length(groups_trial{iGroup})';
    x = x(1:meanTrlDuration);
    
    
    %col = linspace(0,1,numel(x));  % This is the color
    %col = [0 groups_vel{iGroup}'];  % This is the color for velocity
    %col = [0 0 groups_acc{iGroup}'];  % This is the color for accelatation
    col = [0 groups_change{iGroup}']; col = col(1:meanTrlDuration);  % This is the color for error(PD) change
    %col = [0 abs(groups_change{iGroup}')]; col = col(1:meanTrlDuration);  % This is the color for error(PD) change rectified values

    col = col - col(2);
    
    col = [0 diff(col)];
  
    col = abs(col);
    
    plot(x,col,groups_color{iGroup},'LineWidth',2)
    
end
xlim([3 500])
hold off
legend('SB','SR','WB','WR')

