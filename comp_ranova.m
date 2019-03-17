%% statistics in matlab
statisticsPath = uigetdir('','choose the statistics folder, e.g. C:\...\read_c3d_davos\statistics');
load([statisticsPath '\' 'statistics_header_practice.mat']);

condition_sleep = statistics_header(:,1);
condition_practice = statistics_header(:,2);
condition_all = statistics_header(:,3);

SB = [];
SR = [];
WB = [];
WR = [];

for sbj = 1:length(condition_all)
    if strcmp(condition_all(sbj),'SB')
        SB = [SB; sbj];
    elseif  strcmp(condition_all(sbj),'SR')
         SR = [SR; sbj];
    elseif  strcmp(condition_all(sbj),'WB')
        WB = [WB; sbj];
    elseif  strcmp(condition_all(sbj),'WR')
         WR = [WR; sbj];
    end
end

parameters = {'trialDuration','stat_pred_PD','stat_max_PD','stat_direction_change','stat_first_reaction', ...
'stat_binary_time_elapse','stat_area_curve','stat_direction_change_raw','stat_first_reaction_raw'};

%% practice * retest
acc_trials = 30; %mean over trials
for stats = 1:length(parameters)
    load([statisticsPath '\' parameters{stats} '_practice.mat']);
    param_practice = eval(parameters{stats});
    load([statisticsPath '\' parameters{stats} '_transfer.mat']);
    param_savings = eval(parameters{stats});
    
%     param_practice = param_practice([SR,WR],:);
%     param_savings = param_savings([SR,WR],:);
    
    t = table(condition_sleep, condition_practice, nanmean(param_practice(:,end-acc_trials:end),2), nanmean(param_savings(:,1:acc_trials),2), ... %first trials practice! change also in the mean
        'VariableNames',{'cond_sleep','cond_practice','last_trials','first_trials'});
%     t(strcmp(t.cond_practice,'blocked'),:) = [];
    Time = [1 2]';
    
    rm = fitrm(t,'last_trials-first_trials ~ cond_sleep*cond_practice + cond_sleep + cond_practice','WithinDesign',Time);
%     rm = fitrm(t,'last_trials-first_trials ~ cond_sleep','WithinDesign',Time);
    
    [stat_results,A,C,D] = ranova(rm);
    assignin('base',['anova_' parameters{stats}],stat_results);
    %% compute readaptation using the difference between trials     

%     diff_savings = diff(param_savings,[],2);
%     t_diff = table(condition_sleep, condition_practice, nanmean(diff_savings(:,1:acc_trials),2), ...
%         'VariableNames',{'cond_sleep','cond_practice','readaptation'});
%     
%     rm_diff = fitlm(t_diff,'readaptation ~ cond_sleep*cond_practice + cond_sleep + cond_practice');
%     [stat_results_diff] = anova(rm_diff);
%     assignin('base',['anova_' parameters{stats}],stat_results_diff);
    
    %% export values for mean
    mean_group = [];
    stdErr_group = [];
    mean_group_diff = [];
    stdErr_group_diff = [];
    
    mean_group(1,1) = nanmean(nanmean(param_practice(SB,1:acc_trials),2),1);
    mean_group(2,1) = nanmean(nanmean(param_practice(SR,1:acc_trials),2),1);
    mean_group(3,1) = nanmean(nanmean(param_practice(WB,1:acc_trials),2),1);
    mean_group(4,1) = nanmean(nanmean(param_practice(WR,1:acc_trials),2),1);

    mean_group(1,2) = nanmean(nanmean(param_practice(SB,end-acc_trials:end),2),1);
    mean_group(2,2) = nanmean(nanmean(param_practice(SR,end-acc_trials:end),2),1);
    mean_group(3,2) = nanmean(nanmean(param_practice(WB,end-acc_trials:end),2),1);
    mean_group(4,2) = nanmean(nanmean(param_practice(WR,end-acc_trials:end),2),1);
    
    mean_group(1,3) = nanmean(nanmean(param_savings(SB,1:acc_trials),2),1);
    mean_group(2,3) = nanmean(nanmean(param_savings(SR,1:acc_trials),2),1);
    mean_group(3,3) = nanmean(nanmean(param_savings(WB,1:acc_trials),2),1);
    mean_group(4,3) = nanmean(nanmean(param_savings(WR,1:acc_trials),2),1);

    
    stdErr_group(1,1) = squeeze(std(nanmean(param_practice(SB,1:acc_trials),2),1))/sqrt(numel(SB));
    stdErr_group(2,1) = squeeze(std(nanmean(param_practice(SR,1:acc_trials),2),1))/sqrt(numel(SR));
    stdErr_group(3,1) = squeeze(std(nanmean(param_practice(WB,1:acc_trials),2),1))/sqrt(numel(WB));
    stdErr_group(4,1) = squeeze(std(nanmean(param_practice(WR,1:acc_trials),2),1))/sqrt(numel(WR));

    stdErr_group(1,2) = squeeze(std(nanmean(param_practice(SB,end-acc_trials:end),2),1))/sqrt(numel(SB));
    stdErr_group(2,2) = squeeze(std(nanmean(param_practice(SR,end-acc_trials:end),2),1))/sqrt(numel(SR));
    stdErr_group(3,2) = squeeze(std(nanmean(param_practice(WB,end-acc_trials:end),2),1))/sqrt(numel(WB));
    stdErr_group(4,2) = squeeze(std(nanmean(param_practice(WR,end-acc_trials:end),2),1))/sqrt(numel(WR));
    
    stdErr_group(1,3) = squeeze(std(nanmean(param_savings(SB,1:acc_trials),2),1))/sqrt(numel(SB));
    stdErr_group(2,3) = squeeze(std(nanmean(param_savings(SR,1:acc_trials),2),1))/sqrt(numel(SR));
    stdErr_group(3,3) = squeeze(std(nanmean(param_savings(WB,1:acc_trials),2),1))/sqrt(numel(WB));
    stdErr_group(4,3) = squeeze(std(nanmean(param_savings(WR,1:acc_trials),2),1))/sqrt(numel(WR));

    
    %%
    mean_group_diff(1,1) = nanmean(nanmean(param_savings(SB,1:acc_trials),2)./nanmean(param_practice(SB,end-acc_trials:end),2),1);
    mean_group_diff(2,1) = nanmean(nanmean(param_savings(SR,1:acc_trials),2)./nanmean(param_practice(SR,end-acc_trials:end),2),1);
    mean_group_diff(3,1) = nanmean(nanmean(param_savings(WB,1:acc_trials),2)./nanmean(param_practice(WB,end-acc_trials:end),2),1);
    mean_group_diff(4,1) = nanmean(nanmean(param_savings(WR,1:acc_trials),2)./nanmean(param_practice(WR,end-acc_trials:end),2),1);
    
    stdErr_group_diff(1,1) = nanstd(nanmean(param_savings(SB,1:acc_trials),2)./nanmean(param_practice(SB,end-acc_trials:end),2),1)/sqrt(numel(SB));
    stdErr_group_diff(2,1) = nanstd(nanmean(param_savings(SR,1:acc_trials),2)./nanmean(param_practice(SR,end-acc_trials:end),2),1)/sqrt(numel(SR));
    stdErr_group_diff(3,1) = nanstd(nanmean(param_savings(WB,1:acc_trials),2)./nanmean(param_practice(WB,end-acc_trials:end),2),1)/sqrt(numel(WB));
    stdErr_group_diff(4,1) = nanstd(nanmean(param_savings(WR,1:acc_trials),2)./nanmean(param_practice(WR,end-acc_trials:end),2),1)/sqrt(numel(WR));
    
    groups_color = {'r', 'g', 'b', 'k'};
    x_offset = [1, 1.1, 1.2, 1.3];
    figure;
    figureTitle = parameters{stats};
    titleIdx = strfind(figureTitle,'_');
    figureTitle(titleIdx) = ' ';
    title([figureTitle ' with std Err'])
    hold on
    for iGroup = 1:size(mean_group_diff,1)

        x = x_offset(iGroup);
        y = mean_group_diff(iGroup);
        sem = stdErr_group_diff(iGroup); 

        errorbar(x,y,sem,'Color',groups_color{iGroup},'LineWidth',2)

    end
    legend('SB','SR','WB','WR')
    
    %%
    assignin('base',['mean_' parameters{stats}],mean_group);
    assignin('base',['stdErr_' parameters{stats}],stdErr_group);
    
    %% plot mean with stdErr
    groups_color = {'r', 'g', 'b', 'k'};
    figureTitle = parameters{stats};
    titleIdx = strfind(figureTitle,'_');
    figureTitle(titleIdx) = ' ';
    
    figure;
    title([figureTitle ' with std Err'])
    hold on
    for iGroup = 1:size(mean_group,1)

        x = 1:3;
        y = mean_group(iGroup,:);
        sem = stdErr_group(iGroup,:); 

        errorbar(x,y,sem,'Color',groups_color{iGroup},'LineWidth',2)
    end
    legend('SB','SR','WB','WR')
    
    h = waitforbuttonpress;
    if h==1
        close all
    end
end
    

%% practice * transfer
% 
% acc_trials = 30; %mean over 6 trials (one block)
% for stats = 1:length(parameters)
%     load([statisticsPath '\' parameters{stats} '_practice.mat']);
%     param_practice = eval(parameters{stats});
%     load([statisticsPath '\' parameters{stats} '_transfer.mat']);
%     param_transfer = eval(parameters{stats});
%     
%     
%     t = table(condition_sleep, condition_practice, nanmean(param_practice(:,end-acc_trials:end),2), nanmean(param_transfer(:,1:acc_trials),2), ...
%         'VariableNames',{'cond_sleep','cond_practice','first_trials_prct','first_trials_trnsf'});
%     Time = [1 2]';
%     
%     rm = fitrm(t,'first_trials_prct-first_trials_trnsf ~ cond_sleep*cond_practice + cond_sleep + cond_practice','WithinDesign',Time);
%     [stat_results,A,C,D] = ranova(rm);
%     assignin('base',['anova_' parameters{stats}],stat_results);
%     
%     %% export values for mean
%     mean_group(1,1) = nanmean(nanmean(param_practice(SB,end-acc_trials:end),2),1);
%     mean_group(2,1) = nanmean(nanmean(param_practice(SR,end-acc_trials:end),2),1);
%     mean_group(3,1) = nanmean(nanmean(param_practice(WB,end-acc_trials:end),2),1);
%     mean_group(4,1) = nanmean(nanmean(param_practice(WR,end-acc_trials:end),2),1);
%     mean_group(1,2) = nanmean(nanmean(param_transfer(SB,1:acc_trials),2),1);
%     mean_group(2,2) = nanmean(nanmean(param_transfer(SR,1:acc_trials),2),1);
%     mean_group(3,2) = nanmean(nanmean(param_transfer(WB,1:acc_trials),2),1);
%     mean_group(4,2) = nanmean(nanmean(param_transfer(WR,1:acc_trials),2),1);
%     
%     stdErr_group(1,1) = squeeze(std(nanmean(param_practice(SB,end-acc_trials:end),2),1))/sqrt(numel(SB));
%     stdErr_group(2,1) = squeeze(std(nanmean(param_practice(SR,end-acc_trials:end),2),1))/sqrt(numel(SR));
%     stdErr_group(3,1) = squeeze(std(nanmean(param_practice(WB,end-acc_trials:end),2),1))/sqrt(numel(WB));
%     stdErr_group(4,1) = squeeze(std(nanmean(param_practice(WR,end-acc_trials:end),2),1))/sqrt(numel(WR));
%     stdErr_group(1,2) = squeeze(std(nanmean(param_transfer(SB,1:acc_trials),2),1))/sqrt(numel(SB));
%     stdErr_group(2,2) = squeeze(std(nanmean(param_transfer(SR,1:acc_trials),2),1))/sqrt(numel(SR));
%     stdErr_group(3,2) = squeeze(std(nanmean(param_transfer(WB,1:acc_trials),2),1))/sqrt(numel(WB));
%     stdErr_group(4,2) = squeeze(std(nanmean(param_transfer(WR,1:acc_trials),2),1))/sqrt(numel(WR));
%     
%     mean_group_diff(1,1) = nanmean(nanmean(param_transfer(SB,1:acc_trials),2)./nanmean(param_practice(SB,end-acc_trials:end),2),1);
%     mean_group_diff(2,1) = nanmean(nanmean(param_transfer(SR,1:acc_trials),2)./nanmean(param_practice(SR,end-acc_trials:end),2),1);
%     mean_group_diff(3,1) = nanmean(nanmean(param_transfer(WB,1:acc_trials),2)./nanmean(param_practice(WB,end-acc_trials:end),2),1);
%     mean_group_diff(4,1) = nanmean(nanmean(param_transfer(WR,1:acc_trials),2)./nanmean(param_practice(WR,end-acc_trials:end),2),1);
%     
%     stdErr_group_diff(1,1) = nanstd(nanmean(param_transfer(SB,1:acc_trials),2)./nanmean(param_practice(SB,end-acc_trials:end),2),1)/sqrt(numel(SB));
%     stdErr_group_diff(2,1) = nanstd(nanmean(param_transfer(SR,1:acc_trials),2)./nanmean(param_practice(SR,end-acc_trials:end),2),1)/sqrt(numel(SR));
%     stdErr_group_diff(3,1) = nanstd(nanmean(param_transfer(WB,1:acc_trials),2)./nanmean(param_practice(WB,end-acc_trials:end),2),1)/sqrt(numel(WB));
%     stdErr_group_diff(4,1) = nanstd(nanmean(param_transfer(WR,1:acc_trials),2)./nanmean(param_practice(WR,end-acc_trials:end),2),1)/sqrt(numel(WR));
%     
%     groups_color = {'r', 'g', 'b', 'k'};
%     x_offset = [1, 1.1, 1.2, 1.3];
%     figure;
%     figureTitle = parameters{stats};
%     titleIdx = strfind(figureTitle,'_');
%     figureTitle(titleIdx) = ' ';
%     title([figureTitle ' with std Err'])
%     hold on
%     for iGroup = 1:size(mean_group_diff,1)
% 
%         x = x_offset(iGroup);
%         y = mean_group_diff(iGroup);
%         sem = stdErr_group_diff(iGroup); 
% 
%         errorbar(x,y,sem,'Color',groups_color{iGroup},'LineWidth',2)
% 
%     end
%     legend('SB','SR','WB','WR')
% 
%     assignin('base',['mean_' parameters{stats}],mean_group);
%     assignin('base',['stdErr_' parameters{stats}],stdErr_group);
%     
%     %% plot mean with stdErr
%     groups_color = {'r', 'g', 'b', 'k'};
%     figureTitle = parameters{stats};
%     titleIdx = strfind(figureTitle,'_');
%     figureTitle(titleIdx) = ' ';
%     
%     figure;
%     title([figureTitle ' with std Err'])
%     hold on
%     for iGroup = 1:size(mean_group,1)
% 
%         x = 1:2;
%         y = mean_group(iGroup,:);
%         sem = stdErr_group(iGroup,:); 
% 
%         errorbar(x,y,sem,'Color',groups_color{iGroup},'LineWidth',2)
% 
%     end
%     legend('SB','SR','WB','WR')
%     
%     waitforbuttonpress;
%     
% end



