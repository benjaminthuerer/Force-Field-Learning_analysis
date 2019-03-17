%% plot statistics for each trial so that it is seen if an effect is lasting over many trials

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

parameters = {'trialDuration','stat_pred_PD','stat_max_PD', ...
'stat_area_curve','stat_direction_change_raw','stat_first_reaction_raw'};

%% practice * retest
% acc_trials = 6; %mean over trials
for stats = 1:length(parameters)
    stats_param = [];
    for trls = 1:30
        acc_trials = trls;
        load([statisticsPath '\' parameters{stats} '_practice.mat']);
        param_practice = eval(parameters{stats});
        load([statisticsPath '\' parameters{stats} '_transfer.mat']);
        param_savings = eval(parameters{stats});

        t = table(condition_sleep, condition_practice, nanmean(param_practice(:,end-acc_trials:end),2), nanmean(param_savings(:,1:acc_trials),2), ... %first trials practice! change also in the mean
            'VariableNames',{'cond_sleep','cond_practice','last_trials','first_trials'});
        Time = [1 2]';

        rm = fitrm(t,'last_trials-first_trials ~ cond_sleep*cond_practice + cond_sleep + cond_practice','WithinDesign',Time);

        [stat_results,A,C,D] = ranova(rm);
        
        stats_param(:,trls) = stat_results.pValue;
    end
    assignin('base',['pvalues_' parameters{stats}],stats_param);
    disp(['parameter ' parameters{stats} ' done']);
end

for stats = 1:length(parameters)
    yLine = 0.05;
    figure;
    x = eval(['pvalues_' parameters{stats}]);
    plot(x(1,:),'k');
    hold on
    plot(x(2,:),'b');
    plot(x(3,:),'r');
    plot(x(4,:),'g');
    plot(repmat(yLine,1,30),'--r')
    legend('time','sleep*time','practice*time','sleep*practice*time');
    title([parameters{stats}]);
end