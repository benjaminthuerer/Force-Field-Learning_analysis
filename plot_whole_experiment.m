% plot mean progress over the entire experiment

statistics_parameter = 'stat_max_PD';

szearios = {'practice','savings','transfer'};
path = 'D:\backup_benny\read_c3d_davos\statistics\';

for i = 1:3
    stat_pd = [];
    stat_max_PD = [];
    stat_mean_PD = [];
    rectify_area_curve = [];
    idx_channel = [];
    stat_area_curve = [];
    
    szenario = szenarios{i};
    load([path 'trialPD_' szenario])
    for sbj = 1:size(trialPD,1)
            for trl = 1:size(trialPD,2)
    %             stat_pred_PD120(sbj,trl) = max(trialPD(sbj,trl,2:120),[],3); % prediction in the first 120 ms
    %             stat_pred_PD150(sbj,trl) = max(trialPD(sbj,trl,2:150),[],3); % prediction in the first 150 ms
    %             stat_pred_PD180(sbj,trl) = max(trialPD(sbj,trl,2:180),[],3); % prediction in the first 180 ms
                stat_max_PD(sbj,trl) = max(abs(trialPD(sbj,trl,2:end-1)),[],3); % maximum PD over the whole trial
                stat_mean_PD(sbj,trl) = squeeze(mean(trialPD(sbj,trl,2:end-1),3)); % mean PD over the whole trial
                rectify_area_curve = squeeze(abs(trialPD(sbj,trl,2:end-1)));
                idx_channel = rectify_area_curve <= 0.0045;
                rectify_area_curve(idx_channel) = 0;
                stat_area_curve(sbj,trl) = sum(rectify_area_curve)/meanTrlDuration; % area under the curve
            end
    end
    stat_pd = eval(statistics_parameter);
    assignin('base',[szenario '_pd'],stat_pd);
end

WB_practice_mean = squeeze(mean(practice_pd(WB,:),1));
WR_practice_mean = squeeze(mean(practice_pd(WR,:),1));
SR_practice_mean = squeeze(mean(practice_pd(SR,:),1));
SB_practice_mean = squeeze(mean(practice_pd(SB,:),1));
WB_savings_mean = squeeze(mean(savings_pd(WB,:),1));
WR_savings_mean = squeeze(mean(savings_pd(WR,:),1));
SR_savings_mean = squeeze(mean(savings_pd(SR,:),1));
SB_savings_mean = squeeze(mean(savings_pd(SB,:),1));
WB_transfer_mean = squeeze(mean(transfer_pd(WB,:),1));
WR_transfer_mean = squeeze(mean(transfer_pd(WR,:),1));
SR_transfer_mean = squeeze(mean(transfer_pd(SR,:),1));
SB_transfer_mean = squeeze(mean(transfer_pd(SB,:),1));

figure
subplot(1,3,1)
plot(smooth(smooth(WB_practice_mean)),'r')
hold on
plot(smooth(smooth(WR_practice_mean)),'b')
plot(smooth(smooth(SR_practice_mean)),'g')
plot(smooth(smooth(SB_practice_mean)),'k')
subplot(1,3,2)
plot(smooth(WB_savings_mean),'r')
hold on
plot(smooth(WR_savings_mean),'b')
plot(smooth(SR_savings_mean),'g')
plot(smooth(SB_savings_mean),'k')
subplot(1,3,3)
plot(smooth(WB_transfer_mean),'r')
hold on
plot(smooth(WR_transfer_mean),'b')
plot(smooth(SR_transfer_mean),'g')
plot(smooth(SB_transfer_mean),'k')
legend({'WB','WR','SR','SB'})