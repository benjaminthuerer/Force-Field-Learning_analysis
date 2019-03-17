% plot average over groups (progress over trials)

trials_indices = [7:12, 14:19, 21:26, 28:33, 35:40];

% plot prediction
WB_PD_smooth = smooth(mean(stat_pred_PD(WB,trials_indices),1),10);
% WB_PD_smooth = smooth(WB_PD_smooth,10);
WR_PD_smooth = smooth(mean(stat_pred_PD(WR,trials_indices),1),10);
% WR_PD_smooth = smooth(WR_PD_smooth,10);
SB_PD_smooth = smooth(mean(stat_pred_PD(SB,trials_indices),1),10);
% SB_PD_smooth = smooth(SB_PD_smooth,10);
SR_PD_smooth = smooth(mean(stat_pred_PD(SR,trials_indices),1),10);
% SR_PD_smooth = smooth(SR_PD_smooth,10);

figure;
plot(WB_PD_smooth, 'b')
hold on
plot(WR_PD_smooth,'k')
plot(SB_PD_smooth,'g')
plot(SR_PD_smooth,'r')
legend('WB','WR','SB','SR')

% plot max PD
WB_PD_smooth = smooth(mean(stat_max_PD(WB,trials_indices),1),10);
% WB_PD_smooth = smooth(WB_PD_smooth,10);
WR_PD_smooth = smooth(mean(stat_max_PD(WR,trials_indices),1),10);
% WR_PD_smooth = smooth(WR_PD_smooth,10);
SB_PD_smooth = smooth(mean(stat_max_PD(SB,trials_indices),1),10);
% SB_PD_smooth = smooth(SB_PD_smooth,10);
SR_PD_smooth = smooth(mean(stat_max_PD(SR,trials_indices),1),10);
% SR_PD_smooth = smooth(SR_PD_smooth,10);

figure;
plot(WB_PD_smooth, 'b')
hold on
plot(WR_PD_smooth,'k')
plot(SB_PD_smooth,'g')
plot(SR_PD_smooth,'r')
legend('WB','WR','SB','SR')

% plot max diff
WB_diff_smooth = smooth(mean(stat_max_diff(WB,trials_indices),1),10);
% WB_diff_smooth = smooth(WB_diff_smooth,10);
WR_diff_smooth = smooth(mean(stat_max_diff(WR,trials_indices),1),10);
% WR_diff_smooth = smooth(WR_diff_smooth,10);
SB_diff_smooth = smooth(mean(stat_max_diff(SB,trials_indices),1),10);
% SB_diff_smooth = smooth(SB_diff_smooth,10);
SR_diff_smooth = smooth(mean(stat_max_diff(SR,trials_indices),1),10);
% SR_diff_smooth = smooth(SR_diff_smooth,10);

plot(WB_diff_smooth, 'b')
hold on
plot(WR_diff_smooth,'k')
plot(SB_diff_smooth,'g')
plot(SR_diff_smooth,'r')
legend('WB','WR','SB','SR')