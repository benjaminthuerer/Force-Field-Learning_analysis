
    mean_group_p = nanmean(param_practice(:,:),2);
    mean_group_s = nanmean(param_savings(:,1:12),2);
    
    [corr_r,corr_p] = corr(mean_group_p,mean_group_s,'type','Spearman')
    corr_r^2