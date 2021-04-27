
fsize=[20 15];
cols=jet(21);
for ei=1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_betaPeaks' ],'rnames','r','lags','keptT','r_perm','p','lag_pos','lag_neg','lag_both');
        [lag_both_sorted, ris_sorted]=sort(lag_both,'ascend');
        rn=max(find(~isnan(lag_both_sorted)));

        figure('unit','centimeter','position',[0 0 fsize])
        hold on
        for ri=1:rn;
            
            plot(lags,r(ris_sorted(ri),:),'linewidth',2,'color',cols(lag_both_sorted(ri)+11,:));
        end
        hold off
        legend(rnames(ris_sorted(1:rn)),'location','bestoutside');
        legend boxoff
title(exp)
        grid on
    end
end