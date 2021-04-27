clear all
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
[networks,~,networksi]=unique(roi_table.network);
lags_tested={-10:10, -40:40, -60:60};
% seed='DMN1';
% ris_seed=find(networksi==find(strcmp(networks,seed)));
ris_seed=11;
binSize=30;

fsize=[45 10];
for ei=1:5;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        for seedi=1:length(ris_seed);
            ri=ris_seed(seedi);
            seed=rnames{ri};
            load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/' seed '_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
            if seedi==1; rz_ll=atanh(r); else
                rz_ll=rz_ll+atanh(r); end
        end
        rz_ll=(rz_ll)/length(ris_seed);
        
        for ni=1:length(networks);
            ris=find(strcmp(networks{ni},roi_table.network));
            rz_ll_network(ni,:)=nanmean(nanmean(rz_ll(ris,:,lags==0,:),4),1);
        end
        figure('unit','centimeter','position',[0 0 fsize]);
        plot(rz_ll_network','linewidth',2);
        title(exp);
        xlabel('Time (TR)');
        ylabel('R(z)');
        grid on
        legend(strrep(networks,'_',' '));
        set(gca,'fontsize',20);
        hold on;
        line([0 size(rz_ll,2)],[0 0],'color','k');
        hold off
        legend boxoff
        xlim([0 700])
    end
    clear rz_ll_network
end




