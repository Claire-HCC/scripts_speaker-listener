clear all;
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';

binSize_tested=[10 20 30 40]; % tr;
lags_tested={-10:10 -10:-4, -10:-3, -10:-1};
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(table2array(roi_table(:,1)));

ri=10; rname=rnames{ri};
for ei=3;
    exp=experiments{ei};
    % mkdir([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/']);
    for lagi=3;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for binSizei=1:4;
            binSize=binSize_tested(binSizei);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'herd','herd_null','pfdr','t','p','rnames');
            herd(ri)
            herd_z(:,binSizei)=atanh(herd);
            herd_null_z(:,:,binSizei)=atanh(herd_null);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
            r2_sl(:,:,binSizei)=r2;
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
            r2_ll(:,:,:,binSizei)=r2;
            r2_ll_m(:,:,binSizei)=nanmean(r2,3);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL.mat'],'r2');
            r2_sl_perm(:,:,:,binSizei)=r2;
            r2_sl_perm_m(:,:,binSizei)=nanmean(r2,3);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/perm/binSize' num2str(binSize) '_lag0-0_permSL.mat'],'r2');
            r2_ll_perm(:,:,:,binSizei)=r2;
            r2_ll_perm_m(:,:,binSizei)=nanmean(r2,3);
            
            figure; subplot(1,2,1);
            scatter(r2_sl(ri,:,binSizei),r2_ll_m(ri,:,binSizei),40,'k','fileld')
            subplot(1,2,2);
            plot([r2_sl(ri,:,binSizei);r2_ll_m(ri,:,binSizei)]');
            title({[exp ', ' rname],['binSize' num2str(binSize),', lag' num2str(min(lags)) '-' num2str(max(lags))],sprintf('herd=%.2f, null=%.2f',herd_z(ri,binSizei),nanmean(herd_null_z(ri,:,binSizei)))});
            
        end
    end
end

for binSizei=1:4;
    subplot(2,2,binSizei)
    binSize=binSize_tested(binSizei);
    temp=herd_null_z(ri,:,binSizei);
    histogram(temp(:),-1:0.1:1,'normalization','probability');
    hold on;
    line([herd_z(ri,binSizei) herd_z(ri,binSizei)],[0 0.5])
    hold off
    title({[exp ', ' rname],['binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))]});
    grid on;
    ylim([ 0 0.5])
end

keptT=find(~isnan(squeeze(r2_sl_perm(ri,:,1,4))));
figure;
subplot(2,2,1);
quantilePlot((squeeze(r2_sl_perm(ri,keptT,:,1)))',keptT,'k',0.3,[0 1]);
% plot(keptT,squeeze(r2_sl_perm(ri,keptT,:,1))');
hold on
plot(keptT,squeeze(r2_sl(ri,keptT,1)))
grid on
ylim([0 0.8])
title({[exp ' SL, ' rname],['binSze10, ' 'lag' num2str(min(lags)) '-' num2str(max(lags))]})

subplot(2,2,2);
quantilePlot((squeeze(r2_sl_perm(ri,keptT,:,2)))',keptT,'k',0.3,[0 1]);
% plot(keptT,squeeze(r2_sl_perm(ri,keptT,:,1))');
hold on
plot(keptT,squeeze(r2_sl(ri,keptT,2)))
grid on
ylim([0 0.8])
title({[exp ' SL, ' rname],['binSze20, ' 'lag' num2str(min(lags)) '-' num2str(max(lags))]})

subplot(2,2,3);
quantilePlot((squeeze(r2_sl_perm(ri,keptT,:,3)))',keptT,'k',0.3,[0 1]);
% plot(keptT,squeeze(r2_sl_perm(ri,keptT,:,1))');
hold on
plot(keptT,squeeze(r2_sl(ri,keptT,3)))
grid on
ylim([0 0.8])
title({[exp ' SL, ' rname],['binSze30, ' 'lag' num2str(min(lags)) '-' num2str(max(lags))]})

subplot(2,2,4);
quantilePlot((squeeze(r2_sl_perm(ri,keptT,:,4)))',keptT,'k',0.3,[0 1]);
% plot(keptT,squeeze(r2_sl_perm(ri,keptT,:,1))');
hold on
plot(keptT,squeeze(r2_sl(ri,keptT,4)))
grid on
ylim([0 0.8])
title({[exp ' SL, ' rname],['binSze30, ' 'lag' num2str(min(lags)) '-' num2str(max(lags))]})
legend('','permSL','real');
legend boxoff
