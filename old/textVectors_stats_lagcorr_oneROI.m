clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
lags=-10:-3;
binSize=1;

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
exp=experiments{2};
rnames_selected={'HG_R','vV2','pANG_L'};
for ei=3%1:4;
    exp=experiments{ei};
    
    %   mkdir([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/textVectors/']);
    
    load([expdir exp '/sound/text_trVectors_binSize' num2str(binSize)   '.mat'],'text_trVectors');
    textvarNames=text_trVectors.Properties.VariableNames;
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
    r2_sl=r2;
    r2_sl_m=nanmean(r2_sl,3);
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))],'r2','rnames');
    r2_ll=r2;
    r2_ll_m=nanmean(r2,3);
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_minus1Leave1out_bined/binSize' num2str(binSize) '_lag0-0.mat'],'r2');
    r2_ll0=nanmean(r2,4);
    r2_ll0_m=nanmean(r2_ll0,3);
    
    keptT=[min(find(~isnan(r2_sl(1,:,1)))):max(find(~isnan(r2_sl(1,:,1))))];
    [~,tn,listenerN]=size(r2_sl);
    
    ris=find(ismember(rnames,rnames_selected));
    for i=1:length(ris);
        ri=ris(i);
        cols=distinguishable_colors(length(textvarNames));
        fsize=[25 10];
        figure('unit','centimeter','position',[0 0 fsize]);
        subplot(1,2,1);
        hold on;
        r2=r2_sl_m;
        r=lagcorr_claire(table2array(text_trVectors(keptT,:)),repmat(r2(ri,keptT)',1,length(textvarNames)),-15:15);
        title([upper(exp(1)) exp(2:end) ', ' rnames{ri}]);
        for vi=1:length(textvarNames);
            plot(-15:15,r(:,vi),'linewidth',2,'color',cols(vi,:));
        end
        hold off
        grid on
        legend(textvarNames);
        legend boxoff
        title({rnames{ri},'Speaker-Listener'})
        ylim([-0.2 0.2]);
        
        subplot(1,2,2);
        hold on
        r2=r2_ll0_m;
        r=lagcorr_claire(table2array(text_trVectors(keptT,:)),repmat(r2(ri,keptT)',1,length(textvarNames)),-15:15);
        title([upper(exp(1)) exp(2:end) ', ' rnames{ri}]);
        for vi=1:length(textvarNames);
            plot(-15:15,r(:,vi),'linewidth',2,'color',cols(vi,:));
        end
        hold off
        grid on
        legend(textvarNames);
        legend boxoff
        title('Listener-Listener')
        ylim([-0.2 0.2]);
        
    end
end

