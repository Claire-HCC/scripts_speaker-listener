clear all;
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';

binSize_tested=[1 5 10 20 30 40]; % tr;
lags_tested={-10:10 -10:-4, -10:-3, -10:-1};
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(table2array(roi_table(:,1)));

ri=11; rname=rnames{ri};
for ei=3:4;
    exp=experiments{ei};
    
    for lagi=3;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for binSizei=3;%1:length(binSize_tested)
            binSize=binSize_tested(binSizei);
            
            load(['Y:\claire\speaker-listener\' exp '\sound\listenerEvents.mat']);
            onsets=[1 ; find(diff(listenerEvents_trVector)>=1)+1 ];
            % each time as the center of the sliding bin;
          
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'herd','herd_null','pfdr','t','p','rnames');
            herd_z(:,1)=atanh(herd);
            herd_null_z=atanh(herd_null);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
            r2_sl=r2;
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL.mat'],'r2');
            r2_sl_perm=r2;
            r2_sl_perm_m=nanmean(r2,3);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
            r2_ll=r2;
            r2_ll_m=nanmean(r2,3);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/perm/binSize' num2str(binSize) '_lag0-0_permSL.mat'],'r2');
            r2_ll_perm=r2;
            r2_ll_perm_m=nanmean(r2,3);
            
            t=nan;
            figure;
            for ti=1:size(r2_sl,2);
                [~,~,~,stats]=ttest(squeeze(r2_sl_perm(ri,ti,:)),squeeze(r2_sl(ri,ti)));
                t(ti)=-stats.tstat;
            end
            area((listenerEvents_trVector==0)*max(t),'facecolor',[0.5 0.5 0.5],'linestyle','none');
            hold on
            plot(t,'r','linewidth',2);
            hold off;
            
            ylabel({'real vs. null SL coupling','t-value'});
            xlabel('time (TR)')
            title({exp,rname,['binSize' num2str(binSize) ', lag' num2str(min(lags)) '-' num2str(max(lags))]});
            xlim([0 size(r2_sl,2)])
            set(gca,'fontsize',14)
            
            erp_win=-5:10;
            % remove short event
            
            figure;
            erp=nan([length(onsets), length(erp_win)]);
            for evi=1:length(onsets);
                if sum(ismember(onsets(evi)+erp_win,1:size(r2_sl,2)))==length(erp_win);
                    erp(evi,:)=r2_sl(11,onsets(evi)+erp_win);
                end
            end
            plot(erp_win,nanmean(erp),'r','linewidth',2);
            hold on
            
               erp_perm=nan([length(onsets), length(erp_win), size(r2_sl_perm,3)]);
            for evi=1:length(onsets);
                if sum(ismember(onsets(evi)+erp_win,1:size(r2_sl,2)))==length(erp_win);
                erp_perm(evi,:,:)=r2_sl_perm(11,onsets(evi)+erp_win,:);
            end
            end
            quantilePlot(squeeze(nanmean(erp_perm,1))',erp_win,'k',0.3,[0 1])
            legend('real','null');
            legend boxoff
            
            
            %             erp_m=nanmean(erp);
            %             erp_perm_m=squeeze(nanmean(erp_perm))';
            %             t=nan;
            %             for ti=1:length(erp_win);
            %                 [~,~,~,stats]=ttest(erp_perm_m(:,ti),erp_m(ti));
            %                 t(ti)=-stats.tstat;
            %             end
            %
            %             plot(erp_win,t,'k','linewidth',2);
            %             hold off
            grid on
            title({exp,rname,['binSize' num2str(binSize) ', lag' num2str(min(lags)) '-' num2str(max(lags))]});
            %  ylabel({'real vs. null SL coupling','t-value'});
            ylabel({'Speaker-Listener coupling','{R-squared}'});
            xlabel('time to event onset (TR)')
            set(gca,'fontsize',14)
        end
    end
end
