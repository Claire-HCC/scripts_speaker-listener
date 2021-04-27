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

ri=2; rname=rnames{ri};
for ei=1:4;
    exp=experiments{ei};
    % mkdir([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/']);
    
    for lagi=3;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for binSizei=1%:4;
            binSize=binSize_tested(binSizei);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'herd','herd_null','pfdr','t','p','rnames');
            herd_z=atanh(herd);
            herd_null_z=atanh(herd_null);
            ts=t;
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
            r2_sl=r2;
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
            r2_ll=r2;
            r2_ll_m=nanmean(r2,3);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL.mat'],'r2');
            r2_sl_perm=r2;
            r2_sl_perm_m=nanmean(r2,3);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/perm/binSize' num2str(binSize) '_lag0-0_permSL.mat'],'r2');
            r2_ll_perm_m=r2;
            
            keptT=find(~isnan(r2_sl(ri,:)) & ~isnan(r2_ll_m(ri,:)));
            fsize=[9 10 ];
            figure('unit','centimeter','position',[0 0 fsize]);
            hold on;
            scatter(zscore(r2_sl(ri, keptT),0,2),zscore(r2_ll_m(ri, keptT),0,2),40,'r','fileld','MarkerFaceAlpha',.2);
            
            for perm=1:size(r2_sl_perm,3);
                coefs(perm,:) = polyfit(zscore(r2_sl_perm(ri, keptT,perm),0,2),zscore(r2_ll_perm_m(ri, keptT,perm),0,2),1);
                Y(perm,1:length(keptT)) = polyval(coefs(perm,:),zscore(r2_sl(ri, keptT),0,2));
            end
            [x,I]=sort(zscore(r2_sl(ri, keptT),0,2));
            Y=Y(:,I);
            ciplot_claire(Y,x,'k',0.3);
            hold on;
            coefs = polyfit(zscore(r2_sl(ri, keptT),0,2),zscore(r2_ll_m(ri, keptT),0,2),1);
            Y = polyval(coefs,zscore(r2_sl(ri, keptT),0,2));
            plot(zscore(r2_sl(ri, keptT),0,2),Y,'r','linewidth',2);
            title({[upper(exp(1)) exp(2:end) ', ' rname],sprintf('R=%.2f, t(%d)=%.2f',herd(ri),sum(~isnan(herd_null_z(1,:)))-1,ts(ri))}); % ['binSize' num2str(binSize),', lag' num2str(min(lags)) '~' num2str(max(lags))]
            xlabel({'Speaker-Listener coupling','(Normalized R-squared)'});
            ylabel({'Listener-Listener coupling','(Normalized R-squared)'});
            set(gca,'fontsize',14);
            h = findobj(gca,'Type','line');
            legend(h,'Real', 'Null','location','northwest');
            legend boxoff
            
           fsize=[30 9];
            figure('unit','centimeter','position',[0 0 fsize]);
           % quantilePlot((squeeze(r2_sl_perm(ri,keptT,:)))',keptT,'r',0.3,[0 1]);
           ciplot_claire((squeeze(r2_sl_perm(ri,keptT,:)))',keptT,'r',0.3);
            h = findobj(gca,'Type','line');
            set(h,'linewidth',1.5,'linestyle',':')
            hold on
            plot(keptT,r2_sl(ri, keptT),'r','linewidth',2);
            plot(keptT,r2_ll_m(ri, keptT),'b','linewidth',2);
            hold off
            xlim([0 max(keptT)])
            h = findobj(gca,'Type','line');
            legend(h,'Listener-Listener Coupling','Real Speaker-Listener Couplung', 'Null Speaker-Listener Couplung','location','northwest');
            legend boxoff
            title({[upper(exp(1)) exp(2:end) ', ' rname],sprintf('R=%.2f, t(%d)=%.2f',herd(ri),sum(~isnan(herd_null_z(1,:)))-1,ts(ri))}); % ['binSize' num2str(binSize),', lag' num2str(min(lags)) '~' num2str(max(lags))]
            ylabel('R-squared');
            xlabel('Time (TR)')
            hold off;
            set(gca,'fontsize',14);
            
           
        end
    end
end

