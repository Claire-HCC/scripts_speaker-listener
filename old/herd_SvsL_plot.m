clear all;
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';

binSize_tested=[1 5 10 20 30 40]; % tr;
lags_tested={-10:10 -10:-4, -10:-3, -10:-1, -20:-3};
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(table2array(roi_table(:,1)));

ri=49; rname=rnames{ri};
for ei=3;
    exp=experiments{ei};
    % mkdir([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/']);
    
    for lagi=[3];%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for binSizei=[ 1 3]%:4;
            binSize=binSize_tested(binSizei);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd_SvsL/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'herd','herd_LL','pfdr','t','p','rnames');
            herd_z=atanh(herd);
            herd_LL_z=atanh(herd_LL);
            ts=t;
            
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
            
            fsize=[9 10 ];
            figure('unit','centimeter','position',[0 0 fsize]);
            hold on;
            scatter(zscore(r2_sl_m(ri, keptT),0,2),zscore(r2_ll_m(ri, keptT),0,2),40,'r','fileld','MarkerFaceAlpha',.2);
            
            for s=1:listenerN;
                coefs = polyfit(zscore(r2_ll(ri, keptT,s),0,2),zscore(r2_ll0(ri, keptT,s),0,2),1);
                Y_pseudo(s,1:length(keptT)) = polyval(coefs,zscore(r2_sl_m(ri, keptT),0,2));
                
                coefs = polyfit(zscore(r2_sl(ri, keptT,s),0,2),zscore(r2_ll0(ri, keptT,s),0,2),1);
                Y(s,1:length(keptT)) = polyval(coefs,zscore(r2_sl_m(ri, keptT),0,2));
                
            end
            [x,I]=sort(zscore(r2_sl_m(ri, keptT),0,2));
            Y_pseudo=Y_pseudo(:,I);
            Y=Y(:,I);
            
            ciplot_claire(Y_pseudo,x,'k',0.3);
            hold on;
            ciplot_claire(Y,x,'r',0.3);
            hold off
            title({[upper(exp(1)) exp(2:end) ', ' rname],sprintf('R=%.2f, t(%d)=%.2f',herd(ri),sum(~isnan(herd_LL_z(1,:)))-1,ts(ri))}); % ['binSize' num2str(binSize),', lag' num2str(min(lags)) '~' num2str(max(lags))]
            xlabel({'Speaker-Listener coupling','(Normalized R-squared)'});
            ylabel({'Listener-Listener coupling','(Normalized R-squared)'});
            set(gca,'fontsize',14);
            h = findobj(gca,'Type','line');
            legend(h,'Real', 'Pseudo','location','northwest');
            legend boxoff
            
            fsize=[30 9];
            figure('unit','centimeter','position',[0 0 fsize]);
            % quantilePlot((squeeze(r2_ll(ri,keptT,:)))',keptT,'r',0.3,[0 1]);
            ciplot_claire((squeeze(r2_ll(ri,keptT,:)))',keptT,'r',0.3);
            h = findobj(gca,'Type','line');
            set(h,'linewidth',1.5,'linestyle',':')
            hold on
            ciplot_claire((squeeze(r2_sl(ri,keptT,:)))',keptT,'r',0.3);
            hold on
            ciplot_claire((squeeze(r2_ll0(ri,keptT,:)))',keptT,'b',0.3);
            hold off
            xlim([0 max(keptT)])
            h = findobj(gca,'Type','line');
            legend(h,'Listener-Listener Coupling','Real Speaker-Listener Couplung', 'Pseudo-Speaker-Listener Couplung','location','northwest');
            legend boxoff
            title({[upper(exp(1)) exp(2:end) ', ' rname],sprintf('R=%.2f, t(%d)=%.2f',herd(ri),sum(~isnan(herd_LL_z(1,:)))-1,ts(ri))}); % ['binSize' num2str(binSize),', lag' num2str(min(lags)) '~' num2str(max(lags))]
            ylabel('R-squared');
            xlabel('Time (TR)')
            hold off;
            set(gca,'fontsize',14);

        end
    end
end

