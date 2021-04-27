

clear all
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');

tic % 15 min
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data

lags_tested={-10:-1,-10:-4};
binSize_tested=[ 10 15 30]; % tr;

for ei=1:4;
    exp=experiments{ei};
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    for binSizei=1:length(binSize_tested);
        binSize=binSize_tested(binSizei);
        
        for lagi=1:length(lags_tested);
            lags=lags_tested{lagi};
            
            load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b_l','b_ls','b_s','F_l','F_ls','F_s2l','F_l2l','F_s','r2_ls','r2_l','r2_s','rnames');
            real=F_s2l;
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rnames{1} '.mat' ],'gdata');
            listenersN=size(gdata,3);
            
            for s=1:listenersN;
                load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' num2str(s) ],'b_l','b_ls','b_s','F_l','F_ls','F_s2l','F_l2l','F_s','r2_ls','r2_l','r2_s','rnames');
                null(:,:,s)=F_s2l;
            end
            
            keptT=find(sum(isnan(real),1)==0 & sum(nansum(isnan(null),3)~=listenersN,1)~=0);
            p=nan(size(real));
            p(:,keptT)=nansum(null(:,keptT,:)>real(:,keptT),3)/listenersN;
            [~,ris]=sort(nansum(p==0,2),'descend');
            rnames(ris(1:5))
            figure;
          
            plot(nansum(p==0,2));
              title([exp ', binSize=' num2str(binSize) ', lags=' num2str(min(lags)) '-' num2str(max(lags))]);
            %             for ri=1:length(rnames);
            %                 for ti=1:length(keptT);
            %                     [~,p(ri,keptT(ti))]=ttest(squeeze(null(ri,keptT(ti),:)),real(ri,keptT(ti)),'tail','left');
            %                 end
            %             end
            % sig_fdr=fdr0(p(:),0.05);
            % sig_fdr=reshape(sig_fdr,size(p));
            set(gca,'xtick',1:length(rnames),'xticklabels',rnames)
            xtickangle(70)
            grid on
            save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/granger_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats'],'F_s2l','p','lags','rnames','keptT');
        end
        
        clear null
    end
end

