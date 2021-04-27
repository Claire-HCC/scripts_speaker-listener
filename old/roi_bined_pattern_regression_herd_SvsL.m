clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

binSize_tested=[1 5 10 20 30 ]; % tr;
lags_tested={-10:10, -10:-3, -10:-1,-20:-3};
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
rnames_table=table2array(roi_table(:,3));
roi_ids=cell2mat(table2array(roi_table(:,1)));
crop_start=10;

for ei=3:4;
    exp=experiments{ei};
   % mkdir([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd_SvsL/']);
    
    for binSizei=1:length(binSize_tested);
        binSize=binSize_tested(binSizei);
        
        for lagi=2;%1:length(lags_tested);
            lags=lags_tested{lagi};
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
            r2_sl=r2;
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
            r2_ll=r2;
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_minus1Leave1out_bined/binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
            r2_ll0=r2;
            r2_ll0_m=nanmean(r2_ll0,4);
            
            keptT=[min(find(~isnan(r2_sl(1,:,1)))):max(find(~isnan(r2_sl(1,:,1))))];
            
            [~,tn,listenerN]=size(r2_ll);
            
            herd=nan([length(rnames) listenerN]);
            herd_LL=nan([length(rnames) listenerN]);
            for s=1:size(r2_ll,3);
                herd(:,s)=corr_col(r2_sl(:,keptT,s)',r2_ll0_m(:,keptT,s)');
                herd_LL(:,s)=corr_col(r2_ll(:,keptT,s)',r2_ll0_m(:,keptT,s)');
            end
            
            herd_z=(atanh(herd));
            herd_LL_z=(atanh(herd_LL));
           
            %% pos
            for ri=1:length(rnames_table);
                rname=rnames_table(ri);
                [~,p(ri,1),~,stats]=ttest2(herd_LL_z(ri,:),herd_z(ri,:),'tail','left','Vartype','unequal');
                t(ri,1)=-stats.tstat;
            end
           
            % p(sig_classificaition==0,:)=NaN;
            pfdr=nan(size(p));
            if sum(~isnan(p))>0;
                [~,~,pfdr(~isnan(p))]=fdr(p(~isnan(p)));
            end
            
            disp([exp ', ' binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))])
            rnames_table(pfdr<.05)
            
            save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd_SvsL/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'herd','herd_LL','pfdr','t','p','rnames');
            
            nii=roiTable2wholeBrainNii_mor([roi_ids(pfdr<.05), t(pfdr<.05)]);
            save_nii(nii,[expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd_SvsL/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_tmap_pfdr.nii']);
            
            
            %% neg
            for ri=1:length(rnames_table);
                rname=rnames_table(ri);
                [~,p(ri,1),~,stats]=ttest2(herd_LL_z(ri,:),herd_z(ri,:),'tail','right','Vartype','unequal');
                t(ri,1)=-stats.tstat;
            end
            
            % p(sig_classificaition==0,:)=NaN;
            pfdr=nan(size(p));
            if sum(~isnan(p))>0;
                [~,~,pfdr(~isnan(p))]=fdr(p(~isnan(p)));
            end
            
            disp([exp ', ' binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))])
           % rnames_table(pfdr<.05)
            
            save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd_SvsL/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_neg.mat'],'herd','herd_LL','pfdr','t','p','rnames');
            
            nii=roiTable2wholeBrainNii_mor([roi_ids(pfdr<.05), t(pfdr<.05)]);
            save_nii(nii,[expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd_SvsL/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_tmap_pfdr_neg.nii']);
            
            %% both
            for ri=1:length(rnames_table);
                rname=rnames_table(ri);
                [~,p(ri,1),~,stats]=ttest2(herd_LL_z(ri,:),herd_z(ri,:),'Vartype','unequal');
                t(ri,1)=-stats.tstat;
            end
            
            % p(sig_classificaition==0,:)=NaN;
            pfdr=nan(size(p));
            if sum(~isnan(p))>0;
                [~,~,pfdr(~isnan(p))]=fdr(p(~isnan(p)));
            end
            
            disp([exp ', ' binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))])
            % rnames_table(pfdr<.05)
            
            save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd_SvsL/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_both.mat'],'herd','herd_LL','pfdr','t','p','rnames');
            
            nii=roiTable2wholeBrainNii_mor([roi_ids(pfdr<.05), t(pfdr<.05)]);
            save_nii(nii,[expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd_SvsL/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_tmap_pfdr_both.nii']);
            
        end
    end
end

