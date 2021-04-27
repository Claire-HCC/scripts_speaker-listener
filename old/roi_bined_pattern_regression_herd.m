clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

binSize_tested=[10 20 30 40]; % tr;
lags_tested={-10:10 -10:-4, -10:-1, -20:-4, -30:-4,-10:-3};
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
rnames_table=table2array(roi_table(:,3));
roi_ids=cell2mat(table2array(roi_table(:,1)));
crop_start=10;

for ei=1:4;
    exp=experiments{ei};
    % mkdir([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/']);
    
    for binSizei=1:length(binSize_tested);
        binSize=binSize_tested(binSizei);
        
        for lagi=6;%1:length(lags_tested);
            lags=lags_tested{lagi};
            
            load([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_each/lag-10-10_classification' ],'pfdr');
            sig_classificaition=pfdr<.05;
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
            r2_sl=r2;
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
            r2_ll=r2;
            r2_llm=nanmean(r2_ll,3);
            
            keptT=[min(find(nansum(r2_sl,1)>0)):max(find(nansum(r2_sl)>0))];
            
            herd(:,1)=corr_col(r2_sl(:,keptT)',r2_llm(:,keptT)');
            
            [~,tn,listenerN]=size(r2_ll);
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL'],'r2');
            r2_sl=r2;
            
            load(sprintf('%s/%s/fmri/pattern/regression/%s/roi/%s/LL_leave1out_bined/perm/binSize%d_lag0-0_permSL',expdir,exp,timeUnit,froidir,binSize),'r2');
            r2_ll=r2;
            
            herd_null=nan([length(rnames) listenerN]);
            for perm=1:size(r2_ll,3);
                herd_null(:,perm)=corr_col(r2_sl(:,keptT,perm)',r2_ll(:,keptT,perm)');
            end
            
            herd_z=real(atanh(herd));
            herd_null_z=real(atanh(herd_null));
            
            for ri=1:length(rnames_table);
                rname=rnames_table(ri);
                [~,p(ri,1),~,stats]=ttest(herd_null_z(ri,:),herd_z(ri),'tail','left');
                t(ri,1)=-stats.tstat;
            end
            
            %%¢ø¢÷¢û
            %            % p(sig_classificaition==0,:)=NaN;
            %             pfdr=nan(size(p));
            %             if sum(~isnan(p))>0;
            %                 [~,~,pfdr(~isnan(p))]=fdr(p(~isnan(p)));
            %             end
            %
            %             disp([exp ', ' binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))])
            %             rnames_table(pfdr<.05)
            %
            %             save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'herd','herd_null','pfdr','t','p','rnames','sig_classificaition');
            %
            %             nii=roiTable2wholeBrainNii_mor([roi_ids(pfdr<.05), t(pfdr<.05)]);
            %             save_nii(nii,[expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_tmap_pfdr.nii']);
            %
            
%             %% neg
%             for ri=1:length(rnames_table);
%                 rname=rnames_table(ri);
%                 [~,p(ri,1),~,stats]=ttest(herd_null_z(ri,:),herd_z(ri),'tail','right');
%                 t(ri,1)=-stats.tstat;
%             end
%             
%             % p(sig_classificaition==0,:)=NaN;
%             pfdr=nan(size(p));
%             if sum(~isnan(p))>0;
%                 [~,~,pfdr(~isnan(p))]=fdr(p(~isnan(p)));
%             end
%             
%             disp([exp ', ' binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))])
%             rnames_table(pfdr<.05)
%             
%             save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_neg.mat'],'herd','herd_null','pfdr','t','p','rnames','sig_classificaition');
%             
%             nii=roiTable2wholeBrainNii_mor([roi_ids(pfdr<.05), t(pfdr<.05)]);
%             save_nii(nii,[expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_tmap_pfdr_neg.nii']);
            
                        %% both
            for ri=1:length(rnames_table);
                rname=rnames_table(ri);
                [~,p(ri,1),~,stats]=ttest(herd_null_z(ri,:),herd_z(ri));
                t(ri,1)=-stats.tstat;
            end
            
            % p(sig_classificaition==0,:)=NaN;
            pfdr=nan(size(p));
            if sum(~isnan(p))>0;
                [~,~,pfdr(~isnan(p))]=fdr(p(~isnan(p)));
            end
            
            disp([exp ', ' binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))])
            rnames_table(pfdr<.05)
            
            save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_both.mat'],'herd','herd_null','pfdr','t','p','rnames','sig_classificaition');
            
            nii=roiTable2wholeBrainNii_mor([roi_ids(pfdr<.05), t(pfdr<.05)]);
            save_nii(nii,[expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_tmap_pfdr_both.nii']);
            
        end
    end
end

