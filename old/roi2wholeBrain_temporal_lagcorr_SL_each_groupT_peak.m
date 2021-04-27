%% ttest across subject
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -40:40};
seed='vPCUN';
for ei=3:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
            load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2wholeBrain/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','keptT','keptvox');
            rz=atanh(r);
            rzm=nanmean(rz,3);
            p=nan(size(rzm));
            t=p;
            for vi=1:length(keptvox);
                [~,p(vi,:),~,stats]=ttest(squeeze(rz(vi,:,:))',0,'tail','right');
                t(vi,:)=stats.tstat;
            end

            pfwe=p*(sum(~isnan(p(:))));
            
            pfdr=nan(size(p(:)));
            [ ~,~,pfdr(~isnan(p(:)))]=fdr(p(~isnan(p(:))));
            pfdr=reshape(pfdr,size(p));
            
            vis=find(sum(~isnan(rzm),2)~=0);
            peaks=nan([length(keptvox) 1 ]);
            peakLags=peaks;
            peaks_pfwe=nan([length(keptvox) 1 ]);
            peakLags_pfwe=peaks_pfwe;
            peaks_pfdr=nan([length(keptvox) 1 ]);
            peakLags_pfdr=peaks_pfdr;
            rz_tempfwe=rzm;
            rz_tempfwe(pfwe>.05)=NaN;
            
            rz_tempfdr=rzm;
            rz_tempfdr(pfdr>.05)=NaN;
            for i=1:length(vis);
                vi=vis(i);
                [~, peakLagi]=max(rzm(vi,:),[],2);;%max(abs(rzm(vi,:)),[],2);
                peakLags(vi,1)=(lags(peakLagi));
                peaks(vi,1)=rzm(vi,peakLagi);
                
                if min(pfdr(vi,:))<.05;
                    [~, peakLagi]=max(rz_tempfdr(vi,:),[],2);;%max(abs(rz_tempfdr(vi,:)),[],2);
                    peakLags_pfdr(vi,1)=(lags(peakLagi));
                    peaks_pfdr(vi,1)=rz_tempfdr(vi,peakLagi);
                end
                
                if min(pfwe(vi,:))<.05;
                    [~, peakLagi]=max(rz_tempfwe(vi,:),[],2);%max(abs(rz_tempfwe(vi,:)),[],2);
                    peakLags_pfwe(vi,1)=(lags(peakLagi));
                    peaks_pfwe(vi,1)=rz_tempfwe(vi,peakLagi);
                end
            end
            
              mat=nan([length(keptvox) 1]);
            mat(keptvox)=peakLags+0.00000001; % otherwise, peaklag0 won't show up in xjview
            mat(keptvox(min(p,[],2)<.05))=NaN;
            nii=mat2nii(mat);
            nii.img(1,1,[1 2])=[-1 1];
            save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2wholeBrain/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_posPeakLags_p05.nii']);
            
            
              mat=nan([length(keptvox) 1]);
            mat(keptvox)=peakLags+0.00000001; % otherwise, peaklag0 won't show up in xjview
            mat(keptvox(isnan(peakLags)))=NaN;
            nii=mat2nii(mat);
            nii.img(1,1,[1 2])=[-1 1];
            save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2wholeBrain/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_posPeakLags.nii']);
            
            
            mat=nan([length(keptvox) 1]);
            mat(keptvox)=peakLags_pfwe+0.00000001; % otherwise, peaklag0 won't show up in xjview
            mat(keptvox(isnan(peakLags_pfwe)))=NaN;
            nii=mat2nii(mat);
            nii.img(1,1,[1 2])=[-1 1];
            save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2wholeBrain/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_posPeakLags_pfwe.nii']);
            
            mat=nan([length(keptvox) 1]);
            mat(keptvox)=peakLags_pfdr+0.00000001;
            mat(keptvox(isnan(peakLags_pfdr)))=NaN;
            nii=mat2nii(mat);
            nii.img(1,1,[1 2])=[-1 1];
            save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2wholeBrain/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_posPeakLags_pfdr.nii']);
            
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2wholeBrain/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT_peaks' ],...
                't','r','rzm','lags','keptT','p','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_pfdr','peakLags_pfdr','peakLags','peaks');
        end
    end


