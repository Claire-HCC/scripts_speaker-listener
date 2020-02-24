
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -30:30};

for ei=1:4;%1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'sig_fdr');
        regression_mask=sig_fdr;
        
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/perm/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'r');
        r_perm=r;
        
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags','keptT');
        p=mean(r_perm>r,3);
        %  p(regression_mask==0,:)=NaN;
        p(isnan(r))=NaN;
        
        p_fdr=nan(size(p(:)));
        if sum(~isnan(p))>0
            [~,~,p_fdr(~isnan(p))]=fdr(p(~isnan(p)));
            p_fdr=reshape(p_fdr,size(p));
        end
        
        ris=find(sum(~isnan(r),2)~=0);
        peak=nan([length(rnames) 1 ]);
        peakLags=peak;
        [peak(ris) peakLagi]=max(r(ris,:),[],2);
        peakLags(ris,1)=(lags(peakLagi))';
        
        peak_pfdr=nan([length(rnames) 1 ]);
        peakLags_pfdr=peak_pfdr;
        r_temp=r;
        r_temp(p_fdr>.05)=NaN;
        ris=find(sum(~isnan(r_temp),2)~=0);
        [peak_pfdr(ris) peakLagi_pfdr]=max(r_temp(ris,:),[],2);
        peakLags_pfdr(ris,1)=(lags(peakLagi_pfdr))';
        
        subplot(2,4,ei);
        imagesc(r);
        subplot(2,4,ei+1);
        imagesc(r_temp);
       
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peak_pfdr)),   peak_pfdr(~isnan( peak_pfdr))]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/'  froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_peak_pfdr.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags)),   peakLags(~isnan(peakLags))+0.00000001]);
        nii.img(1,1,1)=1;
        save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/'  froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_peakLags.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfdr)),   peakLags_pfdr(~isnan( peakLags_pfdr))+0.00000001]);
        nii.img(1,1,1)=1;
        save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/'  froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_peakLags_pfdr.nii']);
        
        save([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_peaks' ],'rnames','r','lags','keptT','r_perm','p','p_fdr','peak','peakLags','peakLags_pfdr');
    end
end

