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
interROItypes={'roi2rois','roi'};
seeds={'vPCUN','pANG_L','HG_L'};
interSubjTypes={'LL_leave1out','SL_each'};

for ir=1%:length(interROItypes);
    interROItype=interROItypes{ir};
    
    
    if strcmp(interROItype,'roi');
        seeds_={''};
    else
        seeds_= cellfun(@(x) [x '_'],seeds,'UniformOutput',0);
    end
    
    for is=1%:length(interSubjTypes);
        interSubjType=interSubjTypes{is};
        
        for ei=3;%1:4;
            exp=experiments{ei};
            
            for sdi=1%:length(seeds_);
                seed_=[seeds_{sdi}];
                
                for lagi=1%:length(lags_tested);
                    lags=lags_tested{lagi};
                    
                    load([expdir exp '\fmri\temporal\lagcorr\tr\' interROItype '\' froidir '\' interSubjType '\' seed_ 'lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT_peaks'],'rnames','r','lags','keptT','rzm','pfdr','peakLags_pfdr')
                    
%                     
%                     for li=1:length(lags);
%                         nii=roiTable2wholeBrainNii_mor([roi_ids rzm(:,li)]);
%                         save_nii(nii,[expdir exp '\fmri\temporal\lagcorr\tr\' interROItype '\' froidir '\' interSubjType '\' seed_ 'lag' num2str(min(lags)) '-' num2str(max(lags))   '_R_lag' num2str(lags(li)) '.nii']);
%                     end
%                     
%                     rzm(pfdr>.05)=0;
%                     for li=1:length(lags);
%                         nii=roiTable2wholeBrainNii_mor([roi_ids rzm(:,li)]);
%                         save_nii(nii,[expdir exp '\fmri\temporal\lagcorr\tr\' interROItype '\' froidir '\' interSubjType '\' seed_ 'lag' num2str(min(lags)) '-' num2str(max(lags))   '_R_pfdr_lag' num2str(lags(li)) '.nii']);
%                     end
                    
                    for li=1:length(lags);
                        nii=roiTable2wholeBrainNii_mor([roi_ids(peakLags_pfdr==lags(li)) ones(sum(peakLags_pfdr==lags(li)),1)]);
                        % otherwise freesurfer label file can not be
                        % generated
                        nii.img(1:2,1:2,1:2)=1;
                        save_nii(nii,[expdir exp '\fmri\temporal\lagcorr\tr\' interROItype '\' froidir '\' interSubjType '\' seed_ 'lag' num2str(min(lags)) '-' num2str(max(lags))   '_peakLag_pfdr_lag' num2str(lags(li)) '.nii']);
                    end
                end
            end
        end
    end
end


