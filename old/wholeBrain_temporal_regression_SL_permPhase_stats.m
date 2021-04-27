clear all
close all

% loc='cluster';
set_parameters;
timeUnit='tr' ;
lags_tested={-4:4, -10:10};


for ei=3:4;%1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain//SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase.mat' ],'r2');
        r2_perm=r2;
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain//SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat' ],'r2','b','keptvox');
        
        p=nanmean(r2_perm>r2,2);
        pfdr=nan(size(p));
        [~,~,pfdr]=fdr(p);
        
        save([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'r2','p','pfdr','r2_perm','keptvox','-v7.3');
        
        mat=nan(voxn,1);
        mat(keptvox(pfdr<.05))=r2(pfdr<.05);
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_r2_pfdr.mat' ])
    end
end


