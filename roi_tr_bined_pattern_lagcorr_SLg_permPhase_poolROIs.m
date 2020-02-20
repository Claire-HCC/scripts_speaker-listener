
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
binSize=40;
lags_tested={-10:10,  -30:30};
permN=1000;

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        load([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'r');
        [~,tn,~]=size(r);
        r_perm=nan([length(rnames) tn  length(lags) permN]);
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            if exist([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/perm/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname '.mat']);
                load([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/perm/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname],'r','lags','keptT');
                r_perm(ri,:,:,:)=r;
            end
        end
        r=r_perm;
        save([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/perm/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'rnames','r','lags','keptT','-v7.3');
        
    end
end

