loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

lags_tested={-10:10,  -40:40};
permN=1000;

for ei=1:2;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        r_perm=nan([length(rnames)  length(lags) permN]);
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            if exist([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname '.mat']);
                load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname],'r','lags','keptT');
                r_perm(ri,:,:)=r;
            end
            
            r=r_perm;
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'rnames','r','lags','keptT','-v7.3');
        end
    end
end

