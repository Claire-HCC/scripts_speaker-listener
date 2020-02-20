
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

lags_tested={-10:10,  -30:30, -4:4};
permN=1000;

for ei=1:4;%2:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        b_perm=nan([length(rnames)  length(lags)+1 permN]);
        r2_perm=nan([length(rnames)   permN]);
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            if exist([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname '.mat']);
                load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname],'b','r2','keptT');
                b_perm(ri,:,:)=b;
                r2_perm(ri,:)=r2;
            end
            
            b=b_perm;
            r2=r2_perm;
            save([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/perm/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'rnames','b','r2','keptT');
        end
    end
end

