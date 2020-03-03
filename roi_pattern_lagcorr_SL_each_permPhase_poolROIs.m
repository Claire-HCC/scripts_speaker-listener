
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

lags_tested={-10:10,  -30:30};
permN=1000;

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            if exist([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLeach/perm/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' rname '.mat']);
                load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLeach/perm/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' rname],'r','lags','keptT');
                if ri==1;
                    r_perm=r;
                elseif length(size(r))==4;
                    r_perm(ri,:,:,:)=r;
                else
                    r_perm(ri,:,:,:)=NaN;
                end
            end
            r=r_perm;
            save([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLeach/perm/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'rnames','r','lags','keptT');
        end
    end
end

