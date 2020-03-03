
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
lags_tested={-10:10,  -40:40};
permN=1000;
for ei=1:4;
    exp=experiments{ei};
    
    for lagi=2%:length(lags_tested);
        lags=lags_tested{lagi};
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            if exist([expdir '/' exp '/fmri/pattern/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname '.mat']);
                load([expdir '/' exp '/fmri/pattern/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname],'r','lags','keptT');
                if ri==1;
                    r_perm=r(:,:,1:permN);
                else
                    r_perm(ri,:,:)=r(:,:,1:permN);
                end
            else
                 r_perm(ri,:,:)=NaN;
            end
        end
        r=r_perm;
        save([expdir '/' exp '/fmri/pattern/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'rnames','r','lags','keptT');
        clear r_perm;
    end
end

