
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
lags_tested={-10:10,  -40:40};
seeds={'HG_L','vPCUN'};

for ei=[1 2 4];%1:4;
    exp=experiments{ei};
    
    for lagi=1;%1:2%:length(lags_tested);
        lags=lags_tested{lagi};
        
        for sdi=1:length(seeds);
            seed=seeds{sdi};
            for ri=1:length(rnames);
                rname=rnames{ri};
                if exist([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/perm/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname '.mat']);
                    load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/perm/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname],'r','lags','keptT');
                    if ri==1;
                        r_perm=r;
                    else
                        r_perm(ri,:,:,:)=r;
                    end
                else
                    r_perm(ri,:,:,:)=NaN;
                end
            end
            r=r_perm;
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/perm/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'rnames','r','lags','keptT','-v7.3');
            clear r_perm;
        end
    end
end

