function wholeBrainPermPahse2roi_L_batch(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
permN=10000;
tic % 15 min
rname=rnames{ri};

for ei=[3 4];%1:2;%1:4;
    exp=experiments{ei};
    % mkdir(sprintf('%s/%s/fmri/timeseries/%s/roi/%s/permPhase/',expdir,exp,timeUnit,froidir));
    
    if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat'])>0 & exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/permPhase/listenerAll_zscore_' rname '.mat'])==0;
        
        fr = sprintf('%s/roi_mask/%s/mat/%s',expdir,froidir,rname);
        load(fr,'roimask');
        
        for perm=1:permN;
            f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/perm/listenerAll_zscore_permPhase%04d.mat',expdir,exp,timeUnit,perm);
            load(f,'gdata','keptvox');
            
            if sum(roimask(keptvox))>10;
                gdata_temp(:,:,:,perm)=gdata(logical(roimask(keptvox)),:,:);
            end
        end
        gdata=gdata_temp;
        f= sprintf('%s/%s/fmri/timeseries/%s/roi/%s/permPhase/listenerAll_zscore_%s.mat',expdir,exp,timeUnit,froidir,rname);
        save(f,'gdata','-v7.3');
        clear gdata_temp
        
    end
end