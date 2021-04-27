function wholeBrainPermPahse2roi_batch(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
permN=10000;
tic % 15 min
rname=rnames{ri};

lags=-10:10;
for ei=[4];%[3 4];%1:2;%1:4;
    exp=experiments{ei};
    % mkdir(sprintf('%s/%s/fmri/timeseries/%s/roi/%s/permPhase/',expdir,exp,timeUnit,froidir));
    
    if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'])>0 & exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/permPhase/speaker_zscore_' rname '.mat'])==0;
        
        fr = sprintf('%s/roi_mask/%s/mat/%s',expdir,froidir,rname);
        load(fr,'roimask');
        
        for perm=1:permN;
            f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/perm/speaker_zscore_permPhase%04d.mat',expdir,exp,timeUnit,perm);
            load(f,'data','keptvox');
            
            if sum(roimask(keptvox))>10;
                data_temp(:,:,perm)=data(logical(roimask(keptvox)),:);
            end
        end
        data=data_temp;
        f= sprintf('%s/%s/fmri/timeseries/%s/roi/%s/permPhase/speaker_zscore_%s.mat',expdir,exp,timeUnit,froidir,rname);
        save(f,'data','-v7.3');
        clear data_temp
        
    end
end