clear all;

loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

tic % 15 min

for ei=[3 4];%1:2;%1:4;
    exp=experiments{ei};
    
    %  mkdir(sprintf('%s/%s/fmri/timeseries/%s/roi/%s/permPhase/',expdir,exp,timeUnit,froidir));
    for perm=1:1000;
        f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/perm/speaker_grayMasked_permPhase%04d.mat',expdir,exp,timeUnit,perm);
        load(f,'data','keptvox');
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            fr = sprintf('%s/roi_mask/%s/mat/%s',expdir,froidir,rname);
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname '.mat'])>0 & exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/permPhase/speaker_' rname '.mat'])==0;
                load(fr,'roimask');
                
                if sum(roimask(keptvox))>10;
                    data_temp.(rname)(:,:,perm)=data(logical(roimask(keptvox)),:);
                end
            end
        end
    end
    
    rnames=fieldnames(data_temp)
    for ri=1:length(rnames);
        rname=rnames{ri};
        data=data_temp.(rname);
        
        f= sprintf('%s/%s/fmri/timeseries/%s/roi/%s/permPhase/speaker_grayMasked_%s.mat',expdir,exp,timeUnit,froidir,rname);
        save(f,'data');
        
        data=zcore(data,0,2)
        f= sprintf('%s/%s/fmri/timeseries/%s/roi/%s/permPhase/zscore_speaker_grayMasked_%s.mat',expdir,exp,timeUnit,froidir,rname);
        save(f,'data');
        
    end
    clear data data_temp
end