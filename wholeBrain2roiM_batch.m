clear all;

% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

tic % 15 min
for ei=2:4;
    exp=experiments{ei};
    
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    gdata_listener=gdata;
    [~,tn,subjn]=size(gdata);
    clear gdata
    
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/speaker.mat',expdir,exp,timeUnit);
    load(f,'data','keptvox');
    data_speaker=data;
    clear data
    
    gdata=nan([length(rnames) tn subjn]);
    data=nan([length(rnames) tn]);
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        
        fr = sprintf('%s/roi_mask/%s/mat/%s',expdir,froidir,rname);
        load(fr,'roimask');
        
        if sum(roimask(keptvox))>10;
            gdata(ri,:,:)=nanmean(gdata_listener(logical(roimask(keptvox)),:,:));
            data(ri,:,:)=nanmean(data_speaker(logical(roimask(keptvox)),:,:));
        end
    end
    save([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_roisM' ],'data');
    save([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_roisM' ],'gdata');
    
    clear data gdata
end

