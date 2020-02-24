clear all;

loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

tic % 15 min

for ei=3:4;%1:4;
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir ]);
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll_zscore.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    gdata_listener=gdata;
    
    clear gdata
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/speaker_zscore.mat',expdir,exp,timeUnit);
    load(f,'data','keptvox');
    data_speaker=data;
        
    for ri=1:length(rnames);
        rname=rnames{ri};
        fr = sprintf('%s/roi_mask/%s/mat/%s',expdir,froidir,rname);
        load(fr,'roimask');
        
        if sum(roimask(keptvox))>10;
            clear gdata
            gdata=gdata_listener(logical(roimask(keptvox)),:,:);
            if ~isempty(gdata);

                save([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname ],'gdata');
            end
            
            clear data
            data=data_speaker(logical(roimask(keptvox)),:,:);
            if ~isempty(data);
                save([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname ],'data');
            end
            
            clear data
        end
    end
end
beep;