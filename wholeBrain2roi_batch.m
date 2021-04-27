function wholeBrain2roi_batch(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='aal';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

exp=exp_parameters.experiments{ei}

mkdir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir ]);
f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
load(f,'gdata','keptvox');
keptvox_wholeBrain_L=keptvox;
gdata_listener=gdata;
clear gdata

for ni=1:length(rnames);
    rname=rnames{ni};
    fr = sprintf('%s/roi_mask/%s/mat/%s',expdir,froidir,rname);
    load(fr,'roimask');
    
    if sum(roimask(keptvox_wholeBrain_L))>10;
        clear gdata
        gdata=gdata_listener(logical(roimask(keptvox_wholeBrain_L)),:,:);
        
        if ~isempty(gdata);
            keptvox=keptvox_wholeBrain_L(ismember(keptvox_wholeBrain_L,find(roimask)));
            save([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname ],'gdata','keptvox');
        end
        
    end
end

