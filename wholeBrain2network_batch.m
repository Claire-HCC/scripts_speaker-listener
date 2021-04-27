function wholeBrain2network_batch(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','DMNb'};
networks=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
networks=strrep({networks.name},'.mat','');

tic % 15 min

exp=exp_parameters.experiments{ei}

mkdir([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir ]);
f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
load(f,'gdata','keptvox');
keptvox_wholeBrain_L=keptvox;
gdata_listener=gdata;
clear gdata

for ni=1:length(networks);
    network=networks{ni};
    fr = sprintf('%s/roi_mask/%s/mat/%s',expdir,froidir,network);
    load(fr,'roimask');
    
    if sum(roimask(keptvox_wholeBrain_L))>10;
        clear gdata
        gdata=gdata_listener(logical(roimask(keptvox_wholeBrain_L)),:,:);
        
        if ~isempty(gdata);
            keptvox=keptvox_wholeBrain_L(ismember(keptvox_wholeBrain_L,find(roimask)));
            save([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network ],'gdata','keptvox');
        end
        
    end
end

