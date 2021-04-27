function voxs2voxs_temporal_corr_LL_selfself

loc='cluster';
set_parameters;
timeUnit='tr' ;
crop_start=25;
crop_end=20;
exp='pieman_rest';
perc=0.15;
ei=find(ismember(exp_parameters.experiments,exp));
%     try
%     rmdir([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/'],'s')
%     end

mkdir([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/']);

% f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll_isc%dPercMasked.mat',expdir,'pieman_old',timeUnit,perc*100);
% load(f,'mask');

f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
load(f,'gdata','keptvox');
% 
% % masked with pieman_old;
% gdata=gdata(logical(mask(keptvox)) ,:,:);
% keptvox=keptvox(logical(mask(keptvox)) ,:,:);

gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;

[~,tn,listenerN]=size(gdata);
keptT=(crop_start+1):(tn-crop_end);

rzm=nan([length(keptvox) length(keptvox)   ]);
for si=1:listenerN;
    
    x=zscore(gdata(:,keptT,si),0,2)';
    
    % minus 0.00001 to avoud inf due to atanh
    rzm=nansum(cat(3,rzm,atanh(corr(x,x)-0.00001)),3);
end
rzm=rzm/sum(squeeze(~isnan(gdata(1,1,:))));
save([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc'],'rzm','keptvox','keptT','-v7.3');

% save([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_isc' num2str(100*perc) 'PercMasked'],'rzm','keptvox','keptT','-v7.3');

