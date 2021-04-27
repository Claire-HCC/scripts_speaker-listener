loc='cluster';
set_parameters;
timeUnit='tr' ;
crop_start=25;
crop_end=20;
exp='pieman_rest';
perc=0.3;
ei=find(ismember(exp_parameters.experiments,exp));
%     try
%     rmdir([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/'],'s')
%     end

mkdir([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/']);

f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll_isc%dPercMasked.mat',expdir,'pieman_old',timeUnit,perc*100);
load(f,'mask');

f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
load(f,'gdata','keptvox');

% masked with pieman_old;
gdata=gdata(logical(mask(keptvox)) ,:,:);
keptvox=keptvox(logical(mask(keptvox)) ,:,:);

gdata(:,:,exp_parameters.subjects_excluded{ei})=[];

[~,tn,listenerN]=size(gdata);
keptT=(crop_start+1):(tn-crop_end);

g1=datasample(1:listenerN,round(listenerN/2));
g2=find(~ismember(1:listenerN,g1));

rzm_g1=nan([length(keptvox) length(keptvox)   ]);
rzm_g2=nan([length(keptvox) length(keptvox)   ]);
for si=1:listenerN;
    
    x=zscore(gdata(:,keptT,si),0,2)';
    
    if ismember(si,g1);
        % minus 0.00001 to avoud inf due to atanh
        rzm_g1=nansum(cat(3,rzm_g1,atanh(corr(x,x)-0.00001)),3);
    else
        rzm_g2=nansum(cat(3,rzm_g2,atanh(corr(x,x)-0.00001)),3);
    end
end
rzm_g1=rzm_g1/length(g1);
rzm_g2=rzm_g2/length(g2);

save([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_isc' num2str(100*perc) 'PercMasked_g2'],'rzm_g1','rzm_g2','keptvox','keptT','-v7.3');

