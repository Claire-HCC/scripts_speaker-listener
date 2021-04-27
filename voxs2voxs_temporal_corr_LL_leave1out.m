function voxs2voxs_temporal_corr_LL_leave1out

loc='cluster';
set_parameters;
timeUnit='tr' ;
crop_start=25;
crop_end=20;
eis=[1 2 4 11 12 9 10 13 14];

for eii=2:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    
    try
    rmdir([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_leave1out/'],'s')
    
    end
    mkdir([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_leave1out/']);
    
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll_iscmasked.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    gdata(:,:,subjects_excluded{ei})=NaN;
    
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    
    r=nan([length(keptvox) length(keptvox) listenerN   ]);
    
    for si=1:listenerN;
        othersi=1:listenerN;
        othersi=othersi(othersi~=si);
        
        y=nanmean(zscore(gdata(:,keptT,othersi),0,2),3)';
        x=zscore(gdata(:,keptT,si),0,2)';
        r(:,:,si)=corr(x,y);
      %  disp(si)
    end
    rzm=nanmean(atanh(r),3);
    save([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_leave1out/isfc'],'rzm','keptvox','keptT','-v7.3');
end

