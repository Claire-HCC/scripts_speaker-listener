clear all
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

crop_start=10;
lags_tested={-10:10, -30:30};
speakerSeed='vPCUN';
permN=1000;

for ei=2:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for perm=1:permN;
            
            load(sprintf('%s/%s/fmri/temporal_lagcorr/%s/roi2wholeBrain/%s/SLg/perm/%s_lag%d-%d_permPhase%04d',expdir,exp,timeUnit,froidir,speakerSeed,min(lags),max(lags),perm),'r','lags','keptvox','keptT');
            r_perm(:,:,perm)=r;
        end
        r=r_perm;
        save(sprintf('%s/%s/fmri/temporal_lagcorr/%s/roi2wholeBrain/%s/SLg/perm/%s_lag%d-%d_permPhase',expdir,exp,timeUnit,froidir,speakerSeed,min(lags),max(lags)),'r','lags','keptvox','keptT','-v7.3');
       clear r_perm 
    end
end
    
