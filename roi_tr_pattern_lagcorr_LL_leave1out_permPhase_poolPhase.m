
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

lags_tested={-10:10,  -30:30};
permN=1000;

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load(sprintf('%s/%s/fmri/pattern_lagcorr/%s/roi/%s/LL_leave1out/perm/lag%d-%d_permPhase%04d',expdir,exp,timeUnit,froidir,min(lags),max(lags),1),'r','lags','rnames','keptT');
        [~,~,listenerN]=size(r);
        r_perm=nan([length(rnames) length(lags) listenerN permN]);
        
        for perm=1:permN;
            
            rname=rnames{ri};
            if exist(sprintf('%s/%s/fmri/pattern_lagcorr/%s/roi/%s/LL_leave1out/perm/lag%d-%d_permPhase%04d.mat',expdir,exp,timeUnit,froidir,min(lags),max(lags),perm));
                load(sprintf('%s/%s/fmri/pattern_lagcorr/%s/roi/%s/LL_leave1out/perm/lag%d-%d_permPhase%04d',expdir,exp,timeUnit,froidir,min(lags),max(lags),perm),'r','lags','rnames','keptT');
                r_perm(:,:,:,perm)=r;
            end
        end
        
        r=r_perm;
        save([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/LL_leave1out/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'rnames','r','lags','keptT');
    end
end

