function roi_permPhase_L(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
permN=10000;
rname=rnames{ri};

for ei=1:4;%:4;
    exp=experiments{ei};
    
    if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
        
        [roi_voxn tn listenerN]=size(gdata);
        gdata_perm=nan([size(gdata) permN]);
        
        for perm=1:permN;
            rng(perm)
            
            for si=1:listenerN;
                if si==1;
                    [temp, rp]= (phase_rand2(gdata(:,:,si)',1));
                    gdata_perm(:,:,si,perm)=temp';
                else
                    gdata_perm(:,:,si,perm)= (phase_rand2(gdata(:,:,si)',1,rp))';
                end
            end
        end
    end
    gdata=gdata_perm;
    f= sprintf('%s/%s/fmri/timeseries/%s/roi/%s/permPhase/listenerAll_zscore_%s.mat',expdir,exp,timeUnit,froidir,rname);
    save(f,'gdata','-v7.3');
end


