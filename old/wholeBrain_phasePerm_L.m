function wholeBrain_phasePerm_L(perm_start)

loc='cluster';
set_parameters;

for ei=[3 4];%1:2;%1:4;
    exp=experiments{ei};
    
    % mkdir(sprintf('%s/%s/fmri/timeseries/tr/wholeBrain/perm/',expdir,exp)
    
    f= sprintf('%s/%s/fmri/timeseries/tr/wholeBrain/listenerAll_zscore.mat',expdir,exp);
    load(f,'gdata','keptvox');
    gdata_real=gdata;
    subjn=size(gdata_real,3);
    
    for perm=(perm_start*1000-1000+1):(perm_start*1000);
        outf= sprintf('%s/%s/fmri/timeseries/tr/wholeBrain/perm/listenerAll_zscore_permPhase%04d.mat',expdir,exp,perm);
        if ~exist(outf);
            rng(perm)
            for si=1:subjn;
                
                if si==1;
                    [ temp, rp]=phase_rand2(gdata_real(:,:,si)',1);
                else;
                    [ temp]=phase_rand2(gdata_real(:,:,si)',1,rp);
                end
                gdata(:,:,si)=temp';
            end
            
            save(outf,'gdata','keptvox','-v7.3');
            clear gdata
        end
    end
end

