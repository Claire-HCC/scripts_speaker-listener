function wholeBrain_phasePerm

loc='cluster';
set_parameters;


for ei=[3 4];%1:2;%1:4;
    exp=experiments{ei};
    
    % mkdir(sprintf('%s/%s/fmri/timeseries/tr/wholeBrain/perm/',expdir,exp)
    
    f= sprintf('%s/%s/fmri/timeseries/tr/wholeBrain/speaker_zscore.mat',expdir,exp);
    load(f,'data','keptvox');
    data_real=data;
    
    for perm=1:1000;
        
        rng(perm)
        data=phase_rand2(data_real',1);
        data=data';
        
        outf= sprintf('%s/%s/fmri/timeseries/tr/wholeBrain/perm/speaker_zscore_permPhase%04d.mat',expdir,exp,perm);
        save(outf,'data','keptvox','-v7.3');
        clear data
    end
end

