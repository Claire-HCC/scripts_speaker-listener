function wholeBrain_phasePerm(perm)

loc='cluster';
set_parameters;
rng(perm)

for ei=[1 2 4];%1:2;%1:4;
    exp=experiments{ei};
    clear data
    f= sprintf('%s/%s/fmri/timeseries/tr/wholeBrain/speaker.mat',expdir,exp);
    load(f,'data','keptvox');
    data_real=data;
    
    data=phase_rand2(data_real',1);
    data=data';
    
    outf= sprintf('%s/%s/fmri/timeseries/tr/wholeBrain/perm/speaker_permPhase%04d.mat',expdir,exp,perm);
    save(outf,'data','keptvox','-v7.3');

    outf= sprintf('%s/%s/fmri/timeseries/tr/wholeBrain/perm/zscore_speaker_permPhase%04d.mat',expdir,exp,perm);
    data=zscore(data,0,2);
    save(outf,'data','keptvox','-v7.3');
end

