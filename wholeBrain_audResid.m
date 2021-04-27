clear all;
loc='mypc';
set_parameters;
timeUnit=['tr' ];

for ei=5;%6;
    exp=exp_parameters.experiments{ei};
    
     f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    [~,tn,listenerN]=size(gdata);
    
    load([expdir exp '/sound/' exp '_listener_audenv'],'aud')
    
    % get hrf
    xBF.dt=exp_parameters.tr(ei);
    xBF.name='hrf (with time and dispersion derivatives)';
    bf = spm_get_bf(xBF);
    
    X=[];
    for hrfi=1:3;
        X(:,hrfi)=conv(aud-mean(aud),bf.bf(:,hrfi));
    end
    % I have cropped 3 TRs from the fMRI data to account for hrf delay..
    X=X(4:(tn+3),:);
    X=[(ones(tn,1)) X];
    
    betas=[];
    resid=[];
    for si = 1:listenerN;
        for vi=1:length(keptvox);
            y=squeeze(gdata(vi,:,si))';
            [b]=regress(y,X);
            Y_aud =X(:,2:end)*b(2:end);
            resid(vi,:,si) = y - Y_aud;
            betas(vi,:,si)=b(2:end);
        end
    end
    
    gdata=resid;
    f_out = sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll_audResid.mat',expdir,exp,timeUnit);
    save(f_out,'gdata','betas','keptvox','-v7.3');
end






