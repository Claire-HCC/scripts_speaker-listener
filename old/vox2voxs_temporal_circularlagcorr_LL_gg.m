function vox2voxs_temporal_circularlagcorr_LL_gg(batchi)
loc='cluster';
set_parameters;

vn=500;
vis_start=[0:vn:8000];
vi_start=vis_start(batchi);

timeUnit='tr' ;
crop_start=25;
crop_end=20;

for ei=12;% forgot vox 1:11 finished
    exp=experiments{ei};
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/perm/']);
    
    % peakLag0 isc mask
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2vox/LL_leave1out/isc'   ],'peakLags','peaks','keptvox');
    maskPeakLag0=zeros(voxn,1);
    thr=sort(peaks(peakLags==0),'descend');
    thr=thr(round(length(keptvox)*0.3));
    maskPeakLag0(keptvox(peaks>thr & peakLags'==0))=1;
    clear peakLags peaks keptvox
    
    %
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    [~,tn,listenerN]=size(gdata);
    gdata=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
    keptvox=keptvox(ismember(keptvox,find(maskPeakLag0)),:,:);
    gdata(:,:,subjects_excluded{ei})=NaN;
    
    keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    
    % random subject group
    if batchi==1;
        subjs_g1=[];
        subjs_g2=[];
        subjs_shuffled=randperm(listenerN);
        subjs_shuffled(ismember(subjs_shuffled,subjects_excluded{ei}))=[];
        subjs_g1=subjs_shuffled(1:round(length(subjs_shuffled)/2));
        subjs_g2=subjs_shuffled((1+round(length(subjs_shuffled)/2)):end);
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/perm/subjs_gg.mat'], 'subjs_g1', 'subjs_g2');
    else
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/perm/subjs_gg.mat'], 'subjs_g1', 'subjs_g2');
    end
    
    for vi=(vi_start+[1:vn]);
        fout=([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/perm/isfc_voxind' num2str(keptvox(vi)) '.mat'  ]);
        
        
        if ~exist(fout);
            disp(vi)
            y=nanmean(zscore(gdata(:,keptT, subjs_g2),0,2),3)';
            x=nanmean(zscore(gdata(vi,keptT, subjs_g1),0,2),3)';
            x=repmat(x,1,size(y,2));
            r=circularlagcorr_byCol(x,y,lags)';
            
            voxind=keptvox(vi);
            save(fout,'r','keptvox','keptT','voxind');
        end
    end
end


