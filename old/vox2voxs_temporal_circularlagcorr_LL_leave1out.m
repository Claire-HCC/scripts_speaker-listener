function vox2voxs_temporal_circularlagcorr_LL_leave1out(vi_start)

vn=500;
vis_start=[0:vn:8000];
vi_start=vis_start(vi_start);

loc='cluster';
set_parameters;
timeUnit='tr' ;
crop_start=25;
crop_end=20;

for ei=1;% forgot vox 1:11 finished
    
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2vox/LL_ggut/isc'   ],'peakLags','peaks','keptvox');
    maskPeakLag0=zeros(voxn,1);
    thr=sort(peaks(peakLags==0),'descend');
    thr=thr(round(length(keptvox)*0.3));
    maskPeakLag0(keptvox(peaks>thr & peakLags'==0))=1;
    clear peakLags peaks keptvox
    
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_ggut/perm/']);
    
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    gdata=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
    keptvox=keptvox(ismember(keptvox,find(maskPeakLag0)),:,:);
    gdata(:,:,subjects_excluded{ei})=NaN;
    
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    
    for vi=(vi_start+[1:vn]);
        fout=([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_ggut/perm/isfc_voxind' num2str(keptvox(vi)) '.mat'  ]);
        disp(vi)
        if ~exist(fout);
            r=nan([length(keptvox) length(lags) listenerN   ]);
            for si=1:listenerN;
                othersi=1:listenerN;
                othersi=othersi(othersi~=si);
                
                y=nanmean(zscore(gdata(:,keptT,othersi),0,2),3)';
                x=zscore(gdata(vi,keptT,si),0,2)';
                x=repmat(x,1,size(y,2));
                r(:,:,si)=circularlagcorr_byCol(x,y,lags)';
            end
            rzm=nanmean(atanh(r),3);
            voxind=keptvox(vi);
            save(fout,'rzm','keptvox','keptT','voxind');
        end
    end
end


