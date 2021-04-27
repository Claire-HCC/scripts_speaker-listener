loc='cluster';
set_parameters;
exp='simulation';

simN=30;
crop_start=25;
crop_end=20;

filelist = dir(sprintf('%s/simulation/sim%03d/*mat',expdir,simN));
filelist = filelist(~endsWith({filelist.name}, '_r.mat'));
filelist = filelist(~endsWith({filelist.name}, '_peaks.mat'));
filelist={filelist.name}';

for fi=1:length(filelist);
    if exist(sprintf('%s/simulation/sim%03d/%s',expdir,simN ,[strrep(filelist{fi},'.mat','') '_r.mat']))==0;
        
        load(sprintf('%s/simulation/sim%03d/%s',expdir,simN,filelist{fi}),'v_hrf_tr');
        [levn, tn,~]=size(v_hrf_tr);
        
        keptT=(crop_start+1):(tn-crop_end);
        v_hrf_tr=v_hrf_tr(:,keptT,:);
        tn=size(v_hrf_tr,2);
        
        % compute lagcorrelation
        r=[];
        r_timeReversed=[];
        
        for simi=1:(simN);
            for sdi=1:levn;
                for tgi=1:levn;
                    if simi~=simN;
                        r(sdi,tgi,:,simi)=circularlagcorr(v_hrf_tr(sdi,:,simi),v_hrf_tr(tgi,:,simi+1),-floor((tn-1)/2):floor((tn-1)/2));
                        r_timeReversed(sdi,tgi,:,simi)=circularlagcorr(v_hrf_tr(sdi,:,simi),fliplr(v_hrf_tr(tgi,:,simi+1)),-floor((tn-1)/2):floor((tn-1)/2));
                    else
                        r(sdi,tgi,:,simi)=circularlagcorr(v_hrf_tr(sdi,:,simi),v_hrf_tr(tgi,:,1),-floor((tn-1)/2):floor((tn-1)/2));
                        r_timeReversed(sdi,tgi,:,simi)=circularlagcorr(v_hrf_tr(sdi,:,simi),fliplr(v_hrf_tr(tgi,:,1)),-floor((tn-1)/2):floor((tn-1)/2));
                    end
                end
            end
        end
        
        save(sprintf('%s/simulation/sim%03d/%s',expdir,simN ,[strrep(filelist{fi},'.mat','') '_crossStories_r.mat']),...
            'r_timeReversed','r','keptT');
    end
end