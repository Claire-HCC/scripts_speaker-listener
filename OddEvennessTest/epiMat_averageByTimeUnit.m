clearvars
set_parameters;


clear all
close all
loc='mypc';
set_parameters

role='listener';
rnames={'HG_L','vPCUN','STC_L'};

timeUnit='tr';

for ei=1:2;%;%:4;
    exp=experiments{ei};
    load([expdir exp '\sound\listenerEvents.mat']);
    eventn=length(listenerEvents_trVector(listenerEvents_trVector~=0));
    
    for ri=2:length(rnames);
        rname=rnames{ri};
        
        clear data
        
        for s=1:18;
            exp=experiments{ei};
            epi=load(sprintf('%s%s/fmri/mat/wholeBrain/%s%02d.mat',expdir,exp,role,s));
            
            epi=epi.data;
            
            rmask=load([expdir '/roi_mask/mat/' rname '.mat']);
            rmask=rmask.roimask;
            data=(epi(rmask==1,:));
            data_zscore=zscore(data,0,2);
            
            if strcmp(timeUnit,'tr')==0;
                for eventi=1:eventn;
                    i=(listenerEvents_trVector==eventi);
                    data_zscore_timeUnit(:,eventi,s)=mean(data_zscore(:,i),2);
                end
            else
                data_zscore_timeUnit(:,:,s)=data_zscore;
            end
            
        end
        data=data_zscore_timeUnit;
        save([expdir exp '/fmri/mat/roi/' timeUnit '/' timeUnit '_' role '_' rname '.mat'],'data');
        clear  data_zscore_timeUnit;
    end
end


