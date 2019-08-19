
clear all;
tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

iters=1000;
tic % 15 min
win_width=25;
win_step=1;
lags=-20:20;

for ei=1%1:2;
    exp=experiments{ei};
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        fl=[expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname  ];
        
        if exist([fl '.mat'])>0;
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname ],'data');
            
            if size(gdata,1)>10 & size(data,1)>10;
                data=zscore(mean(data,1),0,2)';
                gdata=(zscore(nanmean(gdata,1),0,2));
                
                speaker=data;
                speaker=buffer(speaker,win_width,win_width-win_step);
                speaker=speaker(:,sum(speaker==0)==0);
                
                for iter=1:iters;
                    temp=phase_rand([data  squeeze(gdata)],1);
                    speaker_perm=temp(:,1);
                    gdata_perm(1,:,:)=temp(:,2:end);
                    
                    speaker_perm=buffer(speaker_perm,win_width,win_width-win_step);
                    speaker_perm=speaker_perm(:,win_width:end);
                    
                    for s=1:size(gdata_perm,3);
                        
                        self_perm=squeeze(gdata_perm(:,:,s))';
                        self=squeeze(gdata(:,:,s))';
                        othersi=1:size(gdata,3);
                        othersi(othersi==s)=[];
                        others=nanmean(gdata(:,:,othersi),3)';
                        others_perm=nanmean(gdata_perm(:,:,othersi),3)';
                        
                        self_perm=buffer(self_perm,win_width,win_width-win_step);
                        self_perm=self_perm(:,sum(self_perm==0)==0);
                        self=buffer(self,win_width,win_width-win_step);
                        self=self(:,win_width:end);
                      
                        others=buffer(others,win_width,win_width-win_step);
                        others=others(:,win_width:end);
                   
                        others_perm=buffer(others_perm,win_width,win_width-win_step);
                        others_perm=others_perm(:,win_width:end);
                        
                        isc_SL_permperm(:,:,s,iter)=lagcorr_claire(others_perm,speaker_perm,lags);
                        isc_SL(:,:,s,iter)=lagcorr_claire(others_perm,speaker,lags);
                        isc_LL(:,:,s,iter)=lagcorr_claire(others_perm,self,lags);
                    end
                end
                isc_perm=0.5*log((1+isc_SL)./(1-isc_SL));
                save([expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/iscz_SL_' rname '_perm'],'isc_perm','lags','-v7.3');
                isc_perm=0.5*log((1+isc_LL)./(1-isc_LL));
                save([expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/iscz_LL_' rname  '_perm'],'isc_perm','lags','-v7.3');
                isc_perm=0.5*log((1+isc_SL_permperm)./(1-isc_SL_permperm));
                save([expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/iscz_SL_' rname  '_permperm'],'isc_perm','lags','-v7.3');
            end
        end
    end
    clear isc_SL isc_LL isc_SL_permperm gdata_perm
end

toc
