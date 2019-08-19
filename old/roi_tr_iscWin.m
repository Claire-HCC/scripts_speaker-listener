
clear all;
tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
 rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

tic % 15 min
win_width=25;
win_step=1;
lags=-20:20;

for ei=2%1:2;
    exp=experiments{ei};
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        fl=[expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname  ];
        
        if exist([fl '.mat'])>0;
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname ],'data');
            
            if size(gdata,1)>10 & size(data,1)>10;
                data=zscore(mean(data,1),0,2)';
                speaker=data;
                speaker=buffer(speaker,win_width,win_width-win_step);
                speaker=speaker(:,sum(speaker==0)==0);
                
                gdata=(zscore(nanmean(gdata,1),0,2));
                
                for s=1:size(gdata,3);
                    
                    self=squeeze(gdata(:,:,s))';
                    othersi=1:size(gdata,3);
                    othersi(othersi==s)=[];
                    others=nanmean(gdata(:,:,othersi),3)';
                    
                    
                    self=buffer(self,win_width,win_width-win_step);
                    self=self(:,win_width:end);
                    others=buffer(others,win_width,win_width-win_step);
                    others=others(:,sum(others==0)==0);
                    
                    isc_SL(:,:,s)=lagcorr_claire(self,speaker,lags);
                    isc_LL(:,:,s)=lagcorr_claire(self,others,lags);
                end
                
                isc=0.5*log((1+isc_SL)./(1-isc_SL));
                save([expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/iscz_SL_' rname ],'isc','lags');
                isc=0.5*log((1+isc_LL)./(1-isc_LL));
                save([expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/iscz_LL_' rname  ],'isc','lags');
                
            end
        end
    end
        clear isc_SL isc_LL
end

toc
