
clear all;
tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

win_width=30;
win_step=1;
type='correlation';

for ei=3%1:2;
    exp=experiments{ei};
    
    ris=find(ismember(rnames,{'aCUN'}));
    % [~,ri]=max(herdm);
    
    for rii=1:length(ris);
        ri=ris(rii);
        rname=rnames{ri};
        fl=[expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname  ];
        
        if exist([fl '.mat'])>0;
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname ],'data');
            
            if size(gdata,1)>10 & size(data,1)>10;
                timeN=size(gdata,2);
                subjN=size(gdata,3);
                speaker_temp=data';
                gdata_temp=permute(gdata,[3 2 1]);
                gdata_temp=reshape(gdata_temp,subjN*timeN,size(gdata,1));
                
                temp=[speaker_temp ; gdata_temp];
                
                d=pdist(temp,type);
                mds_dim=2;
                Y=cmdscale(double(d),mds_dim);
                
                speaker_cord=Y(1:timeN,:);
                listener_cord=Y((timeN+1):end,:);
                listener_cord=reshape(listener_cord',mds_dim,subjN,timeN);
                listener_cord=permute(listener_cord,[3 1 2]);
                
                lmax=max(squeeze(max([listener_cord  ]))');
                lmin=min(squeeze(min([listener_cord  ]))');
                clear F
                fn=1;
                
                for ti=1:1:timeN;
                    
                    fsize=[10 10];
                    f=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
                    
                    
                    scatter(listener_cord(ti,1,:),listener_cord(ti,2,:),40,'k','filled');
                    % scatter3(listener_cord(ti,1,:),listener_cord(ti,2,:),listener_cord(ti,3,:),40,'k','filled');
                    hold on
                    scatter(speaker_cord(ti,1),speaker_cord(ti,2),40,'r','filled');
                    % scatter3(speaker_cord(ti,1),speaker_cord(ti,2),speaker_cord(ti,3),40,'r','filled');
                    hold off
                    
                    xlim([lmin(1) lmax(1)]);
                    ylim([lmin(2) lmax(2)]);
                    axis off
                    %                     xlim([-0.8 0.8]);
                    %                     ylim([-0.8 0.8])
                    %                     zlim([-0.8 0.8])
                    F(fn) = getframe(f);
                    fn=fn+1
                    
                    close gcf
                end
                beep
                % pause
                % figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
                % movie(F,1,5)
                myVideo = VideoWriter([expdir '/' exp '/fmri/herding_' rname '_mds_' type '.avi']);
                myVideo.FrameRate = 10;  % Default 30
                open(myVideo);
                writeVideo(myVideo,F);
                close(myVideo)
                
                clear data data2 labels labels3
                %  close all
            end
        end
        
    end
end

