%б@всврв█бH
clear all;
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

binSize=10;% tr;
binStep=1;

type='corr';
col_speaker='r'; %'r'

for ei=3;%1:2;
    exp=experiments{ei};
    
    ris=find(ismember(rnames,{'dPCC'}));
    % [~,ri]=max(herdm);
    
    for rii=1:length(ris);
        ri=ris(rii);
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/' rname '_tr_pattern_binSize' num2str(binSize) '_metric_mds_' type '_' col_speaker '_2d.mat'],'speaker_cord','listener_cord','binN','binSize','binStep','rname')
        
        lmax=max(max(abs(repmat(speaker_cord,1,1,18)-listener_cord),[],3));
        clear F
        fn=1;
        
        ds=diff(speaker_cord);
        dl=diff(listener_cord,1,1);
        for ti=10:1:binN;
            
            fsize=[10 10];
            f=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
            hold on
            
            
            scatter(speaker_cord(ti,1),speaker_cord(ti,2),40,col_speaker,'filled');
            scatter(listener_cord(ti,1,:),listener_cord(ti,2,:),40,'k','filled');
            
            
            if ti>1;
                
                quiver(speaker_cord(ti-1,1),speaker_cord(ti-1,2),ds(ti-1,1),ds(ti-1,2),0,'color',col_speaker);
                % quiver(listener_cord(ti-1,1,:),listener_cord(ti-1,2,:),dl(ti-1,1,:),dl(ti-1,2,:),0,'k');
            end
            if ti>2;
                ang=rad2deg(angle(ds(ti-2,1)+j*ds(ti-1,2)));
                view(ang+270,90);
            else
                view(0,90);
            end
            
            
            hold off
            
            xlim(speaker_cord(ti,1)+[-lmax(1) lmax(1)]);
            ylim(speaker_cord(ti,2)+[-lmax(2) lmax(2)]);
            axis off
            
            F(fn) = getframe(f);
            fn=fn+1
            
            close gcf
        end
        beep
        
        myVideo = VideoWriter([expdir '/' exp '/fmri/' rname '_tr_pattern_metric_binSize' num2str(binSize) '_mds_' type '_' col_speaker '_focused_2d.avi']);
        myVideo.FrameRate = binSize/2;  % Default 30
        open(myVideo);
        writeVideo(myVideo,F);
        close(myVideo)
        
        clear data data2 labels labels3
        %  close all
    end
end

