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

for ei=1;%1:2;
    exp=experiments{ei};
    
    ris=find(ismember(rnames,{'vPCUN'}));
    % [~,ri]=max(herdm);
    
    for rii=1:length(ris);
        ri=ris(rii);
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/animation/' rname '_tr_pattern_binSize' num2str(binSize) '_mds_' type  '_2d.mat'],'speaker_cord','listener_cord','binN','binSize','binStep','rname')
        
        lmax=max(max(abs(repmat(speaker_cord,1,1,size(listener_cord,3))-listener_cord),[],3));
        clear F
        fn=1;
        
        for ti=1:binStep:binN;
            
            fsize=[10 10];
            f=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
            hold on
            
            listenerM_cord=mean(listener_cord(ti,:,:),3);
            ds_sl=pdist([speaker_cord(ti,:); listenerM_cord],'euclidean');
            ds_ll=mean(pdist(squeeze(listener_cord(ti,:,:))','euclidean'));
            
            scatter(listener_cord(ti,1,:),listener_cord(ti,2,:),10,'b','filled');
            viscircles(listenerM_cord,ds_ll,'color','b')
            
            scatter(speaker_cord(ti,1),speaker_cord(ti,2),10,col_speaker,'filled');
            line([speaker_cord(ti,1) listenerM_cord(1)],[speaker_cord(ti,2) listenerM_cord(2)],'color','r','linewidth',2);
            
            v_sl=speaker_cord(ti,:)-listenerM_cord;
            ang=rad2deg(angle(v_sl(1)+j*v_sl(2)));
            view(ang,90);
            
            % focus on speaker
            %     xlim(speaker_cord(ti,1)+[-1.5*lmax(1) 1.5*lmax(1)]);
            %    ylim(speaker_cord(ti,2)+[-1.5*lmax(2) 1.5*lmax(2)]);
            
            % focus on listenerM
            %  xlim(listenerM_cord(1)+[-1.5*lmax(1) 1.5*lmax(1)]);
            %  ylim(listenerM_cord(2)+[-1.5*lmax(2) 1.5*lmax(2)]);
            xlim(listenerM_cord(1)+[-2 2]);
            ylim(listenerM_cord(2)+[-2 2]);
            
            axis square
            axis off
            hold off
            
            F(fn) = getframe(f);
            fn=fn+1
            
            close gcf
        end
        beep
        
        myVideo = VideoWriter([expdir '/' exp '/fmri/animation/' rname '_tr_pattern_binSize' num2str(binSize) '_mds_' type '_2d_focused.avi']);
        myVideo.FrameRate = binSize/2;  % Default 30
        myVideo.Quality=100;
        open(myVideo);
        writeVideo(myVideo,F);
        close(myVideo)
        
        clear data data2 labels labels3
        %  close all
    end
end

