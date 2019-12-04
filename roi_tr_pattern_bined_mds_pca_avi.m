%¡@score*coeff'+mu=input
clear all;
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

binSize=20;% tr;
binStep=1;
type='corr';
col_speaker='r'; %'r'

for ei=3%1:2;
    exp=experiments{ei};
    
    ris=find(ismember(rnames,{'dPCC'}));
    % [~,ri]=max(herdm);
    
    for rii=1:length(ris);
        ri=ris(rii);
        rname=rnames{ri};
        load([expdir '/' exp '/fmri/' rname '_tr_pattern_binSize' num2str(binSize) '_metric_mds_' type '_' col_speaker '_3d.mat'],'speaker_cord','listener_cord','binN','binSize','binStep','rname');
        lmax=max(squeeze(max([listener_cord  ]))');
        lmin=min(squeeze(min([listener_cord  ]))');
        
        clear F
        fn=1;
        for ti=1:1:binN;
            
            fsize=[10 10];
            f=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
            
            scatter3(listener_cord(ti,1,:),listener_cord(ti,2,:),listener_cord(ti,3,:),40,'k','filled');
            hold on
             scatter3(speaker_cord(ti,1),speaker_cord(ti,2),speaker_cord(ti,3),40,'r','filled');
            % scatter3(speaker_cord(ti,1),speaker_cord(ti,2),speaker_cord(ti,3),40,'k','filled');
            hold off
            
            %xlim([lmin(1) lmax(1)]);
            % ylim([lmin(2) lmax(2)]);
            axis off
            xlim([-1 1]);
            ylim([-1 1]);
            zlim([-1 1]);
            % view(90,0); % yz
             view(0,0); % xz
            % view(0,90); % xy
            F(fn) = getframe(f);
            fn=fn+1
            
            close gcf
        end
        beep
        myVideo = VideoWriter([expdir '/' exp '/fmri/' rname '_tr_pattern_binSize' num2str(binSize) '_metric_mds_' type '_' col_speaker '_3d_xz_ans.avi']);
        myVideo.FrameRate = binSize/2;  % Default 30
        open(myVideo);
        writeVideo(myVideo,F);
        close(myVideo)
        
        close all
    end
end


