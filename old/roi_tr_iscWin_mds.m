
clear all;
tic
%loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

win_width=25;
win_step=5;

for ei=2%1:2;
    exp=experiments{ei};
    
    for ri=64;%1:length(rnames);
        rname=rnames{ri};
        fl=[expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname  ];
        
        if exist([fl '.mat'])>0;
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname ],'data');
            
            if size(gdata,1)>10 & size(data,1)>10;
                data=zscore(mean(data,1),0,2)';
                speaker=data;
                
                
                gdata=(zscore(nanmean(gdata,1),0,2));
                
                speaker_temp=[];
                gdata_temp=[];
                timeN=(size(gdata,2)-win_width+1);
                subjN=size(gdata,3);
                for ti=1:win_step:timeN
                    win_s=ti;
                    win_e=ti+win_width-1;
                    speaker_temp=[ speaker_temp; speaker(win_s:win_e)'];
                    gdata_temp=[ gdata_temp; squeeze(gdata(1,win_s:win_e,:))'];
                end
                
                temp=[speaker_temp ; gdata_temp];
                
                d=pdist(temp,'euclidean');
                mds_dim=2;
                [Y,stress,disparities]=mdscale(double(d),mds_dim);
                
                sampleN=(ti-1)/win_step+1;
                speaker_cord=Y(1:sampleN,:);
                listener_cord=Y((sampleN+1):end,:);
                listener_cord=reshape(listener_cord',mds_dim,subjN,sampleN);
                listener_cord=permute(listener_cord,[3 1 2]);
                
                clear F
                fn=1;
                
                for samplei=1:sampleN
                    
                    fsize=[10 10];
                    f=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
                    
                    
                    
                    scatter(listener_cord(samplei,1,:),listener_cord(samplei,2,:),40,'k','filled');
                    %scatter3(listener_cord(samplei,1,:),listener_cord(samplei,2,:),listener_cord(samplei,3,:),40,'k','filled');
                    hold on
                    scatter(speaker_cord(samplei,1),speaker_cord(samplei,2),60,'r','filled');
                    %scatter3(speaker_cord(samplei,1),speaker_cord(samplei,2),speaker_cord(samplei,3),1000,'r','filled','MarkerFaceAlpha',0.3);
                    
                    hold off
                    
                    
                    xlim([min(Y(:,1)) max(Y(:,1))]);
                    ylim([min(Y(:,2)) max(Y(:,2))]);
                    %   xlim([-10 10]);
                    %    ylim([-10 10])
                    F(fn) = getframe(f);
                    fn=fn+1;
                    
                    close gcf
                end
                figure('unit','centimeter','position',[2 2 fsize],'paperposition',[2 2 fsize],'papersize',fsize);
                %  movie(F,1,10)
                %   save( [expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/' rname '_mds.mat'],'Y','stress','disparities','speaker_cord','listener_cord');
                clear data data2 labels labels3
                %  close all
            end
        end
        
    end
end

