% close all
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};

lags_tested={-15:15};

unsig_mask=imread([expdir 'scripts_speaker-listener/unsig_mask.tif']);
exps={'sherlock','merlin','pieman_old','black','forgot','21st_year','bronx','pieman','pieman_oldWord','pieman_rest','crossExps'};
eis=[1:11 ];

fs=[];
for ei=1:2;%
    exp=exps{eis(ei)};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks' ],...
            'lags','keptT','p_peaks','peakLags','peaks','pfdr_peaks','z_peaks','networks','p','pfdr_npeaks','npeaks');
        
        sig=(pfdr_peaks<(.05) & ~isnan(pfdr_peaks) & (abs(npeaks)<peaks | isnan(npeaks) ));
        peakLags_pfdr=peakLags;
        peakLags_pfdr(sig==0)=NaN;
        peakLags_pfdr
        
        
        fsize=[7 6.8];
        figure('unit','centimeter','position',[0 0 fsize],'InvertHardcopy', 'off','color','w');
        
        
        peakLags_bined=discretize(peakLags,[ -15 -10 -6 -3 0 1 4  7  11 16]);
        
        % to mask unsignificant peakLags
        im2=imagesc(peakLags_bined,[1 9]);
        colormap jet
        imAlpha=(sig==0)*0.35;
        set(im2,'AlphaData',imAlpha);
        lm=get(gca,'xlim');
        hold on
        h = image([lm(1)-1 lm(2)+1],[lm(1)-1 lm(2)+1],unsig_mask);
        set(h,'AlphaData',max(unsig_mask,[],3)<=85)
        %  uistack(h,'top');
        
        % the significant ones
        im=imagesc(peakLags_bined,[1 9]);
        colormap jet
        imAlpha=ones(size(peakLags));
        imAlpha(sig==0)=0;
        set(im,'AlphaData',imAlpha)
        
        
        for ni=1:length(networks);
            network=networks{ni};
            
            pos=[ni-0.5 ni-0.5 1 1];
            rectangle('Position',pos ,'edgecolor','w','linewidth',1.5);
            text(pos(1)+pos(3)/2,pos(2)+pos(4)*0.5,network,'HorizontalAlignment','center','color','w','fontweight','bold','fontsize',8);
        end
        hold off
        exp=strrep(exp,'pieman_oldWord','Scrambled Pieman');
        exp=strrep(exp,'pieman_old','Old Pieman');
        exp=strrep(exp,'pieman_rest','Resting');
        
        title(strrep([upper(exp(1)) exp(2:end)],'_',' '));
        
        set(gca,'ytick',[],'xtick',[],'xlim',lm,'ylim',lm,'color','k');
        
        %% without tick label
        yl=ylabel('Seed network');
        %    yl.Position(1)=7.3
        xlabel('Target network');
        set(gca,'fontsize',12)
        
        
        %% with tick label
        %   fsize=[7.5 9.5];
        % set(gcf,'position',[0 0 fsize]);
        %         yyaxis right
        %         yyaxis left
        %         ax=gca;
        %         ax.YAxis(2).Color = 'k';
        %         %  box off
        %         set(ax,'xtick',1:length(networks),'ytick',1:length(networks),....
        %             'xticklabels',cellfun(@(x) [x ' target'],networks,'UniformOutput',0),...
        %             'yticklabels',cellfun(@(x) [x ' seed'],networks,'UniformOutput',0),...
        %             'color','k','xlim',lm,'ylim',lm);
        %         xtickangle(75);
        %         set(gca,'fontsize',10)
        
        print(gcf,['temp' num2str(ei) '.tif'],'-dtiff','-r500');
        fs=[fs  imread(['temp' num2str(ei) '.tif'])];
    end
end

imshow(fs)
imwrite(fs,'temp.tiff')


%% color bar
% cb=interp1(-15:15,smooth(discretize(-15:15,[ -15 -10 -6 -3 0 1 4  7  11 16]),10),-15:0.01:15);
% imagesc(cb); colormap jet
% set(gca,'xtick',[],'ytick',[])