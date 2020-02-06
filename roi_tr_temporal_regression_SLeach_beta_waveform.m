clear all
% close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
crop_start=10;
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};

rnames_selected={'dPCUN'};
figure;
for ei=[1:4];%[1:4];;%
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLeach/regression_SLeach_lag' num2str(min(lags)) '-' num2str(max(lags)) '_classification.mat'  ],'rnames','sig_fdr','b_s','b_l');
%        sig_betaClass=sig_fdr;
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            
            if ismember(rname,rnames_selected);
                %  if sig_betaClass(ri)==1;
                
                b_s_temp=squeeze(b_s(ri,:,:))';
                b_l_temp=squeeze(b_l(ri,:,:))';
                
                subplot(1,4,ei);
                ciplot_claire(squeeze(b_s_temp),lags,'r',0.3);
                hold on;
                ciplot_claire(b_l_temp,lags,'b',0.3)
                title({exp})
            %    line([bPeakLags(ri,1) bPeakLags(ri,1)],[get(gca,'ylim')],'color','k');
                
                grid on
                hold off
            
                if ei==1;
                       xlabel('Speaker precedes-----TR(1.5s)-----Listeners precede','fontsize',10); ylabel('Beta values'); 
                end
            end
            
        end
    end
end

