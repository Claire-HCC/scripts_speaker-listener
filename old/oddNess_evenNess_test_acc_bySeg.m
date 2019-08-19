clear all
% close all
loc='mypc';
set_parameters

iters=1000;
rnames={'HG_L','vPCUN','STC_L'};
role='listener';
exp='merlin';
textSim=repmat([1 -1; -1 1],11,11);

for ri=2%:3%1:length(rnames);
    rname=rnames{ri};
    load([expdir exp '/fmri/mat/roi/segment/segment_' role '_' rname '.mat']);
    
    data=zscore(data,0,2);
    for s=1:18;
        train_subjects=1:18;
        train_subjects=train_subjects(train_subjects~=s);
        train_data=mean(data(:,:,train_subjects),3);
        test_data=data(:,:,s);
        sim_temp=corr(train_data,test_data);
        
        sim(:,:,s)=sim_temp;
        segn=size(sim_temp,1);
        [~,rk]=sort(sim_temp,'descend');
        [~,rk]=sort(rk);
        rk(eye(segn)==1)=NaN;
        
        for segi=1:segn;
            ds(s,segi)=nanmean(rk(textSim(:,segi)==1,segi))-nanmean(rk(textSim(:,segi)==-1,segi));
        end
    end
    
    
    for iter = 1:iters;
        shuff1 = randperm(segn);
        
        for s =1:18;
            sim_temp= sim(shuff1, shuff1,s);
            
            shuffled_mriSim(:,:,iter)=sim_temp;
            
            segn=size(sim_temp,1);
            [~,rk]=sort(sim_temp,'descend');
            [~,rk]=sort(rk);
            rk(eye(segn)==1)=NaN;
            
            for segi=1:segn;
                ds_null(s,segi,iter)=nanmean(rk(textSim(:,segi)==1,segi))-nanmean(rk(textSim(:,segi)==-1,segi));
            end
        end
    end
    d=mean(ds(:));
    d_null=mean(mean(ds_null,1),2);
    % d_null=m_null(:,2)-m_null(:,1);
    
    % calculate p-values
    p= sum(d_null>=d)/iters;
    
    
    disp(sprintf('%s Checkerboard: d=%0.3f (p=%0.3f)',rname,d,p));
    %  clear r_pear r_spear r_pear_null r_spear_null p_pear p_spear p_thr_pear p_thr_spear
end
