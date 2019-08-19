clear all
% close all
loc='mypc';
set_parameters

iters=1000;
rnames={'HG_L','vPCUN','STC_L'};
role='listener';
exp='merlin';
textSim=repmat([1 -1; -1 1],11,11);

for ri=1:3;%:length(rnames);
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
        rs_pear(s)=corr(sim_temp(eye(22)==0),textSim(eye(22)==0),'type','spearman');
    end
    r_pear=mean(rs_pear);
    
    for iter = 1:iters;
        shuff1 = randperm(22);
        shuff2 = randperm(22);
        for s=1:18;
            train_subjects=1:18;
            train_subjects=train_subjects(train_subjects~=s);
            train_data=mean(data(:,:,train_subjects),3);
            test_data=data(:,:,s);
            sim_temp=corr(train_data,test_data);
            sim_temp=sim_temp(shuff1, shuff1);
            shuffled_mriSim(:,:,s) = sim_temp;
            r_null_temp(s)=corr(shuffled_mriSim(eye(22)==0), textSim(eye(22)==0), 'type', 'spearman');
            % Correlate
            
        end
        rs_pear_null(iter,:) = (r_null_temp);
        clear r_null_temp
    end
    r_pear_null=mean(rs_pear_null,2);
    % calculate p-values
    p_pear = sum( r_pear_null>=r_pear)/iters;
    
    disp(sprintf('%s Checkerboard: r=%0.3f (p=%0.3f)',rname,r_pear,p_pear));
end
%  clear r_pear r_spear r_pear_null r_spear_null p_pear p_spear p_thr_pear p_thr_spear

