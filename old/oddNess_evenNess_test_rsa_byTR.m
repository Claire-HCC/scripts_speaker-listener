clear all
% close all
loc='mypc';
set_parameters

iters=1000;
rnames={'HG_L','vPCUN','STC_L'};
role='listener';
exp='sherlock';

timeUnit='tr';

load([expdir exp '\sound\listenerEvents.mat']);
textV=-ones(1,length(listenerEvents_trVector));
textV(mod(listenerEvents_trVector,2)==1)=1;
textV=textV(listenerEvents_trVector>0);
textSim=(textV')*(textV);
eventn=size(textSim,2);

for ri=3%1:3;%:length(rnames);
    rname=rnames{ri};
    load([expdir exp '/fmri/mat/roi/' timeUnit '/' timeUnit '_' role '_' rname '.mat']);
    
    % delete pause scane
        data=zscore(data,0,2);
    data=data(:,listenerEvents_trVector>0,:);

    
    for s=1:18;
        train_subjects=1:18;
        train_subjects=train_subjects(train_subjects~=s);
        train_data=mean(data(:,:,train_subjects),3);
        test_data=data(:,:,s);
        sim_temp=corr(train_data,test_data);
        sim(:,:,s)=sim_temp;
        rs_pear(s)=corr(sim_temp(eye(eventn)==0),textSim(eye(eventn)==0));
    end
    r_pear=mean(rs_pear);
    
    for iter = 1:iters;
        shuff1 = randperm(eventn);
        shuff2 = randperm(eventn);
        for s=1:18;
            train_subjects=1:18;
            train_subjects=train_subjects(train_subjects~=s);
            train_data=mean(data(:,:,train_subjects),3);
            test_data=data(:,:,s);
            sim_temp=corr(train_data,test_data);
            sim_temp=sim_temp(shuff1, shuff1);
            shuffled_mriSim(:,:,s) = sim_temp; 
            r_null_temp(s)=corr(sim_temp(eye(eventn)==0), textSim(eye(eventn)==0), 'type', 'pearson');
            % Correlate
            
        end
        rs_pear_null(iter,:) = (r_null_temp);
        clear r_null_temp
    end
     r_pear_null=mean( rs_pear_null,2);
    % calculate p-values
    p_pear = sum( r_pear_null>=r_pear)/iters;
    
    disp(sprintf('%s Checkerboard: r=%0.3f (p=%0.3f)',rname,r_pear,p_pear));
end
%  clear r_pear r_spear r_pear_null r_spear_null p_pear p_spear p_thr_pear p_thr_spear

