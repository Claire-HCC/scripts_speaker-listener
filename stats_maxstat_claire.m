%% Set paths
clearvars;
loc='mypc';
set_parameters;
input_types= {'trans_filtered_func' , 'trans_res4d'};
input_type=input_types{1};

iterations = 1000;
iterations_PerPiece = 20;
iterations_piece = params.iterations/params.iterations_PerPiece;
%% Max stat procedure
maxstat = NaN(1, iterations);

for ci = 1%:length(conditions);
    cond=conditions{ci};
     maxstat=[];
    for i= 1:iterations_piece;
        load(sprintf([expdir slash 'intersubj' slash input_type slash 'corrmap' slash 'Pair7_' cond '_corr_iter' num2str(iterations) '_%03d.mat'],i));
        maxstat((end+1):(end+iterations_PerPiece)) = max(bootstrap,[],2);
    end
    
    % percentile
    crit_upp = prctile(maxstat, 95);
    fprintf('%s 95 percentile threshold = %1.3f',cond,crit_upp);
    
    %% plot%
    figure('color', 'w', 'name', 'max bootstrapped r'); % , 'position', [2000 600 1330 300]);
    % plot
    hold on;
    [n, x] = hist(maxstat, 30);
    bar(x, n/sum(n));
    % plot crit
    yrange = get(gca, 'ylim');
    plot([crit_upp crit_upp], yrange, 'r', 'linewidth', 2);
    
    text(crit_upp, yrange(2)/2, ['r=' num2str(crit_upp)]);
    
    xlabel('max correlation');
    title(['Pair7 ' cond]);
    % set(gca, 'ylim', [0 .15], 'ytick', [0:.05:.15], 'xlim', [.05 .15], 'xtick', [.05:.025:.15]);
    box off;
    
end
