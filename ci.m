function [CI]=ci(x,ciinv)
% comput ci range for each column of x

SEM = std(x)/sqrt(size(x,1));               % Standard Error
ts = tinv([ (1-(ciinv))/2  1-((1-(ciinv))/2) ],size(x,1)-1);      % T-Score
CI = ts'*SEM;
CI=CI(2,:)-CI(1,:);
end