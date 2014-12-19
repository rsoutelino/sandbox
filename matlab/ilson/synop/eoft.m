function [evalue,evect,amp]=eoft(trmat)

% function [evalue,evect,amp]=eoft(trmat)
% function [evalue,evect,amp]=eoft(trmat)
% computes eof in time (or depth) for general cases
% normally trmat is a matrix containing each time series 
% (or vertical profiles) as a row. 
% That is, for N stations having m data points, 
% trmat is a N by m matrix.


% demeans trmat prior to doing eof

[N m] = size(trmat);
trmatmean = mean(trmat')';
trmat = trmat - trmatmean*ones(1,m);

% computes zero-lag cross-covariance matrix

mcov = trmat(1:N,:)*trmat(1:N,:)';

% computes eigenvalues and eigenvectors

[evect evalue] = eig(mcov);

% normalize eigenvectors

for i = 1:N;
 evect(:,i) = evect(:,i)/norm(evect(:,i));
end

% sort eigenvalues and computes percent variance
% explained by each mode

[evalue, index] = sort(diag(evalue/trace(evalue)));


% computes amplitude functions 

amp = evect'*trmat(1:N,:);










