function [L, jitter]=jitChol(K, maxTries)

if nargin < 2
  maxTries = 5;
end
n=size(K,1); e=1e-10; jitter=0;
L=[];
for i=1:maxTries
    try
        L=chol(K,'lower');
    catch
        K(1:(n+1):end)=K(1:(n+1):end)+e;
        jitter=jitter+e;
        e=e*10;
        continue;
    end
    break;
end

% if isempty(L) %if nothing was assigned in previous step,
%     K(1:(n+1):end)=K(1:(n+1):end)-jitter;
%     e=1e-10; jitter=0; %reset parameters
%     l=max(diag(K)); K=K/l;
%     for i=1:maxTries
%         try
%             L=chol(K,'lower');
%         catch
%             K(1:(n+1):end)=K(1:(n+1):end)+e;
%             jitter=jitter+e;
%             e=e*10;
%             continue;
%         end
%         L=sqrt(l)*L;
%         break;
%     end
% end
%         