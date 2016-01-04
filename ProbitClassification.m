function [mu, SigmaInv,C,Ey,Ey2]=ProbitClassification(t,X,iter)
mu=X\t;
% iterations=100;
% mu_=zeros(length(mu),iterations);
% i=1;
% diff=1;

%Sigma is set no need to change
SigmaInv=X'*X; step=length(SigmaInv);
SigmaInv(1:(step+1):end)=SigmaInv(1:(step+1):end)+1/1000;
SigmaX=SigmaInv\X';
for i=1:iter
    [C, Ey, Ey2]=q_y(X,t,mu);
    mu=q_beta(SigmaX,Ey);
%     mu_(:,i)=mu;
%     if i>1, diff=norm(mu-mu_(:,i-1)); end 
%     i=i+1;
end
% figure; plot(mu_(1,:)); hold on; plot(mu_(2,:),'r'); hold off; disp(i)

function mu=q_beta(SigmaX,Ey)

mu=SigmaX*Ey;

function [C, Ey, Ey2]=q_y(X,t,mu)

mu2=X*mu;
C=sqrt(2*pi)*normcdf(-mu2);
C(t==1)=sqrt(2*pi)-C(t==1);

Ey=(-1).^(1-t).*exp(-0.5*mu2.^2)./C+mu2;
Ey2=(-1).^(1-t).*mu2.*exp(-0.5*mu2.^2)./C+mu2.^2+1;