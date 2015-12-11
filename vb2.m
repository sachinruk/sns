function [muw, sigma2w, Es, sigma2q, F]=vb2(y,X, sigma2q,iter)
%vb but using Hessians
[N,M]=size(X);
% Ew=X\y;
muw=0.1*randn(size(X,2),1);
% Ew=ridge(y,X,2);
% Ew=vbRVM(y,X,100)+0.1*randn(size(X,2),1);
sigma2w=var(muw(muw~=0));
sigma2qw=zeros(size(muw));
% sigma2q=var(y-X*w1);

Es=rand(size(muw));
% iter=1000;
pi=0.5;
F=zeros(iter,1);

yTy=y'*y;
XTy=X'*y;
XTX=X'*X;

for i=1:iter
    Es=qs(y, X, Es, muw, sigma2w, sigma2q, pi);
    [muw, sigma2qw]=qw(y, X, Es, sigma2w, sigma2q,muw, sigma2qw);
    sigma2w=m_sigma2w(muw,sigma2qw,Es);
    sigma2q=m_sigma2q(yTy, XTX, XTy, Es,muw, sigma2qw,N);
    pi=m_pi(Es);
    F(i)=lnlb(y, X, Es,muw, sigma2qw, sigma2w, sigma2q, pi);
end
% figure; plot(F);
if sum(diff(F)<0)
    error('going down');
end
F=lnlb(y, X, Es,muw, sigma2qw, sigma2w, sigma2q, pi);

    
function Es=qs(y, X, Es, Ew, sigma2w, sigma2q,pi)
Esw=Es.*Ew;
XEsw=bsxfun(@times,X,Esw');
XEsw=bsxfun(@minus,sum(XEsw,2),XEsw);
error_m=bsxfun(@minus,y,XEsw);
s_ratio=sigma2q/sigma2w;
proj_x=sum(X.^2)'+s_ratio;
u=log(pi/(1-pi))+0.5*(log(s_ratio)-log(proj_x))...
    +0.5*(sum(X.*error_m)'.^2)./(sigma2q*proj_x);
% u=log(pi/(1-pi))+0.5*(log(s_ratio)-log(proj_x)...
%     +sum(X.*error_m)'.^2)./(sigma2q*proj_x);
Es=1./(1+exp(-u));

function [muw, sigma2qw]=qw(y, X, Es, sigma2w, sigma2q,muw, sigma2qw)
idx=Es>0; X=X(:,idx); Es=Es(idx);
Esw=Es.*muw(idx);
XEsw=bsxfun(@times,X,Esw');
XEsw=bsxfun(@minus,sum(XEsw,2),XEsw);
error_m=bsxfun(@minus,y,XEsw);
s_ratio=sigma2q/sigma2w;
proj_x=sum(X.^2)'+s_ratio;
% Es_sigma2q=Es./sigma2q;

g=(sum(X.*error_m)'-proj_x.*muw(idx)).*Es;
H=bsxfun(@times,bsxfun(@times,X'*X,Es),Es');
H(1:(length(H)+1):end)=proj_x.*Es;
% L=jitChol(H);
L=chol(H,'lower');
muw(idx)=muw(idx)+L'\(L\g);
sigma2qw(idx)=sigma2q./proj_x;

function sigma2w=m_sigma2w(muw,sigma2qw,Es)
sigma2w=sum(Es.*(muw.^2+sigma2qw))/sum(Es);

function pi=m_pi(Es)
pi=mean(Es(:));

function sigma2q=m_sigma2q(yTy, XTX, XTy, Es,muw, sigma2qw,N)
Esw=Es.*muw;
Esw2=Es.*(muw.^2+sigma2qw);
XTXEsw2=bsxfun(@times,bsxfun(@times,XTX,Esw),Esw');
% XEsw=bsxfun(@times,X,Esw');
% XEsw=bsxfun(@minus,sum(XEsw,2),XEsw);
% error_m=y-2*sum(XEsw,2);
% s_ratio=sigma2q/sigma2w;
% proj_x=sum(X.^2)'+s_ratio;
sigma2q=(yTy-2*sum(XTy.*Esw)+sum(diag(XTX).*Esw2)+...
            sum(XTXEsw2(:))-sum(diag(XTXEsw2)))/N;