function [muw, sigma2w, Es, sigma2q]=vb(y,X, sigma2q)

% Ew=X\y;
Ew=0.1*randn(size(X,2),1);
% Ew=ridge(y,X,2);
% Ew=vbRVM(y,X,100)+0.1*randn(size(X,2),1);
sigma2w=var(Ew(Ew~=0));
% sigma2q=var(y-X*w1);

Es=rand(size(Ew));
iter=10000;
pi=0.1;
F=zeros(iter,1);

for i=1:iter
    Es=qs(y, X, Es, Ew, sigma2w, sigma2q, pi);
    [muw, sigma2qw]=qw(y, X, Es, Ew, sigma2w, sigma2q);
    sigma2w=m_sigma2w(muw,sigma2qw,Es);
%     sigma2q=m_sigma2w(mu,sigma2w,Es);
    pi=m_pi(Es);
    F(i)=lnlb(y, X, Es,muw, sigma2qw, sigma2w, sigma2q, pi);
end
figure; plot(F);
    
function Es=qs(y, X, Es, Ew, sigma2w, sigma2q,pi)
Esw=Es.*Ew;
XEsw=bsxfun(@times,X,Esw');
XEsw=bsxfun(@minus,sum(XEsw,2),XEsw);
error_m=bsxfun(@minus,y,XEsw);
s_ratio=sigma2q/sigma2w;
proj_x=sum(X.^2)'+s_ratio;
u=log(pi/(1-pi))+0.5*(log(s_ratio)-log(proj_x)...
    +sum(X.*error_m)'.^2)./(sigma2q*proj_x);
Es=1./(1+exp(-u));

function [muw, sigma2w]=qw(y, X, Es, Ew, sigma2w, sigma2q)
Esw=Es.*Ew;
XEsw=bsxfun(@times,X,Esw');
XEsw=bsxfun(@minus,sum(XEsw,2),XEsw);
error_m=bsxfun(@minus,y,XEsw);
s_ratio=sigma2q/sigma2w;
proj_x=sum(X.^2)'+s_ratio;
muw=sum(X.*error_m)'./proj_x;
sigma2w=sigma2q./proj_x;

function sigma2w=m_sigma2w(muw,sigma2qw,Es)
sigma2w=sum(Es.*(muw.^2+sigma2qw))/sum(Es);

function pi=m_pi(Es)
pi=mean(Es(:));