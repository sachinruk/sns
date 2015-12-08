function [muw, sigma2w, Es, sigma2q]=vb(y,X, sigma2q)

Ew=X\y;
sigma2w=var(Ew(Ew~=0));
% sigma2q=var(y-X*w1);

Es=rand(size(Ew));
iter=100;

for i=1:iter
    Es=qs(y, X, Es, Ew, sigma2w, sigma2q);
    [muw, sigma2qw]=qw(y, X, Es, Ew, sigma2w, sigma2q);
    sigma2w=m_sigma2w(muw,sigma2qw,Es);
%     sigma2q=m_sigma2w(mu,sigma2w,Es);
    pi=m_pi(Es);
end
    
function Es=qs(y, X, Es, Ew, sigma2w, sigma2q,pi)
Esw=Es.*Ew;
XEsw=bsxfun(@times,X,Esw');
XEsw=bsxfun(@minus,sum(XEsw,2),XEsw);
error_m=bsxfun(@minus,y,XEsw);
s_ratio=sigma2q/sigma2w;
proj_x=sum(X.^2)'+s_ratio;
Es=log(pi/(1-pi))+0.5*(log(s_ratio)-log(proj_x)...
    +sum(X.*error_m)'.^2)./(sigma2q*proj_x);

function [muw, sigma2w]=qw(y, X, Es, Ew, sigma2w, sigma2q)
Esw=Es.*Ew;
XEsw=bsxfun(@times,X,Esw');
XEsw=bsxfun(@minus,sum(XEsw,2),XEsw);
error_m=bsxfun(@minus,y,XEsw);
s_ratio=sigma2q/sigma2w;
proj_x=sum(X.^2)'+s_ratio;
muw=X.*error_m./proj_x;
sigma2w=sigma2q./proj_x;

function sigma2w=m_sigma2w(muw,sigma2qw,Es)
sigma2w=sum(Es.*(muw.^2+sigma2qw))/sum(Es);

function pi=m_pi(Es)
pi=mean(Es(:));