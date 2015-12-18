function [muw, sigma2qw, Es, sigma2q, F]=vb3(y,X, iter)

[N,M]=size(X);
% Ew=X\y;
% Ew=ridge(y,X,2);
% Ew=vbRVM(y,X,100)+0.1*randn(size(X,2),1);

% sigma2q=var(y-X*w1);

Es = 0.5*ones(M,1) + randn(M,1); 
Es(Es<=0.2)=0.2;
Es(Es>=0.8)=0.8;
% INITIALIZE FACTOR q(w_m | s_m =1) = N(w_m| muW_m,  sigma2_wm)
muw = 0.1*randn(M,1); 

sigma2w=100;
sigma2q=var(y-X*(Es.*muw));    
% iter=10000;
pi=0.5;
F=zeros(iter,1);
XtXm=sum(X.^2)';
yTy=y'*y;
XTy=X'*y;
XTX=X'*X;
for i=1:iter
    s_ratio=sigma2q/sigma2w;
    proj_x=XtXm+s_ratio;
    t1=log(pi/(1-pi))+0.5*(log(s_ratio)-log(proj_x));
    
    for j=1:M
        Esw=Es.*muw;
        XEsw=bsxfun(@times,X,Esw');
        XEsw=sum(XEsw,2)-XEsw(:,j);
        error_m=y-XEsw;
        proj_error=X(:,j)'*error_m;
        u=t1(j)+0.5*(proj_error^2)/(sigma2q*proj_x(j));
        Es(j)=1./(1+exp(-u));
        muw(j)=proj_error/proj_x(j);        
    end
    sigma2qw=sigma2q./proj_x;
    
    sigma2w=m_sigma2w(muw,sigma2qw,Es);
    sigma2q=m_sigma2q(yTy, XTX, XTy, Es,muw, sigma2qw,N);
    pi=m_pi(Es);    
    
%     F(i)=lnlb(y, X, Es,muw, sigma2qw, sigma2w, sigma2q, pi);
end
% figure; plot(F);
F=lnlb(y, X, Es,muw, sigma2qw, sigma2w, sigma2q, pi);

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