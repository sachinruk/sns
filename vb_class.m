function [muw, sigma2qw, Es, sigma2q, F]=vb_class(y,X, sigma2q, iter,K)

[D,M]=size(X);
N=size(y,2); % y is DxN
% Ew=X\y;
% Ew=ridge(y,X,2);
% Ew=vbRVM(y,X,100)+0.1*randn(size(X,2),1);

% sigma2q=var(y-X*w1);

Es = 0.5*ones(M,K) + 0.001*randn(M,K); 
Es(Es<=0.2)=0.2;
Es(Es>=0.8)=0.8;
% INITIALIZE FACTOR q(w_m | s_m =1) = N(w_m| muW_m,  sigma2_wm)
muw = 0.1*randn(M,N); 

sigma2w=2;
    
% iter=10000;
pi=0.5;
F=zeros(iter,1);

yTy=y'*y;
XTy=X'*y; %M x N
XTX=X'*X;
XtXm=diag(XTX);
for i=1:iter
    s_ratio=sigma2q/sigma2w;
    proj_x=XtXm+s_ratio;
    t1=log(pi./(1-pi))+0.5*N*(log(s_ratio));
    
    for k=1:K
        pi_k=pi_nk(:,k);
        pikXtxm=pi_k*XtXm';
    for m=1:M
        Esw=bsxfun(@times,Es(:,k),muw); %M x N
%         XEsw=bsxfun(@times,X,Esw');
%         XEsw=sum(XEsw,2)-XEsw(:,j);
%         error_m=y-XEsw;
%         proj_error=X(:,j)'*error_m;
        XTXEsw=bsxfun(@times,XTX(:,m),Esw); %MxN matrix
        proj_error=XTy(m,:)-sum(XTXEsw)+XTXEsw(m,:);
        u=t1(m)+sum(0.5*((pi_k.*proj_error)^2)./(sigma2q*(pikXtxm(:,m)+s_ratio))...
            -log(pikXtxm(:,m)+s_ratio));
        Es(m)=1./(1+exp(-u));
        muw(m)=proj_error/proj_x(m);        
    end
    end
    sigma2qw=sigma2q./proj_x;
    
    sigma2w=m_sigma2w(muw,sigma2qw,Es);
    sigma2q=m_sigma2q(yTy, XTX, XTy, Es,muw, sigma2qw,D);
    pi=m_pi(Es);    
    
%     F(i)=lnlb(y, X, Es,muw, sigma2qw, sigma2w, sigma2q, pi);
end
% figure; plot(F);
F=lnlb(y, X, Es,muw, sigma2qw, sigma2w, sigma2q, pi);
    
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