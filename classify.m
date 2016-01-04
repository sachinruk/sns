function [muw, sigma2qw, Es]=classify(z,X, iter)

[N,M]=size(X);
% Ew=X\Ey;
% Ew=ridge(Ey,X,2);
% Ew=vbRVM(Ey,X,100)+0.1*randn(size(X,2),1);

% 1=var(Ey-X*w1);

Es = 0.5*ones(M,1) + randn(M,1); 
Es(Es<=0.2)=0.2;
Es(Es>=0.8)=0.8;
% INITIALIZE FACTOR q(w_m | s_m =1) = N(w_m| muW_m,  sigma2_wm)
muw = 0.1*randn(M,1); 
% Ey=0.5*(2*z-1);

sigma2w=100;
% 1=var(Ey-X*(Es.*muw));       
% iter=10000;
pi=0.5;
% F=zeros(iter,1);
XtXm=sum(X.^2)';
% yTy=Ey'*Ey;
% XTy=X'*Ey;
% XTX=X'*X;
for i=1:iter
    s_ratio=1/sigma2w;
    proj_x=XtXm+s_ratio;
    t1=log(pi/(1-pi))+0.5*(log(s_ratio)-log(proj_x));
    
    for j=1:M
        Esw=Es.*muw;
        XEsw=bsxfun(@times,X,Esw');
        XEsw=sum(XEsw,2)-XEsw(:,j);
        [Ey]=qy(z,X*Esw);
        error_m=Ey-XEsw;
        proj_error=X(:,j)'*error_m;
        u=t1(j)+0.5*(proj_error^2)/(1*proj_x(j));
        Es(j)=1./(1+exp(-u));
        muw(j)=proj_error/proj_x(j);        
    end
    sigma2qw=1./proj_x;
    
    sigma2w=m_sigma2w(muw,sigma2qw,Es);
%     1=m_sigma2q(yTy, XTX, XTy, Es,muw, sigma2qw,N);
    pi=m_pi(Es);    
    
%     F(i)=lnlb(Ey, X, Es,muw, sigma2qw, sigma2w, 1, pi);
end
% figure; plot(F);
% F=lnlb(Ey, X, Es,muw, sigma2qw, sigma2w, 1, pi);

function sigma2w=m_sigma2w(muw,sigma2qw,Es)
sigma2w=sum(Es.*(muw.^2+sigma2qw))/sum(Es);

function pi=m_pi(Es)
pi=mean(Es(:));

function [Ey]=qy(z,mu)
C=sqrt(2*pi)*normcdf(-mu);
C(z==1)=sqrt(2*pi)-C(z==1);
Ey=mu+(2*z-1).*exp(-0.5*mu.^2)./C;

% function 1=m_sigma2q(yTy, XTX, XTy, Es,muw, sigma2qw,N)
% Esw=Es.*muw;
% Esw2=Es.*(muw.^2+sigma2qw);
% XTXEsw2=bsxfun(@times,bsxfun(@times,XTX,Esw),Esw');
% % XEsw=bsxfun(@times,X,Esw');
% % XEsw=bsxfun(@minus,sum(XEsw,2),XEsw);
% % error_m=Ey-2*sum(XEsw,2);
% % s_ratio=1/sigma2w;
% % proj_x=sum(X.^2)'+s_ratio;
% 1=(yTy-2*sum(XTy.*Esw)+sum(diag(XTX).*Esw2)+...
%             sum(XTXEsw2(:))-sum(diag(XTXEsw2)))/N;