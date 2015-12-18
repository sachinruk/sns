function [mu,Sigma,A]=classify_RVM(z,X,iter)

sigma2=0.5; lambda=1e-9; c=1e-9; d=1e-9;
[~,D]=size(X);
Phi2=X'*X;
A=ones(D,1); 
mu=X\z+0.5*randn(D,1);
% l=1;
for i=1:iter
    [Ey]=qy(z,X*mu);
    [Sigma, mu]=weightEst(Ey,X,sigma2,A);
    A=alphaEst(Sigma,mu,lambda);
    sigma2=Esigma(Ey,X,mu,Sigma,Phi2, c,d,D);
    
end

function sigma2=Esigma(y,X,mu,Sigma,Phi2, c,d,n)
a=c+n/2;
PhiTPhiSigma=Phi2.*Sigma;
b=0.5*(norm(y-X*mu)+sum(PhiTPhiSigma(:))+2*d);

sigma2=b/a;

function [Sigma, mu]=weightEst(y,Phi,sigma2,A)
[m,~]=size(Phi);
Ainv=diag(1./A);
% Sigma=Ainv-Ainv*Phi'*((sigma2*eye(m)+Phi*Ainv*Phi')\Phi*Ainv);
% mu=Sigma*(Phi'*y)/sigma2;
d=size(Phi,2);
SigmaInv=Phi'*Phi; SigmaInv(1:(d+1):end)=SigmaInv(1:(d+1):end)+A./sigma2;
L=chol(SigmaInv,'lower');
mu=L'\(L\(Phi'*y));


function A=alphaEst(Sigma,mu,lambda)
A=3./(diag(Sigma)+mu.^2+lambda);

function [Ey]=qy(z,mu)
C=sqrt(2*pi)*normcdf(-mu);
C(z==1)=sqrt(2*pi)-C(z==1);
Ey=mu;
Ey=Ey+(2*z-1).*exp(-0.5*mu.^2)./C;