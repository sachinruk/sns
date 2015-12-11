function [mu,Sigma,A]=vbRVM(y,Phi,iter)

sigma2=0.5; lambda=1e-3; c=1e-9; d=1e-9;
[~,n]=size(Phi);
Phi2=Phi'*Phi;
A=ones(n,1); 
% l=1;
for i=1:iter
    [Sigma, mu]=weightEst(y,Phi,sigma2,A);
%     if (i>1)
%         change_mu(i-1)=norm(mu-mu_old);        
%     end
%     mu_old=mu;
%     for j=1:25
%         lambda=lambdaEst(theta,k,A);
    A=alphaEst(Sigma,mu,lambda);
    sigma2=Esigma(y,Phi,mu,Sigma,Phi2, c,d,n);
%         if (i>1 || j>2)
%             delta_A(l)=norm(A-A_old);
%             delta_lambda(l)=norm(lambda-lambda_old);
%             l=l+1;
%         end
%         A_old=A;
%         lambda_old=lambda;
%     end
%     k1=diag(Sigma)+mu.^2;
%     k2=1+theta;
%     k3=k;
%     a=k1;
%     b=2*k1*k3+2*k2-3;
%     c=-6*k3;
%     delta=b.^2-4*a*c;
%     A_alt=(-b+sqrt(delta))./(2*a);
%     norm(A-A_alt)
end

function sigma2=Esigma(y,Phi,mu,Sigma,Phi2, c,d,n)
a=c+n/2;
PhiTPhiSigma=Phi2.*Sigma;
b=0.5*(norm(y-Phi*mu)+sum(PhiTPhiSigma(:))+2*d);

sigma2=b/a;

function [Sigma, mu]=weightEst(y,Phi,sigma2,A)
[m,~]=size(Phi);
Ainv=diag(1./A);
% Amatrix=inv(sigma2*eye(m)+Phi*Ainv*Phi'); %imagesc(log(abs(Amatrix)));
Sigma=Ainv-Ainv*Phi'*((sigma2*eye(m)+Phi*Ainv*Phi')\Phi*Ainv);
mu=Sigma*Phi'*y/sigma2;

function A=alphaEst(Sigma,mu,lambda)
A=3./(diag(Sigma)+mu.^2+lambda);

% function lambda=lambdaEst(theta,k,alpha)
% lambda=(theta+1)./(k+alpha/2);