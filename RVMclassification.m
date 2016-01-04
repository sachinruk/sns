function [mu, Sigma,C,Ey,Ey2]=RVMclassification(t,X,iter)
mu=X\t;
A=ones(size(mu));
% XtX=X'*X; step=length(XtX);

for i=1:iter
    [C, Ey, Ey2]=q_y(X,t,mu);
%     mu=q_beta(XtX,step,Ey,A);
    [Sigma, mu]=weightEst(Ey,X,A);
    A=alphaEst(Sigma,mu);
end

% function mu=q_beta(XtX,step,Ey,A)
% SigmaInv=XtX;
% SigmaInv(1:(step+1):end)=SigmaInv(1:(step+1):end)+A;
% SigmaX=SigmaInv\X';
% mu=SigmaX*Ey;

function [C, Ey, Ey2]=q_y(X,t,mu)

mu2=X*mu;
C=sqrt(2*pi)*normcdf(-mu2);
C(t==1)=sqrt(2*pi)-C(t==1);

Ey=(2*t-1).*exp(-0.5*mu2.^2)./C+mu2;
Ey2=(2*t-1).*mu2.*exp(-0.5*mu2.^2)./C+mu2.^2+1;

function [Sigma, mu]=weightEst(y,Phi,A)
[m,~]=size(Phi);
Ainv=diag(1./A);
% Amatrix=inv(sigma2*eye(m)+Phi*Ainv*Phi'); %imagesc(log(abs(Amatrix)));
Sigma=Ainv-Ainv*Phi'*((eye(m)+Phi*Ainv*Phi')\Phi*Ainv);
mu=Sigma*Phi'*y;

function A=alphaEst(Sigma,mu)
% A=3./(diag(Sigma)+mu.^2+lambda);
delta=1e-9;
A=(0.5+delta)./(diag(Sigma)+mu.^2+delta);
