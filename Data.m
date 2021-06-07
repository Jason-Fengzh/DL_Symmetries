function [phi,y,angle]= Data(j,K)
n=2*j+1;
X=randn(n);
Y=randn(n);
Z=X+1i*Y;
phi=Z/(norm(Z,'fro'));
% randomly generate the normalized matrix phi
y=zeros(n,n,K);
angle=[];
% the output data is in 3-dimensional array form
for i=1:1:K
    D=zeros(n,n);
    alpha=2*pi*rand(1);
    beta=2*pi*rand(1);
    gamma=2*pi*rand(1);
    anglei=[alpha;beta;gamma];
    angle=[angle,anglei];
    Vbeta=zeros(n,1);
    for m=1:1:n
        Vbeta(m)=exp(1i*beta*(m-j-1));
    end
    for m1=1:1:n
        for m2=1:1:n
            C=ComputeCo(j,m1-j-1,m2-j-1);
            D(m1,m2)=exp(-m1*1i*alpha)*(C.'*Vbeta)*exp(-m2*1i*gamma);
        end
    end
    % generate the Wigner D-matrix
    y(:,:,i)=D*phi;
end