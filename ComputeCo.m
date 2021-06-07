function [Coeffi] = ComputeCo(j,m,m1)
coeffi=zeros(4*j+1,1);
M2=min([j-m1,j+m]);
M1=max([0,m-m1]);
% [M1,M2] is the range of index k
for k=M1:1:M2
    n=2*k-m+m1;
    ck=(-1)^(k-m+m1)*sqrt(factorial(j+m)*factorial(j-m)*factorial(j+m1)*factorial(j-m1))/(factorial(j+m-k)*factorial(j-m1-k)*factorial(k-m+m1)*factorial(k));
    A=zeros(4*j+1,1); % A is the coefficient vector of a fixed k
    for t=0:1:n
        B=zeros(4*j+1,1);
        for r=0:1:2*j-n
            B(1+2*r+2*t)=nchoosek(2*j-n,r)*nchoosek(n,t);
        end
        A=A+B*(-1)^t;   % As mentioned in PDF
    end    
    an=(0.5)^(2*j)*(1i)^n;    
    coeffi=coeffi+A*an*ck;
end
Coeffi=zeros(2*j+1,1);
for s=1:1:2*j+1
    Coeffi(s)=coeffi(2*s-1);
end
    
% We obtain a (4j+1)*1 "coeffi" vector, and the coefficient of exp(ik\beta) equals to coeffi(2j+1+k)

% [Coeffi(4,3,2),ComputeCo(4,3,2)] The results by this method are same as the integral expression