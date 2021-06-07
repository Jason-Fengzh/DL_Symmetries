function [unitphi] = normalize(phi)
% phi is a nxm matrix. we want each column of phi to be unit vector
s=size(phi);
n=s(1);
m=s(2);
unitphi=zeros(n,m);
for i=1:1:m
    x=phi(:,i);
    x=x/norm(x);
    unitphi(:,i)=x;
end
