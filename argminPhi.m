function [phi] = argminPhi(y,z)
% the input y and z should be 3-dim array
% y(:,:,i) is a nxn matrix, 1<=i<=K
% z(:,:,i,j) is a nxn matrix, 1<=j<=q
N1=size(y);
N2=size(z);
n=N1(1);
K=N1(3);
if (length(N2)==4)
    q=N2(4);
end
if (length(N2)==3)
    q=1;
end
% cell structure can store matrix in each cell element, we can also use a 4-dim array
Z1=cell(q,q);
Z2=cell(q,1);
% Z1 and Z2 contain the "matrix element"
for j1=1:1:q
    for j2=1:1:q
        Z1{j1,j2}=zeros(n);
        for i=1:1:K
            Z1{j1,j2}=Z1{j1,j2}+z(:,:,i,j1)'*z(:,:,i,j2);
        end
    end
    Z2{j1}=zeros(n);
    for i=1:1:K
        Z2{j1}=Z2{j1}+z(:,:,i,j1)'*y(:,:,i);
    end
end
% first form the column
A1=cell(1,q);
for i2=1:1:q
    A1{i2}=[];
    for i1=1:1:q
        A1{i2}=[A1{i2};Z1{i1,i2}];
    end
end
% then form the two matrice we want
A=[];
B=[];
for i=1:1:q
    A=[A,A1{i}];
    B=[B;Z2{i}];
end 
C=A\B;
% finally divide the long matrix into phi_1,...,phi_q
phi=zeros(n,n,q);
for j=1:1:q
    phi(:,:,j)=C(1+(j-1)*n:j*n,:);
    phi(:,:,j)=phi(:,:,j)/norm(phi(:,:,j),'fro');  % make their norm equal to 1
end


