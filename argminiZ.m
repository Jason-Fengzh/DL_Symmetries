function [z] = argminiZ(y,phi,lmd)
% y and phi are 2j+1 x 2j+1 complex matrice
N=size(y);
n=N(1);    % n=2j+1
j=(n-1)/2;

cvx_begin quiet

% let q=1
variable Z(n,n,n,n,n,n) complex
variable M(n^3+1,n^3+1) hermitian
variable Z_dim1(n^6,1) complex
variable z(n,n) complex
variable x(n,n,n) complex
variable t

% Fnorm=((norm(y-z*phi,'fro'))^2)/2;
B=y-z*phi;
B1=reshape(B,1,n*n);
J=sub2ind([n,n,n,n,n,n],j+1,j+1,j+1,j+1,j+1,j+1);
[j1,j2]=ind2sub([n^3,n^3],J);

minimize((M(n^3+1,n^3+1)/2+M(j1,j2)/2)*lmd+B1*B1'/2)

    subject to
    
    reshape(z,n^2,1)==LinearL(j)*reshape(x,n^3,1);
    count=0;
    for i2=1:1:n
        for i3=1:1:n
            for i5=1:1:n
                for i6=1:1:n
                    for i1=1:1:n-1
                        for i4=1:1:n-1
                            count=count+1;
                            indexL(count)=sub2ind([n,n,n,n,n,n],i1,i2,i3,i4,i5,i6);
                            indexR(count)=sub2ind([n,n,n,n,n,n],i1+1,i2,i3,i4+1,i5,i6);
                        end
                    end

                end
            end
        end
    end

    for i1=1:1:n
        for i3=1:1:n
            for i4=1:1:n
                for i6=1:1:n
                    for i2=1:1:n-1
                        for i5=1:1:n-1
                            count=count+1;
                            indexL(count)=sub2ind([n,n,n,n,n,n],i1,i2,i3,i4,i5,i6);
                            indexR(count)=sub2ind([n,n,n,n,n,n],i1,i2+1,i3,i4,i5+1,i6);
                        end
                    end                                  
                end
            end
        end
    end

    for i1=1:1:n
        for i2=1:1:n
            for i4=1:1:n
                for i5=1:1:n
                    for i3=1:1:n-1
                        for i6=1:1:n-1
                            count=count+1;
                            indexL(count)=sub2ind([n,n,n,n,n,n],i1,i2,i3,i4,i5,i6);
                            indexR(count)=sub2ind([n,n,n,n,n,n],i1,i2,i3+1,i4,i5,i6+1);
                        end
                    end                                  
                end
            end
        end
    end

    A=sparse(count,n^6);
    row=0;
    for i=1:1:count
        row=row+1;
        A(row,indexL(row))=1;
        A(row,indexR(row))=-1;
    end
    Z0=zeros(count,1);
    Z1=reshape(Z,n^3,n^3);
    x1=reshape(x,n^3,1);

    M == [Z1,x1;x1',t];
    M == hermitian_semidefinite(n^3+1); % PSD
    
    % Z is block toeplitz
    Z_dim1 == reshape(Z,[n^6,1]);
    A*Z_dim1 == Z0;

cvx_end
% z
% normz=norm(z,'fro')