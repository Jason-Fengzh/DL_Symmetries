function [d,index] = idist(phi1,phi2,N,interval)
% input 'interval' should be a 3x2 matrix where the rows describe the range of alpha/ beta/ gamma
% output 'd' is the smallest distance up to some D, 'index' is a 1x3 matrix showing the value of alpha/ beta/ gamma
s=size(phi1);
n=s(1);
j=(n-1)/2;
for i=1:1:3
    step(i)=(interval(i,2)-interval(i,1))/N;
end
minimum=10;
index=ones(1,3);
for i1=1:1:N+1
    for i2=1:1:N+1
        for i3=1:1:N+1
            alpha = interval(1,1)+step(1)*(i1-1);
            beta = interval(2,1)+step(2)*(i2-1);
            gamma = interval(3,1)+step(3)*(i3-1);
            D=zeros(n,n);
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
            theta=-angle(trace(phi1'*(D*phi2)));
            d1=norm(phi1-exp(1j*theta)*D*phi2,'fro');
            if (d1<minimum)
                minimum=d1;
                index=[i1-1,i2-1,i3-1];
            end
        end
    end
end
d=minimum^2;
