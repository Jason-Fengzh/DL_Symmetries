function [x] = argminX(y,Phi,c,lmd)
% y should be (2j+1)^2 x 1 vector, Phi is (2j+1)^2 x c matrix
n=sqrt(length(y)); % n=2j+1
j=(n-1)/2;

cvx_begin quiet

variable x(c,1) complex

normx=0;
for i=1:1:c
    normx=normx+abs(x(i));
end

minimize((y-Phi*x)'*(y-Phi*x)/2+lmd*norm(x,1))

    subject to
    
   

cvx_end
