function [distance,angle] = infidist(phi1,phi2,N,interval,times)
% input 'times' is how many times we repeat the binary search type of procedure
% 'distance' is the inf distance and 'angle' is the value of [alpha,beta,gamma]
d=zeros(1,times);
angle=zeros(3,1);
for k=1:1:times
    if(k==1)
        [d(k),index]=idist(phi1,phi2,N,interval);
    end
    if(k>1)
        interv=interval;
        for i=1:1:3
            interval(i,1)=interv(i,1)+(index(i)-1)*(interv(i,2)-interv(i,1))/N;
            interval(i,2)=interv(i,1)+(index(i)+1)*(interv(i,2)-interv(i,1))/N;
        end
        [d(k),index]=idist(phi1,phi2,N,interval);
    end
end
for j=1:1:3
    angle(j)=interval(j,1)+index(i)*(interval(j,2)-interval(j,1))/N;
end
distance=d(times);


