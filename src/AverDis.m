function [r,nearInx]=AverDis(x,y)
% Function AverDis
%   Calculate the average distance of the adjacent points specified by vector x and y
%   Method: Here a very inefficient method is used, that for each point, we calculate all distance with other points,
%       and find the minor one (to determine which one is the nearest). Better method must could be found, if one
%       has interest to improve it
%
% Input:
%   x, 1 x N verctor, x value of points
%   y, 1 x N verctor, y value of points
% Output:
%   r, double, average distance of all adjacent points
%   nearInx, 2 x N martrix, record the pair of nearest points
%       where nearInx(1,i) is index from 1 to N
%             nearInx(2,i) is the corresponding index for nearInx(1,i)

N=length(x);
N2=length(y);
if (N~=N2)
    error(sprintf('x and y should have the same number of points'))
end

dis=zeros(1,N);
sumdis=0;

for i=1:N
    nearInx(1,i)=i;
    for j=1:N
        if i~=j
            dis(j)=sqrt((x(j)-x(i))^2+(y(j)-y(i))^2); % Calculate distance
        else
            dis(j)=1000000000; % The point-self, set a big value
        end
    end
    [temp nearInx(2,i)]=min(dis); % The distance between ith point and the nearest point
    sumdis=sumdis+temp;
end

r=sumdis/(N-1); % The result
