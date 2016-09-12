function [ dist ] = distance( X,Y )
%DISTANCE compute the Euler distance  between X,Y vectors
%   Detailed explanation goes here
dist=sum((X-Y).^2);
end

