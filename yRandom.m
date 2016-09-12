function [ y ] = yRandom( number,mean,deviation )
%YRANDOM randomly generate the values every time call the function
%   With mean and standard deviation
rng('shuffle');
y=deviation.*randn(number,1)+mean;
end

