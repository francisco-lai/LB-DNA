function [Gr, Ga] = make_greenf(ee, ii, Ln, HH)
%MAKE_GREENF Summary of this function goes here
%   Detailed explanation goes here
    Gr=(ee(ii)*eye(2*Ln)-HH)\eye(2*Ln); Ga=Gr'; %Green's function  matrices
end

