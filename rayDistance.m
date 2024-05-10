function [distance] = rayDistance(origin,termination)

% Name:         rayDistance
% Version:      1.0
% Date:         25 April 2017
% Author:       Ash Bellett
% Description:  Determines distance between two 3D points
% Inputs:       origin: start point of ray
%               termination: end point of ray
% Outputs:      distance: distance from start and end point of ray

direction = termination-origin;
distance = sqrt(direction(1)^2+direction(2)^2+direction(3)^2);

end