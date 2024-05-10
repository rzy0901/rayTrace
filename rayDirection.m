function [direction_x,direction_y,direction_z] = rayDirection(origin,termination)

% Name:         rayDirection
% Version:      1.0
% Date:         25 April 2017
% Author:       Ash Bellett
% Description:  Determines direction between two 3D points
% Inputs:       origin: start point of ray
%               termination: end point of ray
% Outputs:      direction_x: x-direction from start to end point of ray
%               direction_y: y-direction from start to end point of ray
%               direction_z: z-direction from start to end point of ray

direction = termination-origin;
distance = sqrt(direction(1)^2+direction(2)^2+direction(3)^2);
direction = direction/distance;
direction_x = direction(1);
direction_y = direction(2);
direction_z = direction(3);

end