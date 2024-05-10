function drawRay(origin,termination,step_size)

% Name:         drawRay
% Version:      1.0
% Date:         19 February 2017
% Author:       Ash Bellett
% Description:  Draws a ray using start and end point and determines ray
%               distance
% Inputs:       origin: start point of ray
%               termination: end point of ray
%               step_size: number of points to draw line
% Outputs:      distance: distance from start and end point of ray

% direction = termination-origin;
% distance = sqrt(direction(1)^2+direction(2)^2+direction(3)^2);
% direction = direction/distance;
plot3([origin(1) termination(1)],[origin(2) termination(2)],[origin(3) termination(3)],'Color','b','LineWidth',0.005);
% for i = 0:step_size:distance
%     ray_point = origin+i*direction;
%     plot3(ray_point(1),ray_point(2),ray_point(3),'.k','MarkerSize',4)
% end
end