function drawPlane(a,b,c,d,x_bounds,y_bounds,z_bounds,points)

% Name:         drawPlane
% Version:      1.0
% Date:         6 March 2017
% Author:       Ash Bellett
% Description:  Draws a bounded plane of the form ax+by+cz=0.
% Inputs:       a,b,c,d: plane coefficients
%               XYZ_bounds: [minimum_bound maximum_bound] of XYZ direction
% Outputs:      3D plot of plane

values = find([a,b,c,d]~=0);
values = mat2str(values);

switch values
    case '1'
        x = 0;
        y = linspace(y_bounds(1),y_bounds(2),points);
        z = linspace(z_bounds(1),z_bounds(2),points);
        fill3([0 0 0 0],[y_bounds(1) y_bounds(2) y_bounds(2) y_bounds(1)],[z_bounds(1) z_bounds(1) z_bounds(2) z_bounds(2)],[0.9336 0.9102 0.8281]);
%         for i = 1:length(y)
%         for j = 1:length(z)
%             plot3(x,y(i),z(j),'.r','MarkerSize',4)
%         end
%         end
    case '2'
        x = linspace(x_bounds(1),x_bounds(2),points);
        y = 0;
        z = linspace(z_bounds(1),z_bounds(2),points);
        fill3([x_bounds(1) x_bounds(2) x_bounds(2) x_bounds(1)],[0 0 0 0],[z_bounds(1) z_bounds(1) z_bounds(2) z_bounds(2)],[0.9336 0.9102 0.8281]);
%         for i = 1:length(x)
%         for j = 1:length(z)
%             plot3(x(i),y,z(j),'.r','MarkerSize',4)
%         end
%         end
    case '3'
        x = linspace(x_bounds(1),x_bounds(2),points);
        y = linspace(y_bounds(1),y_bounds(2),points);
        z = 0;
        fill3([x_bounds(1) x_bounds(2) x_bounds(2) x_bounds(1)],[y_bounds(1) y_bounds(1) y_bounds(2) y_bounds(2)],[0 0 0 0],[0.6471 0.5176 0.3725]);
%         for i = 1:length(x)
%         for j = 1:length(y)
%             plot3(x(i),y(j),z,'.r','MarkerSize',4)
%         end
%         end
    case '[1 2]'
        x = linspace(x_bounds(1),x_bounds(2),points);
        z = linspace(z_bounds(1),z_bounds(2),points);
        for i = 1:length(x)
        for j = 1:length(z)
            y(i,j) = -(a/b)*x(i);
            plot3(x(i),y(i,j),z(j),'.r','MarkerSize',4)
        end
        end
    case '[1 3]'
        x = linspace(x_bounds(1),x_bounds(2),points);
        y = linspace(y_bounds(1),y_bounds(2),points);
        for i = 1:length(x)
        for j = 1:length(y)
            z(i,j) = -(a/c)*x(i);
            plot3(x(i),y(j),z(i,j),'.r','MarkerSize',4)
        end
        end
    case '[1 4]'
        x = -d/a;
        y = linspace(y_bounds(1),y_bounds(2),points);
        z = linspace(z_bounds(1),z_bounds(2),points);
        fill3([x x x x],[y_bounds(1) y_bounds(2) y_bounds(2) y_bounds(1)],[z_bounds(1) z_bounds(1) z_bounds(2) z_bounds(2)],[0.9336 0.9102 0.8281]);
%         for i = 1:length(y)
%         for j = 1:length(z)
%             plot3(x,y(i),z(j),'.r','MarkerSize',4)
%         end
%         end
    case '[2 3]'
        x = linspace(x_bounds(1),x_bounds(2),points);
        y = linspace(y_bounds(1),y_bounds(2),points);
        for i = 1:length(x)
        for j = 1:length(y)
            z(i,j) = -(b/c)*y(j);
            plot3(x(i),y(j),z(i,j),'.r','MarkerSize',4)
        end
        end
     case '[2 4]'
        x = linspace(x_bounds(1),x_bounds(2),points);
        y = -d/b;
        z = linspace(z_bounds(1),z_bounds(2),points);
        fill3([x_bounds(1) x_bounds(2) x_bounds(2) x_bounds(1)],[y y y y],[z_bounds(1) z_bounds(1) z_bounds(2) z_bounds(2)],[0.9336 0.9102 0.8281]);
%         for i = 1:length(x)
%         for j = 1:length(z)
%             plot3(x(i),y,z(j),'.r','MarkerSize',4)
%         end
%         end
     case '[3 4]'
        x = linspace(x_bounds(1),x_bounds(2),points);
        y = linspace(y_bounds(1),y_bounds(2),points);
        z = -d/c;
        % fill3([x_bounds(1) x_bounds(2) x_bounds(2) x_bounds(1)],[y_bounds(1) y_bounds(1) y_bounds(2) y_bounds(2)],[z z z z],[0.5 0.5 0.5]);
%         for i = 1:length(x)
%         for j = 1:length(y)
%             plot3(x(i),y(j),z,'.r','MarkerSize',4)
%         end
%         end
    case '[1 2 3]'
        x = linspace(x_bounds(1),x_bounds(2),points);
        y = linspace(y_bounds(1),y_bounds(2),points);
        for i = 1:length(x)
        for j = 1:length(y)
            z(i,j) = -(a*x(i)+b*y(j))/c;
            plot3(x(i),y(j),z(i,j),'.r','MarkerSize',4)
        end
        end
    case '[1 2 4]'
        x = linspace(x_bounds(1),x_bounds(2),points);
        z = linspace(z_bounds(1),z_bounds(2),points);
        for i = 1:length(x)
        for j = 1:length(z)
            y(i,j) = -(a*x(i)+d)/b;
            plot3(x(i),y(i,j),z(j),'.r','MarkerSize',4)
        end
        end
    case '[1 3 4]'
        x = linspace(x_bounds(1),x_bounds(2),points);
        y = linspace(y_bounds(1),y_bounds(2),points);
        for i = 1:length(x)
        for j = 1:length(y)
            z(i,j) = -(a*x(i)+d)/c;
            plot3(x(i),y(j),z(i,j),'.r','MarkerSize',4)
        end
        end
    case '[2 3 4]'
        x = linspace(x_bounds(1),x_bounds(2),points);
        y = linspace(y_bounds(1),y_bounds(2),points);
        for i = 1:length(x)
        for j = 1:length(y)
            z(i,j) = -(b*y(j)+d)/c;
            plot3(x(i),y(j),z(i,j),'.r','MarkerSize',4)
        end
        end
    case '[1 2 3 4]'
        x = linspace(x_bounds(1),x_bounds(2),points);
        y = linspace(y_bounds(1),y_bounds(2),points);
        for i = 1:length(x)
        for j = 1:length(y)
            z(i,j) = -(a*x(i)+b*y(j)+d)/c;
            plot3(x(i),y(j),z(i,j),'.r','MarkerSize',4)
        end
        end
end
