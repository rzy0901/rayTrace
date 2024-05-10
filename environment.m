function environment(figure_counter)

constants
parameters
% % Set environment geometry and material properties.
% 
% % Transmitter and receiver positions:
% transmit_pos = [0.5,1.5,1];
% receive_pos = [5.5,5.5,1];
% relay_pos = [0.5771,5.4036,1];
% 
% % Geometric properties of objects:
% % Normal, point and bounds.
% 
% object_geometry(1,1:12) = [[0,0,1],[3,3,0],[0,6],[0,6],[-99,99]];
% %object_geometry(2,1:12) = [[0,1,0],[3,0,1],[0,6],[-99,99],[0,2]];
% object_geometry(3,1:12) = [[0,1,0],[3,6,1],[0,6],[-99,99],[0,2]];
% object_geometry(4,1:12) = [[1,0,0],[0,3,1],[-99,99],[0,6],[0,2]];
% %object_geometry(5,1:12) = [[1,0,0],[6,3,1],[-99,99],[0,6],[0,2]];
% %object_geometry(6,1:12) = [[1,0,0],[3,3,1],[-99,99],[2.5,3.5],[0.8,1.2]];
% 
% % object_geometry(1,1:12) = [[1,0,0],[4,0,0],[-99,99],[0,4],[0,4]];
% % object_geometry(2,1:12) = [[0,1,0],[0,4,0],[0,4],[-99,99],[0,4]];
% % object_geometry(3,1:12) = [[0,0,1],[0,0,0],[0,4],[0,4],[-99,99]];
% 
% % Matrial properties of objects;
% % Reflectance.
% object_material(1,1) = 0.4;
% object_material(2,1) = 0.4;
% object_material(3,1) = 0.4;
% object_material(4,1) = 0.4;
% object_material(5,1) = 0.4;
% object_material(6,1) = 0.4;
% 
% [object_number,~] = size(object_geometry);

%figure_counter = 1;

% Transmitter, relay and receiver:
figure(figure_counter)
% plot3(transmit_pos(1),transmit_pos(2),transmit_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',blue)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on;
% plot3(receive_pos(1),receive_pos(2),receive_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',green)
% % plot3(relay_pos(1),relay_pos(2),relay_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',yellow)

% Objects:
[obstacle_number,~]=size(object_geometry);
% Ax + By + Cz + D = 0;
for i = 1:obstacle_number
    A = object_geometry(i,1);
    B = object_geometry(i,2);
    C = object_geometry(i,3);
    D = -(A*object_geometry(i,4)+B*object_geometry(i,5)+C*object_geometry(i,6));
    x_min = object_geometry(i,7);
    x_max = object_geometry(i,8);
    y_min = object_geometry(i,9);
    y_max = object_geometry(i,10);
    z_min = object_geometry(i,11);
    z_max = object_geometry(i,12);
    drawPlane(A,B,C,D,[x_min,x_max],[y_min,y_max],[z_min,z_max],object_resolution)
    hold on;
end
floor = imread('Floor Texture.jpg');
image([0 6],[0 6],floor);
% % Transmitter, relay and receiver:
% figure(figure_counter+1)
% hold on
% view(30,30)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% plot3(transmit_pos(1),transmit_pos(2),transmit_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',blue)
% plot3(receive_pos(1),receive_pos(2),receive_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',green)
% plot3(relay_pos(1),relay_pos(2),relay_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',yellow)


plot3(transmit_pos(1),transmit_pos(2),transmit_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',blue)
xlabel('x')
ylabel('y')
zlabel('z')
hold on;
plot3(receive_pos(1),receive_pos(2),receive_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',green)
% plot3(relay_pos(1),relay_pos(2),relay_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',yellow)

% % Objects:
% [obstacle_number,~]=size(object_geometry);
% for i = 1:obstacle_number
%     A = object_geometry(i,1);
%     B = object_geometry(i,2);
%     C = object_geometry(i,3);
%     D = -(A*object_geometry(i,4)+B*object_geometry(i,5)+C*object_geometry(i,6));
%     x_min = object_geometry(i,7);
%     x_max = object_geometry(i,8);
%     y_min = object_geometry(i,9);
%     y_max = object_geometry(i,10);
%     z_min = object_geometry(i,11);
%     z_max = object_geometry(i,12);
%     drawPlane(A,B,C,D,[x_min,x_max],[y_min,y_max],[z_min,z_max],object_resolution)
% end