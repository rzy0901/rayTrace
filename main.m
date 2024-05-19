close all; clear; clc;
%% static part
simulation_static
%% dynamic part
figure;
pic_num = 1;
for kk = 1:1:40
    % Plot environment
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
    drawPlane(A,B,C,D,[x_min,x_max],[y_min,y_max],[z_min,z_max])
    hold on;
    end
    view(5,45)
%     view(2)
    floor = imread('Floor Texture.jpg');
    image([0 6],[0 6],floor);
    lgd_tx = plot3(transmit_pos(1),transmit_pos(2),transmit_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',blue);
    hold on;
    lgd_rx = plot3(receive_pos(1),receive_pos(2),receive_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',green);

    % Plot human
    angstep = 0;
    [pos_ped, vel_ped, ax_ped] = move(ped,dt,angstep);
    % 16 joint points
    % [R-ankle L-ankle R-knee L-knee R-hip L-hip R-waist L-waist R-hand L-hand R-elbow L-elbow
    % R-shoulder L-shoulder Base Head]
    x1 = pos_ped(1,:).';
    y1 = pos_ped(2,:).';
    z1 = pos_ped(3,:).';
    % Plot bodies
    scatter3(x1([1:14,16]),y1([1:14,16]),z1([1:14,16]),30,[1 0 0],'filled');
    scatter3(x1(15),y1(15),z1(15),50,[1 0 0],'filled');  % Head
    lw = 2;
    plot3(x1(1:2:7),y1(1:2:7),z1(1:2:7),'Color','r','LineWidth',lw);
    plot3(x1(2:2:8),y1(2:2:8),z1(2:2:8),'Color','r','LineWidth',lw);
    plot3(x1(9:2:13),y1(9:2:13),z1(9:2:13),'Color','r','LineWidth',lw);
    plot3(x1(10:2:14),y1(10:2:14),z1(10:2:14),'Color','r','LineWidth',lw);
    plot3(x1([7,8]),y1([7,8]),z1([7,8]),'Color','r','LineWidth',lw);
    plot3(x1([13,14]),y1([13,14]),z1([13,14]),'Color','r','LineWidth',lw);
    plot3(x1([15,16]),y1([15,16]),z1([15,16]),'Color','r','LineWidth',lw);
    plot3([x1(15),(x1(7)+x1(8))/2],[y1(15),(y1(7)+y1(8))/2],[z1(15),(z1(7)+z1(8))/2],'Color',red,'LineWidth',lw);
    % Plot rays scattered by humans
    for ii = 1:1:16
        plot3([transmit_pos(1) x1(ii) receive_pos(1)],[transmit_pos(2) y1(ii) receive_pos(2)],[transmit_pos(3) z1(ii) receive_pos(3)],...
            'Color',red,'LineWidth',0.05);
    end
    for ii = 1:1:length(ray_connections) 
        flag1 = 0;
        connections = ray_connections{ii};
        for jj =1:1:size(connections,1)
            B = connections(jj,1:3);
            C = connections(jj,4:6);
            for m=1:size(x1,1)
                A = [x1(m) y1(m) z1(m)];
                distance = pointToLineDistance(A,B,C);
                if distance <= 0.4
                    flag1 = 1;
                end
            end
        end
        if flag1 == 0
            if size(connections,1) == 1
                LoS = plot3([connections(1) connections(4)],[connections(2) connections(5)],[connections(3) connections(6)],'Color',yellow,'LineWidth',2);
            else
                for jj = 1:1:size(connections,1)
                    NLoS = plot3([connections(jj,1) connections(jj,4)],[connections(jj,2) connections(jj,5)],[connections(jj,3) connections(jj,6)],'Color',blue,'LineWidth',0.5);
                end
            end
        end
    end
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
%     set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    set(gcf, 'Color', 'w');
    hold off
    axis equal
    try
        legend([lgd_tx, lgd_rx, LoS, NLoS],'Tx','Rx', 'LoS', 'NLoS');
    catch err
        legend([lgd_tx, lgd_rx, NLoS],'Tx','Rx', 'NLoS');
    end
%     axis off;
    drawnow;
    Frame=getframe(gcf);
    Image=frame2im(Frame);
    [Image,map]=rgb2ind(Image,256);
    if pic_num == 1
        imwrite(Image,map,'room1.gif','gif', 'Loopcount',inf,'DelayTime',2*dt);
    else
        imwrite(Image,map,'room1.gif','gif','WriteMode','append','DelayTime',2*dt);
    end
    pic_num = pic_num + 1;
end

function out = isintersected(P1,P2,Q1,Q2)
% 判断线段是否相交
% https://www.cnpython.com/qa/396746 回答二
a = dot(P2-P1,Q1-P1)/(norm(P2-P1)^2);
b = dot(P2-P1,Q2-Q1)/(norm(P2-P1)^2);
C = b*(P2-P1) - (Q2-Q1);
t1 = dot(C,Q1-(1-a)*P1-a*P2)/(norm(C)^2);
t0 = a + t1 * b;
if t1 >= 0 && t1 <= 1 && t0 >= 0 && t0 <= 1
    out = 1;
else
    out = 0;
end
end

% 给定一个三个三维坐标点a，b，c，计算a到直线b，c的垂直距离。
function distance = pointToLineDistance(a, b, c)
    v = c - b;
    w = a - b;
    proj_v_w = dot(w, v) / dot(v, v) * v;
    distance = norm(w - proj_v_w);
end
