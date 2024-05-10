% Close all open figure windows:
close all
% Remove existing variables from memory:
clear
% Refresh command window:
clc
constants
parameters
%% Static part
environment(1) % figure(1)
source
modulator
transmitter
channel
relay
% channel
receiver
sink
pic_num = 1;
% set (gcf,'Position',[100,250,1000,600],'color','w');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'color','w');

axis equal;
% view(2)
  xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
view(5,45)
% view(30,30)
% set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'ztick',[],'zticklabel',[]);
% axis off;

%% Dynamic part, load the determinstic rays calculated before.
%% One person
figure(2);
set (gcf,'Position',[100,250,1000,600],'color','w');
% set (gcf,'color','w');
for kk = 1:1:40% kk = 1: 1: 40
    environment(2) % figure(2)
%     view(2)
    view(5,45);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    %% Walking
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
    plot3([x1(15),(x1(7)+x1(8))/2],[y1(15),(y1(7)+y1(8))/2],[z1(15),(z1(7)+z1(8))/2],'Color','r','LineWidth',lw);
    % Plot rays reflected by humans
    for ii = 1:1:16
        plot3([transmit_pos(1) x1(ii) receive_pos(1)],[transmit_pos(2) y1(ii) receive_pos(2)],[transmit_pos(3) z1(ii) receive_pos(3)],...
            'Color','r','LineWidth',0.05);
    end
    % Plot deterministic rays
%     for ii = 1 :1 : length(start_points)
% %         flag1 = 0; % 判断径是否相交
% %         for jj = 1:1:16
% %             % flag1 = isintersected(start_points(ii,:),end_points(ii,:),transmit_pos,pos_ped(:,jj).') || flag1;
% %             flag1 = isintersected(start_points(ii,:),end_points(ii,:),receive_pos,pos_ped(:,jj).') || flag1;
% %         end
% %         if flag1 == 0
% %             drawRay(start_points(ii,:),end_points(ii,:),ray_resolution);
% %         end
%         drawRay(start_points(ii,:),end_points(ii,:),ray_resolution);
%     end
    % Plot deterministic rays version 2 
    for ii = 1:1:length(ray_connections) % 遍历多径
        flag1 = 0;
        connections = ray_connections{ii};
        for jj =1:1:size(connections,1) % 检测相交
            B = connections(jj,1:3);
            C = connections(jj,4:6);
            for m=1:size(x1,1)
                A = [x1(m) y1(m) z1(m)];
                distance = pointToLineDistance(A,B,C);
                if distance <= 0.2
                    flag1 = 1;
                end
            end
        end
        if flag1 == 0
            if size(connections,1) == 1
                LoS = plot3([connections(1) connections(4)],[connections(2) connections(5)],[connections(3) connections(6)],'Color','y','LineWidth',2);
            else
                for jj = 1:1:size(connections,1)
                    NLoS = plot3([connections(jj,1) connections(jj,4)],[connections(jj,2) connections(jj,5)],[connections(jj,3) connections(jj,6)],'Color','b','LineWidth',0.005);
                end
            end
        end
    end
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
%     set (gcf,'Position',[100,250,1000,600],'color','w');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    hold off
    % legend('Transmitter','Receiver','Relay');
    hold off;
%     set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'ztick',[],'zticklabel',[]);
%     axis equal;
    axis off;
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