close all; clear; clc;
%% Constant
% Speed of light:
light_speed = 299792458;
% Colours:
red = [1,10/51,10/51];
blue = [14/51,32/51,1];
green = [0,40/51,10/51];
yellow = [50/51,40/51,4/17];
%% parameters
% Maximum number of ray collisions:
ray_collisions = 3;
% Modulation scheme:
% 1 = BPSK, 2 = 4-QPSK, 3 = 16-QAM.
modulation = 2;
% Message length:
message_length = 32;           % 32 bits
% Bandwidth:
bandwidth = 1*10^6;             % 1 MHz
% Sample factor:
sample_factor = 10000;          % 10 kS
% Carrier frequency:
frequency_carrier = 2.4*10^9;   % 2.4 GHz
% Transmit signal amplitude:
transmit_amplitude = 1;         % 1 V
% Transmitter gain:
transmit_gain = 10;             % 10 linear == 10 dB
% Receiver gain:
receive_gain = 10;              % 10 linear == 10 dB
% Noise variance:
noise_var = 0.01;
% Set environment geometry and material properties.
% Transmitter and receiver positions:
transmit_pos = [3,5.8,1];
receive_pos = [3,0.2,1];
% Geometric properties of objects: Ax + By + Cz + D = 0;
% Normal, point and bounds.
% [[A B C],[x y z],[x_min x_max],[y_min,y_max],[z_min,z_max]]
object_geometry(1,1:12) = [[0,0,1],[3,3,0],[0,6],[0,6],[0,2.5]];
object_geometry(2,1:12) = [[0,1,0],[3,6,1],[0,6],[0,6],[0,2.5]];
object_geometry(3,1:12) = [[1,0,0],[0,3,1],[0,6],[0,6],[0,2.5]];
object_geometry(4,1:12) = [[1,0,0],[6,3,1],[0,6],[0,6],[0,2.5]];
% Matrial properties of objects;
% Reflectance.
object_material(1,1) = 0.4;
object_material(2,1) = 0.4;
object_material(3,1) = 0.4;
object_material(4,1) = 0.4;
[object_number,~] = size(object_geometry);
% Generate a human.
ped = backscatterPedestrian('Height',1.7,'WalkingSpeed',1,'InitialPosition',[1;3;0],'OperatingFrequency',frequency_carrier);
dt = 0.1; % second
%% Environment
figure(1);
 
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
    drawPlane(A,B,C,D,[x_min,x_max],[y_min,y_max],[z_min,z_max]);
    hold on;
end
floor = imread('Floor Texture.jpg');
image([0 6],[0 6],floor);
view(5,45)
lgd_tx = plot3(transmit_pos(1),transmit_pos(2),transmit_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',blue);
xlabel('X (m)')
ylabel('Y (m)')
zlabel(' (m)')
lgd_rx = plot3(receive_pos(1),receive_pos(2),receive_pos(3),'x','MarkerSize',16,'LineWidth',3,'Color',green);
%% Source
% Create binary message.
for i = 1:message_length
    message(i) = randi([0,1],1);
end
%% modulator
% Convert message to passband signal using specified modulation scheme.
switch modulation
    case 1
        modulated_length = message_length;
        message_modulated = zeros(1,modulated_length);
        for i = 1:modulated_length
            switch message(i)
                case 1
                    message_modulated(i) = 1+0j;
                case 0
                    message_modulated(i) = -1+0j;
                otherwise
                    message_modulated(i) = 0+0j;
            end
        end
    case 2
        modulated_length = message_length/2;
        message_modulated = zeros(1,modulated_length);
        for i = 1:modulated_length
            switch int2str(message(2*i-1:2*i))
                case '0  0'
                    message_modulated(i) = 1/sqrt(2)+1j/sqrt(2);
                case '0  1'
                    message_modulated(i) = -1/sqrt(2)+1j/sqrt(2);
                case '1  1'
                    message_modulated(i) = -1/sqrt(2)-1j/sqrt(2);
                case '1  0'
                    message_modulated(i) = 1/sqrt(2)-1j/sqrt(2);
                otherwise
                    message_modulated(i) = 0+0j;
            end
        end
    case 3
        modulated_length = message_length/4;
        message_modulated = zeros(1,modulated_length);
        for i = 1:modulated_length
            switch int2str(message(4*i-3:4*i))
                case '0  0  0  0'
                    message_modulated(i) = 1/sqrt(11/2);
                case '0  0  0  1'
                    message_modulated(i) = 1j/sqrt(11/2);
                case '0  0  1  0'
                    message_modulated(i) = -1/sqrt(11/2);
                case '0  0  1  1'
                    message_modulated(i) = -1j/sqrt(11/2);
                case '0  1  0  0'
                    message_modulated(i) = 3/(2*sqrt(11/2))+1j*3/(2*sqrt(11/2));
                case '0  1  0  1'
                    message_modulated(i) = -3/(2*sqrt(11/2))+1j*3/(2*sqrt(11/2));
                case '0  1  1  0'
                    message_modulated(i) = -3/(2*sqrt(11/2))-1j*3/(2*sqrt(11/2));
                case '0  1  1  1'
                    message_modulated(i) = 3/(2*sqrt(11/2))-1j*3/(2*sqrt(11/2));
                case '1  0  0  0'
                    message_modulated(i) = 2/sqrt(11/2);
                case '1  0  0  1'
                    message_modulated(i) = (1j*2)/sqrt(11/2);
                case '1  0  1  0'
                    message_modulated(i) = -2/sqrt(11/2);
                case '1  0  1  1'
                    message_modulated(i) = (-1j*2)/sqrt(11/2);
                case '1  1  0  0'
                    message_modulated(i) = 5/(2*sqrt(11/2))+1j*5/(2*sqrt(11/2));
                case '1  1  0  1'
                    message_modulated(i) = -5/(2*sqrt(11/2))+1j*5/(2*sqrt(11/2));
                case '1  1  1  0'
                    message_modulated(i) = -5/(2*sqrt(11/2))-1j*5/(2*sqrt(11/2));
                case '1  1  1  1'
                    message_modulated(i) = 5/(2*sqrt(11/2))-1j*5/(2*sqrt(11/2));
                otherwise
                    message_modulated(i) = 0+0j;
            end
        end
    otherwise
        fprintf('Specify correct modulation scheme.')
end
%% transmitter
frequency_sample = sample_factor*bandwidth;
time = 0:1/frequency_sample:(modulated_length+1)/bandwidth;

message_real = real(message_modulated);
message_imag = imag(message_modulated);

tx_transmit_real = zeros(1,length(time));
tx_transmit_imag = zeros(1,length(time));

for i = 1:modulated_length
    tx_transmit_real = tx_transmit_real + message_real(i)*sinc(bandwidth*time-i);
    tx_transmit_imag = tx_transmit_imag + message_imag(i)*sinc(bandwidth*time-i);
end

signal_transmit = transmit_amplitude*transmit_gain*(tx_transmit_real+1j.*tx_transmit_imag);
tx_signal = signal_transmit;
figure(2);
subplot(2,1,1)
hold on
title('Real transmitted baseband signal')
plot(1000*time,tx_transmit_real,'Color',blue)
plot(1000*(1:modulated_length)/bandwidth,message_real,'.','Color',blue,'MarkerSize',10)
xlabel('Time (ms)')
ylabel('Amplitude')
subplot(2,1,2)
hold on
title('Imaginary transmitted baseband signal')
plot(1000*time,tx_transmit_imag,'Color',red)
plot(1000*(1:modulated_length)/bandwidth,message_imag,'.','Color',red,'MarkerSize',10)
xlabel('Time (ms)')
ylabel('Amplitude')
origin = transmit_pos;
termination = receive_pos;
%% Channel
% Number of possible path combinations (excluding LOS):
path_combination_number = 2^object_number-1;
% Binary vector representing path combinations:
path_combinations = de2bi(0:path_combination_number,object_number);

% Determine maximum number of rays to cast:
ray_number = 0;
for i = 1:ray_collisions
    if i <= object_number
        ray_number = ray_number + factorial(object_number)/factorial(object_number-i);
    end
end
% Account for line-of-sight ray:
ray_number = ray_number + 1;

% Initialise path information matrix:
ray_connections= cell(ray_number,1);
path_matrix = zeros(ray_number,6);
path_matrix(:,3:6) = 1;

% Initialise valid path counter:
m = 1;

% Switch back to ray tracing plot:
figure(1);

% For each path combination:
for i = 1:path_combination_number+1
    % Get a path combination:
    path_combination_current = find(path_combinations(i,:)~=0);
    % As long as path is below maximum collisions numbers:
    if length(path_combination_current) <= ray_collisions
        ray_connection = [];
        % All path permutations of particular path combination:
        path_permutations = unique(perms(path_combination_current),'rows');
        % Get number of permutations of particular path combination:
        [path_permutation_number,~] = size(path_permutations);
        % For each path permutation of particular path combination:
        for j = 1:path_permutation_number
            % Get a path permutation:
            path = path_permutations(j,:);
            % Number of collisions along path:
            path_collisions = length(path);
            % Initialise mirrored points matrix:
            mirror_matrix = zeros(path_collisions+1,3);
            % Set first point to receiver position:
            mirror_matrix(1,:) = termination;
            %%% METHOD OF IMAGES %%%
            % For each object in path:
            for k = 1:path_collisions
                % Get point to mirror:
                x = mirror_matrix(k,1);
                y = mirror_matrix(k,2);
                z = mirror_matrix(k,3);
                % Get object plane equation coefficients:
                a = object_geometry(path(end+1-k),1);
                b = object_geometry(path(end+1-k),2);
                c = object_geometry(path(end+1-k),3);
                d = -(a*object_geometry(path(end+1-k),4)+b*object_geometry(path(end+1-k),5)+c*object_geometry(path(end+1-k),6));
                % Perpendicular distance from point to plane:
                plane_distance = -(a*x+b*y+c*z+d)/(a*a+b*b+c*c);
                % Closest point on plane:
                plane_point = [x,y,z]+plane_distance*[a,b,c];
                % Vector towards plane:
                point_difference = plane_point-[x,y,z];
                % Mirrored point:
                mirror_matrix(k+1,1) = x+2*point_difference(1);
                mirror_matrix(k+1,2) = y+2*point_difference(2);
                mirror_matrix(k+1,3) = z+2*point_difference(3);
            end
            %%% PATH POINTS %%%
            % Initialise path points matrix:
            points_matrix = zeros(path_collisions+2,3);
            % Set origin to transmitter position:
            points_matrix(1,:) = origin;
            % Set termination to receiver position:
            points_matrix(end,:) = termination;
            % Initialise validity flag:
            flag_ray = 0;
            % For each intersection point along the path:
            for k = 1:path_collisions
                % Get ray segment origin:
                A = points_matrix(k,:);
                % Get ray segment termination:
                B = mirror_matrix(end+1-k,:);
                % Distance between origin and mirrored point:
                distance_full = rayDistance(A,B);
                % Direction towards mirrored point:
                [direction_x,direction_y,direction_z] = rayDirection(A,B);
                % Determine distance to object:
                distance = (dot(object_geometry(path(k),1:3),object_geometry(path(k),4:6))-dot(object_geometry(path(k),1:3),A))/dot(object_geometry(path(k),1:3),[direction_x,direction_y,direction_z]);
                % If ray is parallel to plane:
                if distance == 0
                    % Set invalidity flag:
                    flag_ray = 1;
                end
                % Intersection point with object (not considering bounds):
                points_matrix(k+1,:) = A+distance*[direction_x,direction_y,direction_z];
            end
            % If ray path may be valid:
            if flag_ray == 0
                % Initialise occlude flag:
                flag_occlude = 0;
                % For each ray segment:
                for k = 1:path_collisions+1
                    % Get origin of ray segment:
                    A = points_matrix(k,:);
                    % Get termination of ray segment:
                    B = points_matrix(k+1,:);
                    % Determine direction of ray segment:
                    [direction_x,direction_y,direction_z] = rayDirection(A,B);
                    % Determine distance of ray segment:
                    distance = rayDistance(A,B);
                    %%% OCCLUSION CHECK %%%
                    % For each object:
                    for l = 1:object_number
                        % If ray and plane are not parallel and intersection point is not with target object itself or intersection is at receiver:
                        if k == path_collisions+1 || (dot([direction_x,direction_y,direction_z],object_geometry(l,1:3))~=0 && l~=path(k))
                            % Distance to intersection:
                            distance_occlude = (dot(object_geometry(l,1:3),object_geometry(l,4:6))-dot(object_geometry(l,1:3),A))/dot(object_geometry(l,1:3),[direction_x,direction_y,direction_z]);
                            % Intersection point:
                            intersection_occlude = A+distance_occlude*[direction_x,direction_y,direction_z];
                            % Check if intersection lies in obstacle bounds:
                            if intersection_occlude(1) <= object_geometry(l,8)
                            if intersection_occlude(1) >= object_geometry(l,7)
                            if intersection_occlude(2) <= object_geometry(l,10)
                            if intersection_occlude(2) >= object_geometry(l,9)
                            if intersection_occlude(3) <= object_geometry(l,12)
                            if intersection_occlude(3) >= object_geometry(l,11)
                                % Check if intersection occurs along positive direction of ray:
                                if dot(B-A,intersection_occlude-A) > 0 && dot(B-A,intersection_occlude-A) < distance^2
                                    % Assert occlude flag:
                                    flag_occlude = 1;
                                end
                            end
                            end
                            end
                            end
                            end
                            end
                        end
                    end
                    %%% OBJECT BOUNDS CHECK %%%
                    % Line-of-sight path:
                    if isempty(path)
                        LoS = plot3([A(1) B(1)],[A(2) B(2)],[A(3) B(3)],'Color',yellow,'LineWidth',2);
                        ray_connection = [ray_connection;[A B]];
                        path_matrix(m,1) = distance;
                    % Multi-path segments:
                    elseif flag_occlude == 0 && k ~= path_collisions+1
                        % Check whether point lies in object bounds:
                        if B(1) <= object_geometry(path(k),8)
                        if B(1) >= object_geometry(path(k),7)
                        if B(2) <= object_geometry(path(k),10)
                        if B(2) >= object_geometry(path(k),9)
                        if B(3) <= object_geometry(path(k),12)
                        if B(3) >= object_geometry(path(k),11)
                            % Draw ray segment:
                            NLoS = plot3([A(1) B(1)],[A(2) B(2)],[A(3) B(3)],'Color',blue,'LineWidth',0.5);
                            ray_connection = [ray_connection;[A B]];
                            % Update path distance:
                            path_matrix(m,1) = path_matrix(m,1)+distance;
                            % Update attenuation due to collisions:
                            path_matrix(m,4) = path_matrix(m,4)*object_material(path(k));
                        end
                        end
                        end
                        end
                        end
                        end
                    % Terminating at receiver segment:
                    elseif flag_occlude == 0
                        % Check whether point lies in object bounds:
                        if A(1) <= object_geometry(path(end),8)
                        if A(1) >= object_geometry(path(end),7)
                        if A(2) <= object_geometry(path(end),10)
                        if A(2) >= object_geometry(path(end),9)
                        if A(3) <= object_geometry(path(end),12)
                        if A(3) >= object_geometry(path(end),11)
                            % Draw ray segment:
                            NLoS = plot3([A(1) B(1)],[A(2) B(2)],[A(3) B(3)],'Color',blue,'LineWidth',0.005);
                            ray_connection = [ray_connection;[A B]];
                            % Update path distance:
                            path_matrix(m,1) = path_matrix(m,1)+distance;
                            if mod(length(path),2) == 1
                                path_matrix(m,6) = -1;
                            end
                        end
                        end
                        end
                        end
                        end
                        end
                    end
                end
            end
        end
        % Update valid path counter
        ray_connections{m} = ray_connection;
        m = m + 1;
    end
end

legend([lgd_tx, lgd_rx, LoS, NLoS],'Tx','Rx', 'LoS', 'NLoS');
% Clean path information matrix:
index = any(path_matrix(:,1),2);
path_matrix(~index,:) = [];

new_ray_connections = {};
for i = 1:numel(ray_connections)
    if index(i) == 1
        new_ray_connections{end+1,1} = ray_connections{i};
    end
end

ray_connections = new_ray_connections;
clear new_ray_connections;

% Update path information matrix:
path_matrix(:,2) = path_matrix(:,1)/light_speed; % delay
path_matrix(:,3) = (light_speed/(frequency_carrier*sqrt(4*pi)))./path_matrix(:,1); % attenuation for free space model. eq (10) of https://ieeexplore.ieee.org/document/10525191
path_matrix(:,5) = exp((-1j*2*pi*frequency_carrier/light_speed).*path_matrix(:,1)); % phase change
% Get number of valid paths:
[paths_valid,~] = size(path_matrix);

% Initialise ray matrix:
ray_matrix = zeros(paths_valid,3);
% Populate ray matrix:
ray_matrix(:,1) = path_matrix(:,2);
ray_matrix(:,2) = path_matrix(:,3).*path_matrix(:,4);
ray_matrix(:,3) = path_matrix(:,5).*path_matrix(:,6);

% Plot power delay profile:
figure(3)
hold on
plot(1e9*ray_matrix(:,1),10*log10(ray_matrix(:,2)),'x','MarkerSize',10)
title('Power delay profile')
xlabel('Time (ns)')
ylabel('Received signal power (dB)')

% Apply attenuations and delay to transmitted signal:
signal_receive = zeros(1,length(time));
fade_coefficients = zeros(1,paths_valid);
for i = 1:paths_valid
    delay = round(frequency_sample*ray_matrix(i,1));
    fade_coefficients(i) = ray_matrix(i,2).*ray_matrix(i,3);
    signal_receive = signal_receive+fade_coefficients(i)*[zeros(1,delay),signal_transmit(1:end-delay)];
end

% Add AWGN to received signal:
for i = 1:length(time)
    signal_receive(i) = signal_receive(i) + noise_var*randn(1)+1j*noise_var*randn(1);
end
%% receiver
re_rx_attenuation = abs(sum(fade_coefficients));
re_rx_phase_shift = angle(sum(fade_coefficients));

rx_signal_receive = signal_receive*exp(-1j*re_rx_phase_shift);

rx_receive_real = receive_gain*real(rx_signal_receive);
rx_receive_imag = receive_gain*imag(rx_signal_receive);
rx_signal = rx_receive_real+1j*rx_receive_imag;

figure(4)
subplot(2,1,1)
plot(1000*time,rx_receive_real,'Color',blue)
title('Real received baseband signal')
xlabel('Time (ms)')
ylabel('Amplitude')
subplot(2,1,2)
plot(1000*time,rx_receive_imag,'Color',red)
title('Imaginary received baseband signal')
xlabel('Time (ms)')
ylabel('Amplitude')

message_receive_real = zeros(1,modulated_length);
message_receive_imag = zeros(1,modulated_length);

for i = 1:modulated_length
    message_receive_real(i) = rx_receive_real(i*frequency_sample/bandwidth);
    message_receive_imag(i) = rx_receive_imag(i*frequency_sample/bandwidth);
end

rx_message = message_receive_real+1j*message_receive_imag;

figure(5)
hold on
plot(real(message_modulated),imag(message_modulated),'x','Color',blue,'MarkerSize',10,'LineWidth',2)
plot(real(rx_message),imag(rx_message),'x','Color',red,'MarkerSize',10)
title('Signal Constellations')
legend('Transmitted constellation','Received constellation','Location','SouthEast')
axis([min(real(message_modulated))-2,max(real(message_modulated))+2,min(imag(message_modulated))-2,max(imag(message_modulated))+2])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
%% sink
power_transmit = 0;
power_receive = 0;
for i = 1:length(time)
    power_transmit = power_transmit + abs(tx_signal(i)).*abs(tx_signal(i));
    power_receive = power_receive + abs(rx_signal_receive(i)).*abs(rx_signal_receive(i));
end
power_transmit = power_transmit/length(time);
power_receive = power_receive/length(time);

SNR = power_receive/(noise_var*noise_var);
SNR_dB = 10*log10(SNR);
path_loss = power_transmit/power_receive;
path_loss_dB= 10*log10(path_loss);

delay_mean = sum(ray_matrix(:,1).*ray_matrix(:,2))/sum(ray_matrix(:,2));
delay_rms = sqrt(sum(ray_matrix(:,2).*(ray_matrix(:,1)-delay_mean).^2)/sum(ray_matrix(:,2)));

fprintf('Signal-to-noise ratio:\t%.2f\tdB\n',SNR_dB)
fprintf('Path loss (Tx to Rx):\t%.2f\tdB\n',path_loss_dB)
fprintf('Delay spread:\t\t\t%.3f\tns\n',1e9*delay_rms)