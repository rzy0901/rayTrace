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
figure_counter = 1;
figure(figure_counter)

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
                flag_occlude = 0
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
                        drawRay(A,B,ray_resolution)
                        start_points = [start_points; A];
                        end_points = [end_points; B];
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
                            drawRay(A,B,ray_resolution)
                            start_points = [start_points; A];
                            end_points = [end_points; B];
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
                            drawRay(A,B,ray_resolution)
                            start_points = [start_points; A];
                            end_points = [end_points; B];
                            ray_connection = [ray_connection;[A B]];
                            % Update path distance:
                            path_matrix( m,1) = path_matrix(m,1)+distance;
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
path_matrix(:,2) = path_matrix(:,1)/light_speed;
path_matrix(:,3) = (light_speed/(frequency_carrier*sqrt(4*pi)))./path_matrix(:,1);
path_matrix(:,5) = exp((-1j*2*pi*frequency_carrier/light_speed).*path_matrix(:,1));
% Get number of valid paths:
[paths_valid,~] = size(path_matrix);

% Initialise ray matrix:
ray_matrix = zeros(paths_valid,3);
% Populate ray matrix:
ray_matrix(:,1) = path_matrix(:,2);
ray_matrix(:,2) = path_matrix(:,3).*path_matrix(:,4);
ray_matrix(:,3) = path_matrix(:,5).*path_matrix(:,6);

% Plot power delay profile:
% figure(figure_counter+2)
% hold on
% plot(1e9*ray_matrix(:,1),10*log10(ray_matrix(:,2)),'x','MarkerSize',10)
% title('Power delay profile')
% xlabel('Time (ns)')
% ylabel('Received signal power (dB)')

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