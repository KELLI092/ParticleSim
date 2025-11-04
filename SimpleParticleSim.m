clear; close all;

%% Initial conditions and definitions
xbounds = [0 15];
ybounds = [0 15];
dt = 0.000001;
t = 0:dt:0.05;                                     % number of timesteps
N = 1000;                                            % number of particles
r_0 = [rand(N,1)*4-2,rand(N,1)*14+0.5];        % initial positions
V = [rand(N,1)*2+1,rand(N,1)-0.5]*1e3;         % initial velocities

%r_0 = [11 9];
%V = 1e3*[0 1];

num_walls = 3;                                      % number of walls
wall = [5,5;5,10;10,10;10,5];                           % wall coordinates;
                                % creates a connected length of walls
wall_cross_matrix = zeros(N,3,num_walls);           % setup wall vector mat
wall_angle = zeros(num_walls,1);                    % setup wall angle mat

for i = 1:num_walls                                 % populate angle and cross mats
    wall_angle(i) = atan((wall(i+1,2)-wall(i,2)) ...
                         / (wall(i+1,1)-wall(i,1)));
    wall_cross_matrix(:,1,i) = wall(i+1,1) - wall(i,1);
    wall_cross_matrix(:,2,i) = wall(i+1,2) - wall(i,2);
end

bounced = zeros(N,1);                               % recently bounced check mat

r_t = zeros(N,width(r_0),length(t));                % init 3D mat for particle locs

last_rxw = zeros(N,3,num_walls);                    % allocate mem for last_rxw
for wi = 1:num_walls
    last_rxw(:,:,wi) = cross([r_0-(ones(N,2).*wall(wi,:)), ...
                zeros(N,1)],wall_cross_matrix(:,:,wi));
                                % cross product matrix between particles
                                % and walls
end

for i = 1:length(t)
    if i == 1
        r_t(:,:,i) = r_0;
    else
        r_t(:,:,i) = r_t(:,:,i-1) + dt*V;

        rxw = zeros(N,3,num_walls);                 % allocate mem for rxw
        for wi = 1:num_walls
            rxw(:,:,wi) = cross([r_t(:,:,i) - (ones(N,2) .* ...
                wall(wi,:)),zeros(N,1)],wall_cross_matrix(:,:,wi));
        end

        for part = 1:N
            for wi = 1:num_walls
                if rxw(part,3,wi)*last_rxw(part,3,wi) <= 0 && bounced(part) == 0 ...
                        && ((r_t(part,1,i) <= max(wall(wi:wi+1,1)) ...
                        && r_t(part,1,i) >= min(wall(wi:wi+1,1))) ...
                        || (r_t(part,2,i) <= max(wall(wi:wi+1,2)) ...
                        && r_t(part,2,i) >= min(wall(wi:wi+1,2))))
    
                    theta = asin(norm(cross([V(part,:),0],wall_cross_matrix(1,:,wi))) ...
                               /(norm(V(part,:))*norm(wall_cross_matrix(1,:,wi))));
                    if rxw(part,3,wi) >= 0
                        if dot([V(part,:),0],wall_cross_matrix(1,:,wi)) < 0
                            theta = -2*theta;
                        else
                            theta = 2*theta;
                        end
                    else
                        if dot([V(part,:),0],wall_cross_matrix(1,:,wi)) > 0
                            theta = -2*theta;
                        else
                            theta = 2*theta;
                        end
                    end
                    V(part,:) = [V(part,1)*cos(theta)-V(part,2)*sin(theta), ...
                              V(part,1)*sin(theta)+V(part,2)*cos(theta)];
                    bounced(part) = 3;
                elseif bounced(part) ~= 0
                    bounced(part) = bounced(part) - 1;
                end
            end
            if r_t(part,1,i) >= max(xbounds) ...
            || r_t(part,1,i) <= min(xbounds) - 2 ...
            || r_t(part,2,i) >= max(ybounds) ...
            || r_t(part,2,i) <= min(ybounds)
                r_t(part,:,i) = [rand(1,1)*2-2,rand(1,1)*14+0.5];
                V(part,:) = [rand(1,1)*2+1,rand(1,1)-0.5]*1e3;
            end
        end
        last_rxw = rxw;
    end
    
end
