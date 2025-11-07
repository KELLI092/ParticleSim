clear; close all;

%% Initial conditions and definitions
xbounds = [0 15];
ybounds = [0 15];
dt = 0.000001;
t = 0:dt:0.035;                                     % number of timesteps
N = 1000;                                            % number of particles
r_0 = [rand(N,1)*7-2,rand(N,1)*14+0.5];        % initial positions
V = [rand(N,1)+1.5,rand(N,1)-0.5]*1e3;         % initial velocities

%% Set wall conditions
% Set wall coordinates
wall = [7.5,6.5;9,7.5;7.5,8.5;12,8.5;12,6.5;7.5,6.5];

num_walls = length(wall) - 1;                       % number of walls
wall_cross_matrix = zeros(N,3,num_walls);           % setup wall vector mat
wall_angle = zeros(num_walls,1);                    % setup wall angle mat

for i = 1:num_walls                                 % populate angle and cross mats
    wall_angle(i) = atan((wall(i+1,2)-wall(i,2)) ...
                         / (wall(i+1,1)-wall(i,1)));
    wall_cross_matrix(:,1,i) = wall(i+1,1) - wall(i,1);
    wall_cross_matrix(:,2,i) = wall(i+1,2) - wall(i,2);
end

%% Initialize simulation necessary variables
bounced = zeros(N,1);                               % recently bounced check mat

r_t = zeros(N,width(r_0),length(t));                % init 3D mat for particle locs

last_rxw = zeros(N,3,num_walls);                    % allocate mem for last_rxw

for wi = 1:num_walls
    last_rxw(:,:,wi) = cross([r_0-(ones(N,2).*wall(wi,:)), ...
                zeros(N,1)],wall_cross_matrix(:,:,wi));
                                % cross product matrix between particles
                                % and walls
end

%% Simulate for timesteps

for i = 1:length(t)
    % if it is the first timestep, set the position to the initial position
    if i == 1
        r_t(:,:,i) = r_0;
    
    % otherwise, calculate the position
    else
        % Ignoring gravity so position is just r + V*dt
        r_t(:,:,i) = r_t(:,:,i-1) + dt*V;

        % calculate current rxw for wall collision detection
        rxw = zeros(N,3,num_walls);                 % allocate mem for rxw
        for wi = 1:num_walls
            rxw(:,:,wi) = cross([r_t(:,:,i) - (ones(N,2) .* ...
                wall(wi,:)),zeros(N,1)],wall_cross_matrix(:,:,wi));
        end

        % iterate over every particle
        for part = 1:N
            % iterate over every wall
            for wi = 1:num_walls
                % check for whether the particle has crossed the current
                % wall in the last timestep
                if rxw(part,3,wi)*last_rxw(part,3,wi) <= 0 && bounced(part) == 0 ...
                        && ((r_t(part,1,i) <= max(wall(wi:wi+1,1)) ...
                        && r_t(part,1,i) >= min(wall(wi:wi+1,1))) ...
                        || (r_t(part,2,i) <= max(wall(wi:wi+1,2)) ...
                        && r_t(part,2,i) >= min(wall(wi:wi+1,2))))
                    
                    % if it has, check the angle between the particle and
                    % the wall
                    theta = asin(norm(cross([V(part,:),0],wall_cross_matrix(1,:,wi))) ...
                               /(norm(V(part,:))*norm(wall_cross_matrix(1,:,wi))));

                    % do some magic to make sure the angle is right (I kind
                    % of grinded plotting this stuff in desmos until my
                    % math and logic were correct and things worked right)
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

                    % since this is a little bit of a hacky way to
                    % implement collision, the only change to the particle
                    % is not a position update but instead a change in
                    % velocity's angle according to how it bounced.
                    V(part,:) = [V(part,1)*cos(theta)-V(part,2)*sin(theta), ...
                              V(part,1)*sin(theta)+V(part,2)*cos(theta)];

                    % if it did bounce, set the 'bounced' variable for the
                    % particle to 3 so it can't bounce for another 3
                    % timesteps (to make sure it doesn't keep crossing the
                    % same wall and thinks it's bouncing every time it
                    % crosses it).
                    bounced(part) = 2;
                end
            end
            % logic for resetting the particles if it exceeds bounds of the
            % simulation
            if r_t(part,1,i) >= max(xbounds) + 0.25...
            || r_t(part,1,i) <= min(xbounds) - 2 ...
            || r_t(part,2,i) >= max(ybounds) + 0.25...
            || r_t(part,2,i) <= min(ybounds) - 0.25
                r_t(part,:,i) = [rand(1,1)*2-1, rand(1,1)*14+0.5];
                V(part,:) = [rand(1,1)+1.5,rand(1,1)/2-0.25]*1e3;
                bounced(part) = 2;
            end
            % if the particle's bounce counter is non-zero, decrement
            % the bounce counter so it can eventually bounce again.
            if bounced(part) > 0
                bounced(part) = bounced(part) - 1;
            end
        end
        last_rxw = rxw;
    end
    
end
