clear; close all;

%% Initial conditions and definitions
xbounds = [0 3];
ybounds = [0 3];
ship_start = 1;
dt = 5e-7;
start_time = 0;
end_time = 0.02;

%Sim = matfile('Sim.mat');

tlen = round((end_time-start_time)/dt);             % number of timesteps
N = 1e5;                                          % number of particles
r_t = zeros(N,2,2);                        % init 3D mat for particle locs
rands = rand(N,1)*2*pi;
V_inf = 7800;
V_temp = 990;


r_t(:,:,1) = [rand(N,1)*(ship_start+1)-1,rand(N,1)*ybounds(2)];      % initial positions
V = [V_inf+V_temp*cos(rands),V_temp*sin(rands)];    % initial velocities

N_total = 0;
dVx = 0;

clearvars("rands")
%% Set wall conditions
% Set wall coordinates
wall = [0 0;
        0.3 0.15;
        0.5 0.2;
        0.6 0.2;
        0.9 0.45;
        1.2 0.45;
        1.2 0.2;
        1.3 0.2;
        1.7 0.1;
        1.7 -0.1;
        1.3 -0.2;
        0.5 -0.2;
        0.3 -0.15;
        0 0];
%wall(:,2) = wall(:,2) - 2;
num_walls = length(wall) - 1;                       % number of walls
wall = wall + ones(num_walls + 1,2) .* [ship_start, 1.5];
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


last_rxw = zeros(N,num_walls);                    % allocate mem for last_rxw
rxw = zeros(N,num_walls);                         % allocate mem for rxw
%last_rxw2 = zeros(N,3,num_walls);

for wi = 1:num_walls
    last_rxw(:,wi) = (r_t(:,1,1) - wall(wi,1)) .* wall_cross_matrix(:,2,wi) ...
                      - (r_t(:,2,1) - wall(wi,2)) .* wall_cross_matrix(:,1,wi);
    % last_rxw(:,:,wi) = cross([r_t(:,:,1)-(ones(N,2).*wall(wi,:)), ...
    %             zeros(N,1)],wall_cross_matrix(:,:,wi));
                                % cross product matrix between particles
                                % and walls
end

%% Set up rendering

figure
ax = axes;

real_time = 20;

l = tlen;
timesteps = round(l/2):round(0.5*l/(real_time*60)):l;
timesteps = [timesteps, -1];

delete('M.mat')
M = matfile('M.mat','Writable',true);

M.movie(1,length(timesteps)-1) = struct('cdata',[],'colormap',[]);

plot(ax,wall(:,1),wall(:,2),'LineWidth',1)
set(ax,'XLim',xbounds,'YLim',ybounds)
grid on
set(gcf,"Position",[400 200 144*8 128*8],'Resize','off')
h = animatedline(ax,'Color',[0 0 0],'LineStyle','none','Marker','.','MarkerSize',1);

render_counter = 1;

%% Simulate for timesteps

for i = 1:tlen
    % Ignoring gravity so position is just r + V*dt
    r_t(:,:,2) = r_t(:,:,1) + dt*V;

    % calculate current rxw for wall collision detection
    for wi = 1:num_walls
        rxw(:,wi) = (r_t(:,1,2) - wall(wi,1)) .* wall_cross_matrix(:,2,wi) ...
                      - (r_t(:,2,2) - wall(wi,2)) .* wall_cross_matrix(:,1,wi);
        % rxw(:,:,wi) = cross([r_t(:,:,2) - (ones(N,2) .* ...
        %     wall(wi,:)),zeros(N,1)],wall_cross_matrix(:,:,wi));
    end

    % iterate over every particle
    for part = 1:N
        % iterate over every wall
        for wi = 1:num_walls
            % check for whether the particle has crossed the current
            % wall in the last timestep
            if rxw(part,wi)*last_rxw(part,wi) <= 0 && bounced(part) == 0 ...
                    && ((r_t(part,1,2) <= max(wall(wi:wi+1,1))  ...
                    &&   r_t(part,1,2) >= min(wall(wi:wi+1,1))) ...
                    ||  (r_t(part,2,2) <= max(wall(wi:wi+1,2))  ...
                    &&   r_t(part,2,2) >= min(wall(wi:wi+1,2))) ...
                    ||  (r_t(part,1,1) <= max(wall(wi:wi+1,1))  ...
                    &&   r_t(part,1,1) >= min(wall(wi:wi+1,1))) ...
                    ||  (r_t(part,2,1) <= max(wall(wi:wi+1,2))  ...
                    &&   r_t(part,2,1) >= min(wall(wi:wi+1,2))))

                % if it has, check the angle between the particle and
                % the wall
                theta = asin(norm(cross([V(part,:),0],wall_cross_matrix(1,:,wi))) ...
                           /(norm(V(part,:))*norm(wall_cross_matrix(1,:,wi))));

                % do some magic to make sure the angle is right (I kind
                % of grinded plotting this stuff in desmos until my
                % math and logic were correct and things worked right)
                if rxw(part,wi) >= 0
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
                theta = theta + normrnd(0,0.2);

                % since this is a little bit of a hacky way to
                % implement collision, the only change to the particle
                % is not a position update but instead a change in
                % velocity's angle according to how it bounced.
                dVx = dVx - V(part,1);
                V(part,:) = [V(part,1)*cos(theta)-V(part,2)*sin(theta), ...
                          V(part,1)*sin(theta)+V(part,2)*cos(theta)];
                dVx = dVx + V(part,1);
                N_total = N_total + 1;
                r_t(part,:,2) = r_t(part,:,1);
                rxw(part,:) = last_rxw(part,:);
                % if it did bounce, set the 'bounced' variable for the
                % particle to 3 so it can't bounce for another 3
                % timesteps (to make sure it doesn't keep crossing the
                % same wall and thinks it's bouncing every time it
                % crosses it).
                bounced(part) = 0;
            end
        end
        % logic for resetting the particles if it exceeds bounds of the
        % simulation
        if r_t(part,1,2) >= max(xbounds) + 0.125...
        || r_t(part,1,2) <= min(xbounds) - 1 ...
        || r_t(part,2,2) >= max(ybounds) + 0.125...
        || r_t(part,2,2) <= min(ybounds) - 0.125
            r_t(part,:,2) = [rand(1,1)-1, rand(1,1)*ybounds(2)];
            randomnum = rand(1,1)*2*pi;
            V(part,:) = [V_inf+V_temp*cos(randomnum),V_temp*sin(randomnum)];
            bounced(part) = 2;
        end
        % if the particle's bounce counter is non-zero, decrement
        % the bounce counter so it can eventually bounce again.
        if bounced(part) > 0
            bounced(part) = bounced(part) - 1;
        end
    end
    if i == timesteps(render_counter)
        clearpoints(h)
        addpoints(h,r_t(:,1,2),r_t(:,2,2))
        drawnow

        M.movie(1,render_counter) = getframe(gcf);
        render_counter = render_counter + 1;
    end
    last_rxw = rxw;
    r_t(:,:,1) = r_t(:,:,2);

end
%        _
%  _____/ |__
% /         |
% \_________|

%% Write to video

% Clear up memory; data has already been calculated and saved to disk.
clearvars("last_rxw","rxw","V","wall_cross_matrix","bounced","r_t")


vid = VideoWriter('simandrender','MPEG-4');
vid.FrameRate = 60;
open(vid)
writeVideo(vid,M.movie)
close(vid)

rho = 1.06042541856925e-9;
mAv = 2.65686251e-26;
emptyarea = (xbounds(2) - xbounds(1) + 1.125)*(ybounds(2) - ybounds(1) + 0.25) - 0.6675;
dP = dVx * mAv% * N_total;
num_rho = N/emptyarea;
sim_rho = num_rho * mAv;
ratio = rho/sim_rho;
F = dP * ratio / (end_time - start_time)


