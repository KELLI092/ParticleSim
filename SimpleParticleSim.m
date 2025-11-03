clear; close all;

t = 0:0.000001:0.015;
N = 300;
r_0 = [rand(N/2,1)*2,rand(N/2,1)*14+0.5;
       rand(N/2,1)*2+13,rand(N/2,1)*14+0.5];
V = [rand(N/2,1)*2+1,rand(N/2,1)-0.5;
     rand(N/2,1)*-2-1,rand(N/2,1)-0.5]*1e3;

num_walls = 2;
wall = [5,5;10,10];
wall_cross_matrix = ones(N,3);
wall_cross_matrix(:,1) = wall(2,1) - wall(1,1);
wall_cross_matrix(:,2) = wall(2,2) - wall(1,2);
wall_angle = atan((wall(2,2)-wall(1,2))/(wall(2,1)-wall(1,1)));
bounced = zeros(N,1);

r_t = zeros(N,width(r_0),length(t));

last_t = 0;
last_rxw = cross([r_0-(ones(N,2).*wall(1,:)),zeros(N,1)],wall_cross_matrix);
for i = 1:length(t)
    if i == 1
        r_t(:,:,i) = r_0;
    else
        r_t(:,:,i) = r_t(:,:,i-1) + (t(i)-last_t)*V;

        rxw = cross([r_t(:,:,i)-(ones(N,2).*wall(1,:)),zeros(N,1)], ...
            wall_cross_matrix);
        for j = 1:N
            if rxw(j,3)*last_rxw(j,3) <= 0 && bounced(j) == 0 ...
                    && r_t(j,1,i) <= max(wall(:,1)) ...
                    && r_t(j,1,i) >= min(wall(:,1)) ...
                    && r_t(j,2,i) <= max(wall(:,2)) ...
                    && r_t(j,2,i) >= min(wall(:,2))

                theta = asin(norm(cross([V(j,:),0],wall_cross_matrix(1,:))) ...
                           /(norm(V(j,:))*norm(wall_cross_matrix(1,:))));
                if rxw(j,3) >= 0
                    if dot([V(j,:),0],wall_cross_matrix(1,:)) < 0
                        theta = -2*theta;
                    else
                        theta = 2*theta;
                    end
                else
                    if dot([V(j,:),0],wall_cross_matrix(1,:)) > 0
                        theta = -2*theta;
                    else
                        theta = 2*theta;
                    end
                end
                V(j,:) = [V(j,1)*cos(theta)-V(j,2)*sin(theta), ...
                          V(j,1)*sin(theta)+V(j,2)*cos(theta)];
                bounced(j) = 1;
            end
        end
        last_rxw = rxw;
    end

    last_t = t(i);
    
end
