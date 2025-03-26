%% Molecular Dynamics
% Lennard Jones potential (12-6)
% Periodic boundaries
% NEEDS: circles(...), slanCM(...), gif(...)

clear;clc;close all;

fs=30;      % Font size
ls= 20;     % Label size

% Time discretization
dt = 0.001; t0 = 0; tf = 10; 
t = t0:dt:tf; nt = length(t);

N = 50;             % Number of Particles
Lx = 30;            % Width of box
Ly = Lx;            % Height of box
D = 1;              % Diamter of the particles
eps = 5;            % Energy of interactions    

m = 1; gamma = 10; k = 1; % Mass, viscosity, Boltzmann Constant  
T = 100;                  % Temperature
g = sqrt(2*gamma*k*T);    % Varaiance of Langevin Thermostat to make it coherent with the Cannonical Ensemble

% Constants used for Velocity Verlet
a = (1 - gamma*dt/(2*m))/(1 + gamma*dt/(2*m));
b = 1/(1 + gamma*dt/(2*m));

% Matrices to store the positions and velocities of all particles through
% time
x = zeros(N, nt);
y = zeros(N, nt);
vx = zeros(N, nt);
vy = zeros(N, nt);
Fx = zeros(N, nt);
Fy = zeros(N, nt);

% Used my Simulated Annealing algorithm to make sure the particles don't
% start to close to each other
[x0, y0] = packing(N, D/2, Lx, Ly);

% Matrix used to store the particles and the moments in time they jump from
% one side of the box to the other, to make trajectories look as one
stops = logical(zeros(N, nt));

%%
tic % Sart timer. With N=50, T=100, it takes approx 20s

% Inital velocities
vy0 = g*sqrt(dt)*randn(N, 1);
vx0 = g*sqrt(dt)*randn(N, 1);

% Store initial positions
x(:, 1) = x0;
y(:, 1) = y0;
vx(:, 1) = vx0;
vy(:, 1) = vy0;

Px = sum(vx0);
Py = sum(vy0);
P = Px + Py;
vx0 = vx0-P/N;

% New position and velocity functions
np = @(x, v, dW, F) x + dt*b*v + b*dt^2/(2*m)*F + b*dt/(2*m)*g*dW;
nv = @(v, dW, F1, F2) a*v + (dt/(2*m))*(a*F1 + F2) + (b/m)*g*dW;

% Get all the forces sensed by each particle
[Fxi1, Fyi1] = forces(x(:, 1), y(:, 1), N, D, eps, Lx, Ly);

% Store initial values of force
Fx(:, 1) = Fxi1;
Fy(:, 1) = Fyi1;

% Simulate through time
for i = 1:nt-1
    
    pacman = false;                 % Variable that checks if a particle
                                    % jumps from one side to the other

    dWx = sqrt(dt)*randn(N, 1);     % Random normal increment
    dWy = sqrt(dt)*randn(N, 1);
    
    % Forces in time t_i
    Fxi = Fxi1;
    Fyi = Fyi1;
    
    % Get positions in time t_(i+1)
    x(:, i+1) = np(x(:, i), vx(:, i), dWx, Fxi);
    y(:, i+1) = np(y(:, i), vy(:, i), dWy, Fyi);

    % Check boundaries
    xt = x(:, i+1);
    yt = y(:, i+1);

    xd = xt>=Lx;
    if ~isempty(xt(xd))
        x(xd, i+1) = xt(xd) - Lx;
        pacman = true;
    end
    xi = xt<0;
    if ~isempty(xt(xi))
        x(xi, i+1) = xt(xi) + Lx;
        pacman = true;
    end

    yd = yt>=Ly;
    if ~isempty(yt(yd))
        y(yd, i+1) = yt(yd) - Ly;
        pacman =true;
    end
    yi = yt<0;
    if ~isempty(yt(yi))
        y(yi, i+1) = yt(yi) + Ly;
        pacman =true;
    end
    
    if pacman
        yc1 = find(yd);
        yc2 = find(yi);
        xc1 = find(xd);
        xc2 = find(xi);
        
        % Saves index of particles that jumped
        ch = unique([yc1;yc2;xc1;xc2]);
        
        stops(ch, i+1) = true;
    end

    xi1 = x(:, i+1);
    yi1 = y(:, i+1);
    
    % Get forces in time t_(i+1)
    [Fxi1, Fyi1] = forces(xi1, yi1, N, D, eps, Lx, Ly);
    
    % Save forces
    Fx(:, i+1) = Fxi1;
    Fy(:, i+1) = Fyi1;
    
    % Get velocities in time t_(i+1)
    vx(:, i+1) = nv(vx(:, i), dWx, Fxi, Fxi1);
    vy(:, i+1) = nv(vy(:, i), dWy, Fyi, Fyi1);
end
toc     % Stop timer
%% Get times in which jumps occur
v = sqrt(vx.^2 + vy.^2);
nvx = vx./v;
nvy = vy./v;

steps = ones(N, 1);
sc = 1;
ts = [];
for i = 1:1:nt
        
    sec = steps(:, sc);
    soi = stops(:, i);

    if sum(soi)
        sc = sc+1;
        ts = [ts, i];
        sec(soi) = i;
        steps = [steps, sec];
    end
end
ts = [1, ts];
%% ANIMATION
clc;close all;
figure;
tic
R = D*ones(N, 1)/2;
colors = slanCM('glasbey_dark', N);

tray = 0;   % Show trajectories 0-False, 1-True (Slow)
imag = 0;   % Show just last frame of simulation 0-False, 1-True

if imag == 1    
    tin = nt;
else
    tin = 1;
end

for i = tin:15:nt
    
    sc = find(ts<=i, 1, "last");
    
    clf
    hold on;
    % Draw trajectories if tray=1
    if tray
        for j = 1:N
            stepsj = [unique(steps(j, 1:sc)), i];
            scj = length(stepsj);
            for k = 1:scj-1
                startk = stepsj(k);
                stopk = stepsj(k+1)-1;
                plot(x(j, startk:stopk), y(j, startk:stopk), 'Color', colors(j, :), LineWidth=1)
            end
        end
    end
    circles(x(:, i), y(:, i), R, 'FaceColor', zeros(N, 3), 'EdgeColor', 'none');

    ipx = x(:, i)<D/2;
    imx = x(:, i)>Lx-D/2;
    ipy = y(:, i)<D/2;
    imy = y(:, i)>Ly-D/2;
    if ~isempty(x(ipx, i))
        circles(x(ipx, i)+Lx, y(ipx, i), R(ipx), 'FaceColor', zeros(sum(ipx), 3), 'EdgeColor', 'none');
    end
    if ~isempty(x(imx, i))
        circles(x(imx, i)-Lx, y(imx, i), R(imx), 'FaceColor', zeros(sum(imx), 3), 'EdgeColor', 'none');
    end
    if ~isempty(y(ipy, i))
        circles(x(ipy, i), y(ipy, i)+Ly, R(ipy), 'FaceColor', zeros(sum(ipy), 3), 'EdgeColor', 'none');
    end
    if ~isempty(y(imy, i))
        circles(x(imy, i), y(imy, i)-Ly, R(imy), 'FaceColor', zeros(sum(imy), 3), 'EdgeColor', 'none');
    end

    negs = (x(:, i)>Lx-D/2 | x(:, i)<D/2) & (y(:, i)>Ly-D/2 | y(:, i)<D/2);
    if ~isempty(x(negs, i))
        circles(x(negs, i)+sign(Lx/2-x(negs, i))*Lx, y(negs, i)+sign(Ly/2-y(negs, i))*Ly, R(negs), 'FaceColor', zeros(sum(negs), 3), 'EdgeColor', 'none');
    end

    % Uncomment to show direction and magnitude of forces
    % quiver(x(:, i), y(:, i), Fx(:, i), Fy(:, i))
    % Uncomment to show direction of movement
    % quiver(x(:, i), y(:, i), nvx(:, i), nvy(:, i), 0.5)
    hold off
    
    axis equal
    xlim([0, Lx])
    ylim([0, Ly])
    xlabel('$x$', FontSize=fs)
    ylh = ylabel('$y$', FontSize=fs, Rotation=0);
    ylh.Position(1)=ylh.Position(1);

    % Uncomment if want to save a gif of the animation
    % if i == 1
    %     gif('GAS.gif','DelayTime', 1/60)
    % else
    %     gif
    % end
    % drawnow
    pause(1/60)
    
end
toc
%% IMAGES
colors = slanCM('glasbey_dark', N);%rand(N, 3);
clc;close all;
figure;
tiledlayout(2, 2, "TileSpacing","tight", "Padding","tight")
tic
R = D*ones(N, 1)/2;

tis = [1, ceil(0.25*nt), ceil(0.5*nt), nt];
tray = 1;               % Do not change
imag = 0;

if imag == 1
    tin = nt;
else
    tin = 1;
end
for tfi = tis
    nexttile;
    tin = tfi;
    for i = tin:15:tfi 
        sc = find(ts<=i, 1, "last");
        hold on;
        if tray
            for j = 1:N
                stepsj = [unique(steps(j, 1:sc)), i];
                scj = length(stepsj);
                for k = 1:scj-1
                    startk = stepsj(k);
                    stopk = stepsj(k+1)-1;
                    plot(x(j, startk:stopk), y(j, startk:stopk), 'Color', [colors(j, :), 0.5], LineWidth=1)
                end
            end
        else
            clf
        end
        circles(x(:, i), y(:, i), R, 'FaceColor', colors, 'EdgeColor', 'k', 'LineWidth', 0.1);

        ipx = x(:, i)<D/2;
        imx = x(:, i)>Lx-D/2;
        ipy = y(:, i)<D/2;
        imy = y(:, i)>Ly-D/2;
        if ~isempty(x(ipx, i))
            circles(x(ipx, i)+Lx, y(ipx, i), R(ipx), 'FaceColor', colors(ipx, :), 'EdgeColor', 'none');
        end
        if ~isempty(x(imx, i))
            circles(x(imx, i)-Lx, y(imx, i), R(imx), 'FaceColor', colors(imx, :), 'EdgeColor', 'none');
        end
        if ~isempty(y(ipy, i))
            circles(x(ipy, i), y(ipy, i)+Ly, R(ipy), 'FaceColor', colors(ipy, :), 'EdgeColor', 'none');
        end
        if ~isempty(y(imy, i))
            circles(x(imy, i), y(imy, i)-Ly, R(imy), 'FaceColor', colors(imy, :), 'EdgeColor', 'none');
        end

        negs = (x(:, i)>Lx-D/2 | x(:, i)<D/2) & (y(:, i)>Ly-D/2 | y(:, i)<D/2);
        if ~isempty(x(negs, i))
            circles(x(negs, i)+sign(Lx/2-x(negs, i))*Lx, y(negs, i)+sign(Ly/2-y(negs, i))*Ly, R(negs), 'FaceColor', colors(negs, :), 'EdgeColor', 'none');
        end


        % quiver(x(:, i), y(:, i), Fx(:, i), Fy(:, i))
        % quiver(x(:, i), y(:, i), nvx(:, i), nvy(:, i), 0.5)

        rectangle('Position',[0, 0, Lx, Ly])
        hold off

        axis equal
        % xlim([-1, Lx+1])
        % ylim([-1, Ly+1])
        xlim([0, Lx])
        ylim([0, Ly])
        xticks([0, Lx])
        yticks([0, Ly])
        % xlabel('$x$', FontSize=fs)
        % ylh = ylabel('$y$', FontSize=fs, Rotation=0);
        % ylh.Position(1)=ylh.Position(1);
        drawnow
        pause(dt*1)
        title("$t = " + string(t(tin)) + "$")
    end
end
toc
%% KINETIC ENERGY
EK = sum((m/2)*(vx.^2 + vy.^2));
figure;
plot(t, EK, LineWidth=1, HandleVisibility="off")
yline(N*k*T, LineWidth=1.5, DisplayName="$K=NkT$")  
xlabel('$x$', FontSize=fs)
ylabel('$\langle K \rangle$', FontSize=fs, Rotation=0)
legend(FontSize=ls)
axis square














% Functions
function [Fx, Fy] = forces(x, y, N, D, eps, Lx, Ly)

Fdx = @(x, r) 48*eps*(D^6./r.^8).*(D^6./r.^6 - 1/2).*x;
Fdy = @(y, r) 48*eps*(D^6./r.^8).*(D^6./r.^6 - 1/2).*y;

Fx = zeros(N, 1);
Fy = zeros(N, 1);

for i = 1:N

    xip = x(i)-x > Lx/2;
    xim = x(i)-x < -Lx/2;
    xn = x;
    if ~isempty(x(xip))
        xn(xip) = x(xip)+Lx;
    end
    if ~isempty(x(xim))
        xn(xim) = x(xim)-Lx;
    end

    yip = y(i)-y > Ly/2;
    yim = y(i)-y < -Ly/2;
    yn = y;
    if ~isempty(y(yip))
        yn(yip) = y(yip)+Ly;
    end
    if ~isempty(y(yim))
        yn(yim) = y(yim)-Ly;
    end
    
    dx = xn(i) - xn;
    dy = yn(i) - yn;
    dr = sqrt(dx.^2 + dy.^2);


    Fxi = Fdx(dx, dr);
    Fyi = Fdy(dy, dr);

    Fxi(i) = 0;
    Fyi(i) = 0;

    Fx(i) = sum(Fxi);
    Fy(i) = sum(Fyi);
end
end






function [x0, y0] = inPos(N, Lx, Ly, D)

x = Lx*rand(N, 1);
y = Ly*rand(N, 1);

exit = false;

while ~exit
    for i = 1:N
        xip = x(i)-x > Lx/2;
        xim = x(i)-x < -Lx/2;
        xn = x;
        if ~isempty(x(xip))
            xn(xip) = x(xip)+Lx;
        end
        if ~isempty(x(xip))
            xn(xim) = x(xim)-Lx;
        end

        yip = y(i)-y > Ly/2;
        yim = y(i)-y < -Ly/2;
        yn = y;
        if ~isempty(y(yip))
            yn(yip) = y(yip)+Ly;
        end
        if ~isempty(y(yim))
            yn(yim) = y(yim)-Ly;
        end

        dx = xn(i) - xn;
        dy = yn(i) - yn;
    
        r = sqrt(dx.^2 + dy.^2);
        r(i) = 13*D;
        ril = r<D;
        rl = r(ril);

        if ~isempty(rl)
            xl = x(ril);
            xl = Lx*rand(length(xl), 1);
            x(ril) = xl;
            yl = y(ril);
            yl = Ly*rand(length(yl), 1);
            y(ril) = yl;
        end
    end

    check = true;
    for i = 1:N
        xip = x(i)-x > Lx/2;
        xim = x(i)-x < -Lx/2;
        xn = x;
        if ~isempty(x(xip))
            xn(xip) = x(xip)+Lx;
        end
        if ~isempty(x(xip))
            xn(xim) = x(xim)-Lx;
        end

        yip = y(i)-y > Ly/2;
        yim = y(i)-y < -Ly/2;
        yn = y;
        if ~isempty(y(yip))
            yn(yip) = y(yip)+Ly;
        end
        if ~isempty(y(yim))
            yn(yim) = y(yim)-Ly;
        end

        dx = xn(i) - xn;
        dy = yn(i) - yn;
    
        r = sqrt(dx.^2 + dy.^2);
        r(i) = 13*D;
        ril = r<D;
        rl = r(ril);
        if ~isempty(rl)
            check = false;
            break
        end
    end
    clf
    viscircles([x, y],  D/2);
    xlim([0, Lx])
    ylim([0, Ly])
    drawnow
    if check
        exit = true;
    end    
end

x0 = x;
y0 = y;

end

function [XX, YY] = packing(N, r1, Lx, Ly) 

r2 = 0.8;

groups = r1;

nr1 = 1;
nr2 = 1-nr1;

disks = zeros(N, 4);

disks(:, 1) = (Lx*rand(1, N))';
disks(:, 2) = (Ly*rand(1, N))';
disks(:, 3) = sqrt(disks(:, 1).^2 + disks(:, 2).^2);
disks(:, 4) = [ones(1, round(nr1*N))*r1, ones(1, round(nr2*N))*r2]';

itLim = 10000;
T0 = 100;
U = 1e5;
T = T0;
P=0;

found = false;



% figure;
while ~found
    for it = 1:itLim

        i = randi([1, N]);
        xi = disks(i, 1);
        yi = disks(i, 2);

        [dx, dy] = change(disks, i, Lx, Ly, N);
        if isnan(dx) || isnan(dy)
            break
        end
        nx = xi + dx;
        ny = yi + dy;
        
        if yi+dy == yi && dy~=0
            ny = yi + sign(dy)*1e-3;
            % ny==yi
        end
        if xi+dx == xi && dx~=0
            nx = xi + sign(dx)*1e-3;
            % nx==xi
        end

        nr = sqrt(nx^2 + ny^2);
        nR = disks(i, 4);
        nrow = [nx, ny, nr, nR];

        ndisks = disks;
        ndisks(i, :) = nrow;

        U = energy(disks, Lx, Ly, N);
        Un = energy(ndisks, Lx, Ly, N);
       
        dU = Un-U;

        if dU <= 0
            disks = ndisks;
        else
            T = temperature(it, itLim, T0);
            % if U<1e-6   
                P = exp(-dU/T);
            % else
            %     P = exp(-dU/(U*T));
            % end
            rn = rand();
            clc
            if rn < P
                disks = ndisks;
            end
        end
    end
    if U < 1e-2
        found = true;
        break
    end
    ichange = 0;
    if U > 1e-2
        g = 1;
        while g <= length(groups)
            group = disks(disks(:, 4) == 0.8, :);
            ng = length(group);
            get = false;
            for gi = 1:ng
                ichange = groupEnergy(group, Lx, Ly, ng);
                if ichange ~= 0
                    
                    dchange = group(ichange, :);
                    tochange = find(ismember(disks, dchange, 'rows'));
                    xnew = Lx*rand();
                    ynew = Ly*rand();
                    disks(tochange, :) = [xnew, ynew, sqrt(xnew^2 + ynew^2), dchange(4)];
                    get = true;
                    break;
                end
            end
            if get
                break
            end
            g = g+1;
        end
    end
    if U<1e-2
        ichange=100;
    end
    if ichange == 0
        Uij = zeros(1, N);
        n = zeros(1, N);
        for i = 1:N
            for j = 1:N
                if i == j
                    continue
                end
                dij = embdepth(disks, i, j, Lx, Ly);
                if dij > 0
                    n(i) = n(i) + 1;
                end
                Uij(i) = Uij(i) + dij^2;
            end
        end
        ichange = find(Uij==max(Uij));
        dchange = disks(ichange, :);
        xnew = Lx*rand();
        ynew = Ly*rand();
        disks(ichange, :) = [xnew, ynew, sqrt(xnew^2 + ynew^2), dchange(4)];
    end
end
XX = disks(:, 1);
YY = disks(:, 2);
end


function i = groupEnergy(D, Lx, Ly, N)

U = 0;
for i = 1:N
    for j = 0:N
        if i == j
            continue
        end
        dij = embdepth(D, i, j, Lx, Ly);

        U = U + dij^2;
        
        if U > 0
            return
        end
    end
end

i = 0;
return
end

function U = energy(D, Lx, Ly, N)

U = 0;
for i = 1:N
    for j = 0:N
        if i == j
            continue
        end
        dij = embdepth(D, i, j, Lx, Ly);

        U = U + dij^2;

    end
end
end


function [dx, dy] = change(D, i, Lx, Ly, N)
dx = 0;
dy = 0;
for j = 0:N
    
    if i == j
        continue
    end
    
    if j ~= 0
        dij = embdepth(D, i, j, Lx, Ly);
    end
    Dij = dist(D, i, j);

    % if Dij < 1e-1
    %     Dij = 1e-1;
    % end

    if j == 0
        if D(i, 1) + D(i, 4) > Lx
            dxij = -D(i, 1) - D(i, 4) + Lx;
        elseif D(i, 1) - D(i, 4) < 0
            dxij = -D(i, 1) + D(i, 4);
        else
            dxij = 0;
        end
        if D(i, 2) + D(i, 4) > Ly
            dyij = -D(i, 2) - D(i, 4) + Ly;
        elseif D(i, 2) - D(i, 4) < 0
            dyij = -D(i, 2) + D(i, 4);
        else
            dyij = 0;
        end
        dij = sqrt(dxij^2 + dyij^2);
    else
        
        % if D(i, 1) - D(j, 1) == 0
        %     dxij = (0.01)*dij/Dij;
        % else
            dxij = (D(i, 1) - D(j, 1))*dij/Dij;
        % end
        dyij = (D(i, 2) - D(j, 2))*dij/Dij;
    end
    % fprintf("i = " + string(i)+", j = " +string(j)+ "  dx = " + string(dxij) + ", dy = " + string(dyij) + "  dij= "+string(dij)+"\n")
    dx = dx + dxij;
    dy = dy + dyij;
end
% fprintf("\n")
end

function Dij = dist(D, i, j)

if j == 0
    Dij = abs(D(i, 3));
else
    Dij = sqrt( (D(i, 1) - D(j, 1))^2 + (D(i, 2) - D(j, 2))^2   );
end
end



function dij = embdepth(D, i, j, Lx, Ly)

if j == 0
    if D(i, 1) + D(i, 4) > Lx
        dijx = D(i, 1) + D(i, 4) - Lx;
    elseif D(i, 1) - D(i, 4) < 0
        dijx = D(i, 1) - D(i, 4);
    else
        dijx = 0;
    end
    if D(i, 2) + D(i, 4) > Ly
        dijy = D(i, 2) + D(i, 4) - Ly;
    elseif D(i, 2) - D(i, 4) < 0
        dijy = D(i, 2) - D(i, 4);
    else
        dijy = 0;
    end
    dij = sqrt(dijx^2 + dijy^2);
    return
else
    if sqrt((D(j, 1)-D(i, 1)).^2 + (D(j, 2)-D(i, 2)).^2) < D(i, 4) + D(j, 4)
        dij = -sqrt((D(j, 1)-D(i, 1)).^2 + (D(j, 2)-D(i, 2)).^2) + D(i, 4) + D(j, 4);
        return
    else
        dij = 0;
        return
    end
end


end


function T = temperature(it, itMax, T0)

T = T0*0.9985^(it);

end
