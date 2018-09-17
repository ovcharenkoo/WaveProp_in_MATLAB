% Finite-differences in time domain(FDTD) acoustic wave propagation in 3D 
% isotropic medium surrounded by simple sponge boundaries with exponential 
% decay (Cerjan, 1985).
%
% We solve second order wave equation in time domain and displacement 
% formulation getting wavefield in terms of pressure field p.
%
% Acoustic medium is parametrized by acoustic velocity vp only.  We show 
% CFL number and number of points per wavelength prior running loop 
% over time steps.
% 
% Conventional FD star-stencils provides accuracy O(2,2)
% Stencils: dx2: [1 -2 1]/dx2
% --------------------------------------------------------------
% The code is intentionally writen in a single file
% to simplify start up.
%
% The program does not save any files, add such option manually if needed.
% Drawing the wavefield is the most computationally demanding. Increase
% IT_DISPLAY value to reduce output and accelerate computation.
%
% The goal is to provide a simple example of acoustic wave propagation
% in 3D isotropic medium.
%
% --------------------------------------------------------------
% Oleg Ovcharenko and Vladimir Kazei, 2018
%
% oleg.ovcharenko@kaust.edu.sa
% vladimir.kazei@kaust.edu.sa
%
% King Abdullah University of Science and Technology
% Thuwal, Saudi Arabia
% --------------------------------------------------------------

close all;
% Output periodicity in time steps
IT_DISPLAY = 20;

%% MODEL
% Model dimensions
nx = 101;
ny = 101;
nz = 101;
dx = 20;    % [m]
dy = 20;    % [m]
dz = 20;    % [m]

% Elastic parameters
vp = 3300.0 * ones(nz,ny,nx);       % velocity of compressional waves, [m/s]

%% TIME STEPPING
t_total = 0.35;                       % [sec] recording duration
dt = 0.5 * min(dx,dz)/max(vp(:));
nt = round(t_total/dt);             % number of time steps
t = [0:nt]*dt;

CFL = max(vp(:)) * dt / min(dx,dz);
dt2 = dt^2;
%% SOURCE
f0 = 10.0;                          % dominant frequency of the wavelet
t0 = 1.20 / f0;                     % half Ricker wavelet excitation time
factor = 1e8;                       % amplitude coefficient
angle_force = 90.0;                 % spatial orientation
                                    % (not relevant for acoustic case)
ksrc = round(nz/2);                 % source location along OZ
jsrc = round(ny/2);                 % source location along OZ
isrc = round(nx/2);                 % source location along OX

a = pi*pi*f0*f0;
source_term = factor * exp(-a*(t-t0).^2);                             % Gaussian
% source_term =  -factor*2.0*a*(t-t0)*exp(-a*(t-t0)^2);                % First derivative of a Gaussian:
% source_term = factor * (1.0 - 2.0*a*(t-t0).^2).*exp(-a*(t-t0).^2);        % Ricker source time function (second derivative of a Gaussian):
force_x = sin(angle_force * pi / 180) * source_term * dt2;

min_wavelengh = min(vp(vp>0.1))/f0;     % shortest wavelength bounded by velocity in the air

%% ABSORBING BOUNDARY (ABS)
abs_thick = min(floor(0.15*nx), floor(0.15*nz));    % thicknes of the layer
abs_rate = 0.3/abs_thick;                           % decay rate

lmargin = [abs_thick abs_thick abs_thick];
rmargin = [abs_thick abs_thick abs_thick];
weights = ones(nz+2,ny+2,nx+2);
for iz = 1:nz+2
    for iy = 1:ny+2
        for ix = 1:nx+2
            i = 0;
            j = 0;
            k = 0;
            if (ix < lmargin(1) + 1)
                i = lmargin(1) + 1 - ix;
            end
            if (iy < lmargin(2) + 1)
                j = lmargin(2) + 1 - iy;
            end
            if (iz < lmargin(3) + 1)
                k = lmargin(3) + 1 - iz;
            end
            if (nx - rmargin(1) < ix)
                i = ix - nx + rmargin(1);
            end
            if (ny - rmargin(2) < iy)
                j = iy - ny + rmargin(2);
            end
            if (nz - rmargin(3) < iz)
                k = iz - nz + rmargin(3);
            end
            if (i == 0 && j == 0 && k == 0)
                continue
            end
            rr = abs_rate * abs_rate * double(i*i + j*j + k*k);
            weights(iz,iy,ix) = exp(-rr);
        end
    end
end

%% SUMMARY
fprintf('#################################################\n');
fprintf('3D acoustic FDTD wave propagation in isotropic medium \nin displacement formulation with Cerjan(1985) \nboundary conditions\n');
fprintf('#################################################\n');
fprintf('Model:\n\t%d x %d x %d\tgrid nz x ny x nx\n\t%.1e x %.1e x %.1e\t[m] dz x dy x dx\n',nz,ny,nx,dz,dy,dx);
fprintf('\t%.1e x %.1e x %.1e\t[m] model size\n',nx*dx, ny*dy, nz*dz);
fprintf('\t%.1e...%.1e\t[m/s] vp\n', min(vp(:)), max(vp(:)));
fprintf('Time:\n\t%.1e\t[sec] total\n\t%.1e\t[sec] dt\n\t%d\ttime steps\n',t_total,dt,nt);
fprintf('Source:\n\t%.1e\t[Hz] dominant frequency\n\t%.1f\t[sec] index time\n',f0,t0);
fprintf('Other:\n\t%.1f\tCFL number\n', CFL);
fprintf('\t%.2f\t[m] shortest wavelength\n\t%d, %d\t points-per-wavelength OX, OZ\n', min_wavelengh, floor(min_wavelengh/dx), floor(min_wavelengh/dz));
fprintf('#################################################\n');

%% ALLOCATE MEMORY FOR WAVEFIELD
p3 = zeros(nz+2,ny+2,nx+2);            % Wavefields at t
p2 = zeros(nz+2,ny+2,nx+2);            % Wavefields at t-1
p1 = zeros(nz+2,ny+2,nx+2);            % Wavefields at t-2
% Coefficients for derivatives
co_dxx = 1/dx^2;
co_dyy = 1/dy^2;
co_dzz = 1/dz^2;

%% Loop over TIME
tic;
for it = 1:nt
    p3 = zeros(size(p2));
    % Second-order derivatives
    dp_dxx = co_dxx * (p2(2:end-1,2:end-1,1:end-2) - 2*p2(2:end-1,2:end-1,2:end-1) + p2(2:end-1,2:end-1,3:end));
    dp_dyy = co_dyy * (p2(2:end-1,1:end-2,2:end-1) - 2*p2(2:end-1,2:end-1,2:end-1) + p2(2:end-1,3:end,2:end-1));
    dp_dzz = co_dzz * (p2(1:end-2,2:end-1,2:end-1) - 2*p2(2:end-1,2:end-1,2:end-1) + p2(3:end,2:end-1,2:end-1));
    % P(t) = 2*P(t-1) - P(t-2) + G dt2;
    p3(2:end-1,2:end-1,2:end-1) = 2.0*p2(2:end-1,2:end-1,2:end-1) - p1(2:end-1,2:end-1,2:end-1) + (vp.^2).*(dp_dxx + dp_dyy + dp_dzz).*dt2;
    % Add source term
    p3(ksrc,jsrc,isrc) = p3(ksrc,jsrc,isrc) + force_x(it);
    % Exchange data between t-2 (1), t-1 (2) and t (3) and apply ABS
    p1 = p2 .* weights;
    p2 = p3 .* weights;
    % Output
    if mod(it,IT_DISPLAY) == 0
        fprintf('Time step: %d \t %.4f s\n',it, single(t(it)));
        up = permute(p3,[2 3 1]);
        figure(1); clf; subplot(2,1,1); hold on;
        h1 = slice(up, round(nx/2), round(ny/2), round(nz/2));
        set(h1,'edgecolor','none'); alpha(h1,0.6); axis equal tight; colormap jet;
        xlabel('OX'); ylabel('OY'); zlabel('OZ'); axis equal tight;
        title(['Step = ',num2str(it),'/',num2str(nt),', Time: ',sprintf('%.4f',t(it)),' sec']);
        set(gca,'ZDir','reverse'); view(52,24); grid off;
        ax_len = round(0.1 * min(nz,min(ny,nx)));
        x = line([0,ax_len],[0 0],[0,0],'color','r','linewidth',2);
        y = line([0,0],[0 ax_len],[0,0],'color','g','linewidth',2);
        z = line([0 0],[0,0],[0,ax_len],'color','b','linewidth',2);
        subplot(2,1,2); hold off;
        h2 = imagesc([squeeze(p3(round(nz/2),:,:)) squeeze(p3(:,round(ny/2),:)) squeeze(p3(:,:,round(nx/2)))]); axis equal tight;
        xlabel('XY, XZ, YZ middle slices'); colorbar; axis equal tight; drawnow;
    end
end
toc; disp('End');
