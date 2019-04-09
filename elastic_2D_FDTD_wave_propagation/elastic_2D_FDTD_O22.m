% Finite-differences in time domain(FDTD) elastic wave propagation in 2D 
% medium surrounded by simple sponge boundaries with exponential decay
% (Cerjan, 1985).
%
% We solve second order wave equation in time domain and displacement 
% formulation getting wavefield in terms of displacement vector [ux, uz].
%
% Elastic medium is parametrized by density and two velocities rho, vp
% and vs . We show CFL condition and number of points per wavelength
% prior running loop over time steps.
%
% Conventional FD star-stencils deliver accuracy O(2,2)
% Stencils: [1 -2 1]/dx2 and [1 -1 -1 1]/4dxdz
% --------------------------------------------------------------
% The code is intentionally writen in a single file
% to simplify start up.
%
% The program does not save any files, add such option manually if needed.
% Drawing the wavefield is the most computationally demanding. Increase 
% IT_DISPLAY value to reduce output and accelerate computation.
%
% The goal is to provide a simple example of elastic wave propagation
% in 2D isotropic medium.
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
IT_DISPLAY = 10;

%% MODEL
% Model dimensions, [m]
nx = 401;
nz = 401;
dx = 10;
dz = 10;

% Elastic parameters
vp = 3300.0 * ones(nz, nx);         % velocity of compressional waves, [m/s]
vs = vp / 1.732;                    % velocity of shear waves, [m/s]
rho = 2800.0 * ones(size(vp));      % density, [kg/m3]

% Lame parameters
lam = rho.*(vp.^2 - 2*vs.^2);       % first Lame parameter
mu = rho.*vs.^2;                    % shear modulus, [N/m2]

%% TIME STEPPING
t_total = .6;                       % [sec] recording duration
dt = 0.8/(max(vp(:)) * sqrt(1.0/dx^2 + 1.0/dz^2));
nt = round(t_total/dt);             % number of time steps
t = [0:nt]*dt;

CFL = max(vp(:))*dt * sqrt(1.0/dx^2 + 1.0/dz^2);
%% SOURCE
f0 = 10.0;                          % dominant frequency of the wavelet
t0 = 1.20 / f0;                     % excitation time
factor = 1e10;                      % amplitude coefficient
angle_force = 90.0;                 % spatial orientation

jsrc = round(nz/2);                 % source location along OZ
isrc = round(nx/2);                 % source location along OX

a = pi*pi*f0*f0;
dt2rho_src = dt^2/rho(jsrc, isrc);
source_term = factor * exp(-a*(t-t0).^2);                             % Gaussian
% source_term =  -factor*2.0*a*(t-t0)*exp(-a*(t-t0)^2);                % First derivative of a Gaussian:
% source_term = factor * (1.0 - 2.0*a*(t-t0).^2).*exp(-a*(t-t0).^2);        % Ricker source time function (second derivative of a Gaussian):

force_x = sin(angle_force * pi / 180) * source_term * dt2rho_src / (dx * dz);
force_z = cos(angle_force * pi / 180) * source_term * dt2rho_src / (dx * dz);
min_wavelengh = 0.5*min(vs(vs>330))/f0;     % shortest wavelength bounded by velocity in the air

%% ABSORBING BOUNDARY (ABS)
abs_thick = min(floor(0.15*nx), floor(0.15*nz));         % thicknes of the layer
abs_rate = 0.3/abs_thick;      % decay rate

lmargin = [abs_thick abs_thick];
rmargin = [abs_thick abs_thick];
weights = ones(nz+2,nx+2);
for iz = 1:nz+2
    for ix = 1:nx+2
        i = 0;
        j = 0;
        k = 0;
        if (ix < lmargin(1) + 1)
            i = lmargin(1) + 1 - ix;
        end
        if (iz < lmargin(2) + 1)
            k = lmargin(2) + 1 - iz;
        end
        if (nx - rmargin(1) < ix)
            i = ix - nx + rmargin(1);
        end
        if (nz - rmargin(2) < iz)
            k = iz - nz + rmargin(2);
        end
        if (i == 0 && j == 0 && k == 0)
            continue
        end
        rr = abs_rate * abs_rate * double(i*i + j*j + k*k );
        weights(iz,ix) = exp(-rr);
    end
end

%% SUMMARY
fprintf('#################################################\n');
fprintf('2D elastic FDTD wave propagation in isotripic \nmedium in displacement formulation with \nCerjan(1985) boundary conditions\n');
fprintf('#################################################\n');
fprintf('Model:\n\t%d x %d\tgrid nz x nx\n\t%.1e x %.1e\t[m] dz x dx\n',nz, nx, dz,dx);
fprintf('\t%.1e x %.1e\t[m] model size\n',nx*dx, nz*dz);
fprintf('\t%.1e...%.1e\t[m/s] vp\n', min(vp(:)), max(vp(:)));
fprintf('\t%.1e...%.1e\t[m/s] vs\n', min(vs(:)), max(vs(:)));
fprintf('\t%.1e...%.1e\t[kg/m3] rho\n', min(rho(:)), max(rho(:)));
fprintf('Time:\n\t%.1e\t[sec] total\n\t%.1e\tdt\n\t%d\ttime steps\n',t_total,dt,nt);
fprintf('Source:\n\t%.1e\t[Hz] dominant frequency\n\t%.1f\t[sec] index time\n',f0,t0);
fprintf('Other:\n\t%.1f\tCFL number\n', CFL);
fprintf('\t%.2f\t[m] shortest wavelength\n\t%d, %d\t points-per-wavelength OX, OZ\n', min_wavelengh, floor(min_wavelengh/dx), floor(min_wavelengh/dz));
fprintf('#################################################\n');

%% ALLOCATE MEMORY FOR WAVEFIELD
ux3 = zeros(nz+2,nx+2);            % Wavefields at t
uz3 = zeros(nz+2,nx+2);
ux2 = zeros(nz+2,nx+2);            % Wavefields at t-1
uz2 = zeros(nz+2,nx+2);
ux1 = zeros(nz+2,nx+2);            % Wavefields at t-2
uz1 = zeros(nz+2,nx+2);
% Coefficients for derivatives
co_dxx = 1/dx^2;
co_dzz = 1/dz^2;
co_dxz = 1/(4.0 * dx * dz);
co_dzx = 1/(4.0 * dx * dz);
dt2rho=(dt^2)./rho;
lam_2mu = lam + 2 * mu;

%% Loop over TIME
tic;
for it = 1:nt
    ux3 = zeros(size(ux2));
    uz3 = zeros(size(uz2));
    % Second-order derivatives
    % Ux
    dux_dxx = co_dxx * (ux2(2:end-1,1:end-2) - 2*ux2(2:end-1,2:end-1) + ux2(2:end-1,3:end));
    dux_dzz = co_dzz * (ux2(1:end-2,2:end-1) - 2*ux2(2:end-1,2:end-1) + ux2(3:end,2:end-1));
    dux_dxz = co_dxz * (ux2(1:end-2,3:end) - ux2(3:end,3:end) ...
        - ux2(1:end-2,1:end-2) + ux2(3:end,1:end-2));
    % Uz
    duz_dxx = co_dxx * (uz2(2:end-1,1:end-2) - 2*uz2(2:end-1,2:end-1) + uz2(2:end-1,3:end));
    duz_dzz = co_dzz * (uz2(1:end-2,2:end-1) - 2*uz2(2:end-1,2:end-1) + uz2(3:end,2:end-1));
    duz_dxz = co_dxz * (uz2(1:end-2,3:end) - uz2(3:end,3:end) ...
        - uz2(1:end-2,1:end-2) + uz2(3:end,1:end-2));
    % Stress G
    sigmas_ux = lam_2mu .* dux_dxx + lam .* duz_dxz + mu .* (dux_dzz + duz_dxz);
    sigmas_uz = mu .* (dux_dxz + duz_dxx) + lam .* dux_dxz + lam_2mu .* duz_dzz;
    % U(t) = 2*U(t-1) - U(t-2) + G dt2/rho;
    ux3(2:end-1,2:end-1) = 2.0*ux2(2:end-1,2:end-1) - ux1(2:end-1,2:end-1) + sigmas_ux.*dt2rho;
    uz3(2:end-1,2:end-1) = 2.0*uz2(2:end-1,2:end-1) - uz1(2:end-1,2:end-1) + sigmas_uz.*dt2rho;
    % Add source term
    ux3(jsrc, isrc) = ux3(jsrc, isrc) + force_x(it);
    uz3(jsrc, isrc) = uz3(jsrc, isrc) + force_z(it);
    % Exchange data between t-2 (1), t-1 (2) and t (3) and apply ABS
    ux1 = ux2 .* weights;
    ux2 = ux3 .* weights;
    uz1 = uz2 .* weights;
    uz2 = uz3 .* weights;
    % Output
    if mod(it,IT_DISPLAY) == 0
        fprintf('Time step: %d \t %.4f s\n',it, single(t(it)));
        u=sqrt(ux3.^2 + uz3.^2);
        imagesc(u); colorbar; colormap jet; xlabel('m'); ylabel('m');
        title(['Step = ',num2str(it),'/',num2str(nt),', Time: ',sprintf('%.4f',t(it)),' sec']);
        axis equal tight; drawnow;
    end
end
toc; disp('End');
