function StreamlineTplot_2d
%STREAMLINETPLOT_2D Plot 2D results for mantle simulations.
%   Plot temperature profile of LAB for a given time series. Output the
%   cumulative surface distribution of melt migration along with a
%   comparison of the effects of heat transfer between the solid and fluid
%   phase, driven by advection.

output_interval = 1;
endind = 55;

% Root directory of data to be used
root_dir = '~/Desktop/WithNonLinearSolver/';
disp(root_dir);
mu_vals = 5e+21; % Pa s
Tbs = 1300.0;
k_s = [0.02, 0.01, 0.001, 0];

% Define constants
rho_0 = 3300.0; % SI
alpha = 2.5e-5; % thermal expansion, SI
g     = 9.81;   % SI
kappa_0 = 1.E-6;

% from MultipleRuns.py or codes like it in the folder(s) above, we
% establish the Temp scale
% temp_values = [27.+273, Tb+273, 1300.+273, 1500.+273]
% dTemp = temp_values[3] - temp_values[0]
Tscale = 1500-27;
Tval   = 1573;
h      = 1e3; % box dimension in km

hscale = 1e6; % box scale in m
tcont  = 300:100:1600; % in Kelvin
pscale = 1192135725.0; % pressure scale from MultipleRuns.py in Pa

% for streamline calculation, use the following from paraview:
% dimgradP = u*1192135725.0/1e6
% wvel = -(dimgradP-150*9.8*jHat)*1e-15/1e-2
% all SI units

rho_melt = 2800; % kg/m^3

k_over_mu = 1e-13 / 1;
stream_int = 5;
% startx  = (1:stream_int:990)'; % note here units must be in km as displayed in box
seg = 100:50:900;
startx = repmat(seg, 1, length(seg))';
nstream = length(startx) * 0.5;
% starty  = 500 * ones(length(startx), 1);
starty  = reshape(repmat(seg, length(seg), 1), length(seg)^2, 1);
startz  = 100 * ones(length(startx), 1);
prevDist = [];
numt = 1;

combTracers  = [];

% for Tb = Tbs
for vals = combvec(mu_vals, Tbs, k_s)
    cell_vals = num2cell(vals);
    [mu_scale, Tb, k] = cell_vals{:};
    
    mu_str = ['mu=' num2str(mu_scale) '/'];
    Tb_str = ['Tb=' num2str(Tb, '%3.1f') '/'];
    k_str  = ['k='  num2str(k)];
    base  = [root_dir mu_str Tb_str k_str];
    
    coords = h5read([base '/T_solid.h5'], '/Mesh/0/coordinates');
    x = h * coords(1,:); y = h * coords(2,:); z = h * coords(3,:);
    
    % Restructure the data into a matrix using the arrangement of points
    % in the arrays x, y, and z.
    
    x_stride = 0;
    y_stride = 0;
    pos = 1;
    
    while x_stride == 0 || y_stride == 0
        if x_stride == 0 && x(pos) > x(pos+1), x_stride = pos; end
        if y_stride == 0 && y(pos) > y(pos+1), y_stride = pos; end
        pos = pos + 1;
    end
    
    x_step = x_stride;
    y_step = y_stride / x_stride;
    z_step = length(x) / y_stride;
    shape = [x_step, y_step, z_step];
            
    X = reshape(x, shape);
    Y = reshape(y, shape);
    Z = reshape(z, shape);
    
    Ra      = rho_0*alpha*g*Tscale*(hscale^3)/(kappa_0*mu_scale);
    Rafac   = rho_0*alpha*g*Tscale*(hscale^3)/(kappa_0);
    vscale  = rho_0*alpha*g*Tscale*(hscale^2)/mu_scale;
    
    %% Create interpolated data sets
    % Refine a grid in each dimension
    inter = @(M) interpn(M, 2, 'spline');
    XI = inter(X); YI = inter(Y); ZI = inter(Z);
    
    % Go through all appropriate timesteps
    for ind = 1:output_interval:endind
        data_set = ['/VisualisationVector/' num2str(ind)];
        
        temperature = h5read([base '/t6t.h5'],      data_set);
        mu          = h5read([base '/mu.h5'],       data_set);
        vel         = h5read([base '/velocity.h5'], data_set);
        gradp       = h5read([base '/gradp.h5'],    data_set);
        
        % Reshape the corresponding arrays into a matrix with the
        % appropriate size. Scale to real dimensional values while doing so.
        
        T  = Tscale * reshape(temperature, shape);
        if ind == 1, T_init = T; end
        
        MU = mu_scale * reshape(mu, shape);
        
        VX = vscale * reshape(vel(1,:), shape);
        VY = vscale * reshape(vel(2,:), shape);
        VZ = vscale * reshape(vel(3,:), shape);
        
        DPDX = (pscale / hscale) * reshape(gradp(1,:), shape);
        DPDY = (pscale / hscale) * reshape(gradp(2,:), shape);
        DPDZ = (pscale / hscale) * reshape(gradp(3,:), shape);
        
        rho  = rho_0 * (1 - alpha * (T - Tval));
        drho = rho - rho_melt;
        
        WX = -k_over_mu * DPDX;
        WY = -k_over_mu * DPDY;
        WZ = -k_over_mu * (DPDZ - drho * g);
        
        Vmeltx = WX + VX;
        Vmelty = WY + VY;
        Vmeltz = WZ + VZ;
        
        reg = @(x,n) linspace(min(x), max(x), n);
        [XG, YG, ZG] = meshgrid(reg(x, x_step), reg(y, y_step), reg(z, z_step));
        
        %% Output
        figure(1); clf('reset'); hold on;
        
        [faces,LAB,colors] = isosurface(XI, YI, ZI, inter(T), Tval, inter(DPDX));
        patch('Vertices', LAB, ...
              'Faces', faces, ...
              'FaceVertexCData', colors, ...
              'FaceColor','interp', ...
              'edgecolor', 'none');
        
        colormap(jet(10000))
        colorbar
        shading interp
        material metal
        camlight('headlight','infinite');
        lighting gouraud % look into specularity
        daspect([1,1,1])
        view(3)
        
        figure(2); clf('reset'); hold on;
        scatter3(LAB(:,1), LAB(:,2), LAB(:,3));
        clith = LAB(:,3);
        zll   = mean(clith);
        surf(X(:,:,1), Y(:,:,1), zll * ones(x_step, y_step));
        view(3);
        
        % make an array containting the mean depth of the Tval contour as a fn
        % of time
        % find mean viscosity in the convecting interior  - used to find Ra_i
        mu_int = mean(MU(Z < zll));
        Ra_int = Rafac / mu_int;
        
        numt = numt + 1;
        
        figure(1); hold on; grid on; % overlay streamlines and velocity vectors
        scatter3(startx, starty, startz, 'o');
        
        axis([min(x), max(x), min(y), max(y), min(z), max(z)])
        set(gca, 'BoxStyle', 'full', 'Box', 'on')
        %         [sx, sy, sz] = meshgrid(startx, starty, startz);
        han = streamline(XG, YG, ZG, Vmeltx, Vmelty, Vmeltz, startx, starty, startz);
        
        %         set(han, 'color', 'r', 'linewidth', 1.25);
        quiver3(X, Y, Z, VX, VY, VZ, 1, 'k');
        
        trackStream(han, nstream, min(x), max(x));
        
        % Find pressure-gradients along average isotherm depth, zll
        Z_diff = abs(Z - zll);
        close = Z_diff == min(Z_diff(:));
        xpos  = X(close);
        p_x   = DPDX(close); % profile of dpdx along zll
        figure(4); clf('reset');
        % plot(xpos, p_x, 'k-');
        surf(X(:,:,1), Y(:,:,1), reshape(p_x, [x_step, y_step]))
        % Color LAB surface by dpdx
        % Calculate the force on the protrusion due to pressure difference
        % could find the deformation expected from force
        
        %         set(gca,'fontname','Helvetica','fontsize',[14],'ylim',[-150 150])
        %         xlabel('km'); ylabel('Pa/m')
        %         %output to file for each timestep
        %         filename = ['dpdx_zll_t_' num2str(ind)];
        %         dat = [xpos p_x];
        %         WD1 = cd;
        %         cd(base)
        %         eval(['save ' filename ' dat -ascii'])
        %         cd(WD1)
        
        T_0    = 1300+273;
        crust_thickness = 30;
        
        Z_Crust = max(Z(:)) - crust_thickness;
        Z_comp = min(clith(:)); % The greatest depth with no large lateral density change
        Z_Mantle = (Z < Z_Crust) & (Z >= Z_comp);
        dz = Z(1,1,2) - Z(1,1,1);
        dy = Y(1,2,1) - Y(1,1,1);
        dx = X(2,1,1) - X(1,1,1);
        
        delT     = T - T_init;
        rhoarr   = rho_0 * (1 - alpha*(T - T_0));
        integral = cumsum(delT .* Z_Mantle, 3);
        ru_iso   = alpha * dz * 1e3 * integral(:, :, shape(3));
        ru_iso   = ru_iso - ru_iso(1); % Relative to edge
        
        surf(XI(:,:,1), YI(:,:,1), inter(ru_iso));
        
        [dVX_x, dVX_y, dVX_z]  = gradient(VX,dx,dy,dz);
        [dVY_x, dVY_y, dVY_z]  = gradient(VY,dx,dy,dz);
        [dVZ_x, dVZ_y, dVZ_z]  = gradient(VZ,dx,dy,dz);
        
        Sxx = 2 * MU .* dVX_x;
        Syy = 2 * MU .* dVY_y;
        Szz = 2 * MU .* dVZ_z;
        Sxz = MU .* (dVZ_x + dVX_z);
        Syz = MU .* (dVZ_y + dVY_z);
        
        [dSxz_x, dSxz_y, dSxz_z] = gradient(Sxz,dx,dy,dz);
        [dSyz_x, dSyz_y, dSyz_z] = gradient(Syz,dx,dy,dz);
        [dSzz_x, dSzz_y, dSzz_z] = gradient(Szz,dx,dy,dz);
        
        fac = rho_0 * g;
        
        PZ_gravity = cumsum(-fac * (ones(size(Z))*max(Z(:)) - Z), 3);
        PZ_x = cumsum(dSxz_x * dz, 3);
        PZ_y = cumsum(dSyz_y * dz, 3);
        PZ_z = cumsum(dSzz_z * dz, 3);
        
        figure(10);
        lith_interp = scatteredInterpolant(LAB(:,1:2), LAB(:,3), 'linear');
        
        clith_vals = lith_interp(X(:,:,1), Y(:,:,1));
        zs = round(clith_vals / dz) + 1;
        
        for ii = 1:shape(1)
            for jj = 1:shape(2)
                z_indices(ii,jj) = ii + jj * x_step + zs(ii,jj)*x_step*y_step;
            end
        end
        
        PZ = PZ_gravity + PZ_x + PZ_y + PZ_z;
        combined = PZ(z_indices) / fac + ru_iso;
        surf(X(:,:,1), Y(:,:,1), combined);
        %         surf(X(:,:,1), Y(:,:,1), PZ_z(z_indices) / fac);
        %         figure(11);
        %         surf(X(:,:,1), Y(:,:,1), PZ_y(z_indices) / fac);
        %         figure(12);
        %         surf(X(:,:,1), Y(:,:,1), PZ_x(z_indices) / fac);
        %         figure(13);
        %         surf(X(:,:,1), Y(:,:,1), PZ_gravity(z_indices) / fac);
        %         Plith    = g*cumsum(rhoarr,1)*dy*1e3;
        
        drawnow
        %         pause(0.25);
        %         input('continue')
    end
end
end

function newTracers = trackStream(newStreams, n, lower, upper)
numStreams = length(newStreams);
xEndPoints = zeros(1,numStreams);
yEndPoints = zeros(1,numStreams);

for index = 1:numStreams
    xData = get(newStreams(index), 'XData');
    yData = get(newStreams(index), 'YData');
    xEndPoints(index) = xData(end);
    yEndPoints(index) = yData(end);
end

figure(5); clf; hold on;
subplot(2,1,1); hold on; title('Bin Dist at Current step');
% hist(xEndPoints, n); xlim([lower,upper]);
hist3([xEndPoints', yEndPoints'], [20,20]);
view(3);

newTracers = histc(xEndPoints, linspace(lower,upper,n));
% subplot(2,1,2); hold on; title('Cumulative distribution');
% area(combTracers'); xlim([1,n]);

%         if ind == 51
%             figure(6);
%             colormap(hot);shading faceted
%
%             subplot(3,2,1);  %title('Cumulative distribution');
%             xrange = linspace(lower,upper,n);
%             plotx  = xrange(3:end-3);
%             indsout = [1, 3, 5, 7, 9, 11];
%             ploty  = combTracers(indsout,3:end-3);
%             area(plotx,ploty'); xlim([lower,upper]);
%             set(gca,'fontname','Helvetica','fontsize',[14]);
%             set(gca,'xlim',[0 1000],'yTickLabel',' ','box','off')
%
%             clear ploty;
%             ploty  = combTracers(:,3:end-3);
%             subplot(3,2,2);
%             m1 = mean(ploty(1:2,:));
%             m2 = mean(ploty(3:6,:));
%             m3 = mean(ploty(7:end,:));
%             %normalize to 1
%             m1 = m1/max(m3);m2 = m2/max(m3); m3 = m3/max(m3);
%             plot(plotx, movingmean(m1',15),'color',[0 0 0],'linewidth',[2]); hold on
%             plot(plotx, movingmean(m2',15),'color',[1 0 0],'linewidth',[2]);
%             plot(plotx, movingmean(m3',15),'color',[0 0.2 0.8],'linewidth',[2]);
%             xlim([lower,upper]);
%             legend('0-10 my', '10-30 my','30-50 my', 'location','EastOutside')
%             set(gca,'fontname','Helvetica','fontsize',[14]);
%             set(gca,'xlim',[0 1000],'ylim',[0 1])
%
%         end
end