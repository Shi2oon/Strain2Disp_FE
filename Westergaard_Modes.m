function [files,Operation,unit,E,nu] = ...
    Westergaard_Modes(StressIntensityFactor,state,Mode,stPs,file2)
   
if ~exist('Mode','var');        Mode = 'I';         end
if isempty(Mode);               Mode = 'I';         end

% close all; clear; clc
if exist('file2','var')
    newfile = fullfile(file2, [num2str(StressIntensityFactor) '_' Mode ' WS']);
else
    newfile = [ pwd '\2D_' num2str(StressIntensityFactor) '_' Mode ' WS']; 
end
mkdir(newfile);

minGrid = -4E-3; % [m]
if exist('stPs','var')
    gridStep = 8e-3/stPs;
else
    gridStep = 0.2E-3; % [m]
end

maxGrid = 4E-3; % [m]

xvec = minGrid : gridStep : maxGrid; % [m]
yvec = minGrid : gridStep : maxGrid; % [m]

[x,y] = meshgrid(xvec,yvec); % [m]
% StressIntensityFactor = 30; % [MPa m^0.5]
fprintf('preparing synthetic Westergaard Solution Data .. ');
K = StressIntensityFactor * 1E6; % Stress intensity factor [Pa m^0.5]
E = 210E9; % Young's Modulus [Pa]
nu = 0.3; % poisson ratio
mu = E/(2.*(1+nu)); % Shear Modulus [Pa]

    switch state
        case 'plane_strain'
            kappa = 3 - (4 .* nu); % [/]
        case 'plane_stress'
            kappa = (3 - nu)./(1 + nu); % [/]
    end
[theta,r] = cart2pol(x,y);
switch Mode
    case 'I'
        ux = (K./(2.*mu)).*sqrt(r/(2*pi)).*cos(theta/2).*(kappa - 1 + ...
              2.*(sin(theta/2)).^2); % Anderson p99 A2.44a
        uy = (K./(2.*mu)).*sqrt(r/(2*pi)).*sin(theta/2).*(kappa + 1 - ...
              2.*(cos(theta/2)).^2); % Anderson p99 A2.44b
        S11 = (K./sqrt(r.*2*pi)).*cos(theta/2).*(1 - (sin(theta/2).*sin(3*theta/2))); 
        S22 = (K./sqrt(r.*2*pi)).*cos(theta/2).*(1 + (sin(theta/2).*sin(3*theta/2))); 
        S12 = (K./sqrt(r.*2*pi)).*cos(theta/2).*(sin(theta/2).*cos(3*theta/2));

    case 'II'
        ux = (K./(2.*mu)).*sqrt(r/(2*pi)).*sin(theta/2).*(kappa + 1 + ...
              2.*(cos(theta/2)).^2); % Anderson p99 A2.44a
        uy = (K./(2.*mu)).*sqrt(r/(2*pi)).*cos(theta/2).*(kappa - 1 - ...
              2.*(sin(theta/2)).^2); % Anderson p99 A2.44b
        S11 = (-K./sqrt(r.*2*pi)).*sin(theta/2).*(2 + (cos(theta/2).*cos(3*theta/2))); 
        S22 = (K./sqrt(r.*2*pi)).*sin(theta/2).*cos(theta/2).*cos(3*theta/2); 
        S12 = (K./sqrt(r.*2*pi)).*cos(theta/2).*(1-(sin(theta/2).*sin(3*theta/2)));
    case 'fun'
        ux = (rand(1)*(K./(2.*mu)).*sqrt(r/(2*pi)).*sin(theta/2).*(kappa + 1 + ...
              2.*(cos(theta/2)).^2)+rand(1)*(K./(2.*mu)).*sqrt(r/(2*pi)).*cos(theta/2).*(kappa - 1 + ...
              2.*(sin(theta/2)).^2)); % Anderson p99 A2.44a
        uy = (rand(1)*(K./(2.*mu)).*sqrt(r/(2*pi)).*cos(theta/2).*(kappa - 1 - ...
              2.*(sin(theta/2)).^2)+rand(1)*(K./(2.*mu)).*sqrt(r/(2*pi)).*sin(theta/2).*(kappa + 1 - ...
              2.*(cos(theta/2)).^2)); % Anderson p99 A2.44b
end
subplot(1,3,1); contourf(x,y,ux); 
axis image; box off; c=colorbar;    c.Label.String= ['U_{x} [m]'];
subplot(1,3,2);contourf(x,y,uy); 
axis image; box off; c=colorbar;    c.Label.String= ['U_{y} [m]'];
subplot(1,3,3);contourf(x,y,sqrt(uy.^2+ux.^2)); 
axis image; box off; c=colorbar;    c.Label.String= ['U_{mag} [m]'];
set(gcf,'Position',[400 84 1209 1026]);
saveas(gcf, [newfile '\' Mode '_Disp_fields.tiff']);   
saveas(gcf, [newfile '\' Mode '_Disp_fields.fig']); close  
[dux_dx,dux_dy] = gradient(ux,gridStep);
[duy_dx,duy_dy] = gradient(uy,gridStep);
exx = dux_dx;
eyy = duy_dy;
exy = 0.5*(dux_dy + duy_dx);

    %% save displacement fields
    Operation{1,1} = 'DIC';               unit{1,1}  = 'm';
    alldata = [x(:) y(:) ux(:) uy(:)]; % [m]
    if exist('stPs','var')
        files{1,1} = [newfile '\' num2str(StressIntensityFactor) 'MPa_m_DISP_' num2str(stPs) '.mat'];
    else
        files{1,1} = [newfile '\' num2str(StressIntensityFactor) 'MPa_m_DISP.mat'];
    end
    save(files{1,1},'alldata')
    
    Operation{1,2} = 'DIC';               unit{1,2}  = 'mm';
    alldata = [x(:) y(:) ux(:) uy(:)].*1e3; % [m]
    if exist('stPs','var')
        files{1,2} = [newfile '\' num2str(StressIntensityFactor) 'MPa_mm_DISP_' num2str(stPs) '.mat'];
    else
        files{1,2} = [newfile '\' num2str(StressIntensityFactor) 'MPa_mm_DISP.mat'];
    end
    save(files{1,2},'alldata')
    
    Operation{1,3} = 'DIC';               unit{1,3}  = 'um';
    alldata = [x(:) y(:) ux(:) uy(:)].*1e6; % [m]
    if exist('stPs','var')
        files{1,3} = [newfile '\' num2str(StressIntensityFactor) 'MPa_um_DISP_' num2str(stPs) '.mat'];
    else
        files{1,3} = [newfile '\' num2str(StressIntensityFactor) 'MPa_um_DISP.mat'];
    end
    save(files{1,3},'alldata')

    %% save strain fields
    Operation{2,1} = 'Str';               unit{2,1}  = 'm';
    alldata = [x(:) y(:) exx(:) eyy(:) exy(:)]; % [m]
    if exist('stPs','var')
        files{2,1} = [newfile '\' num2str(StressIntensityFactor) 'MPa_m_Strain_' num2str(stPs) '.mat'];
    else
        files{2,1} = [newfile '\' num2str(StressIntensityFactor) 'MPa_m_Strain.mat'];
    end
    save(files{2,1},'alldata')
    
    Operation{2,2} = 'Str';               unit{2,2}  = 'mm';
    alldata = [x(:).*1e3 y(:).*1e3 exx(:) eyy(:) exy(:)]; % [m]
    if exist('stPs','var')
        files{2,2} = [newfile '\' num2str(StressIntensityFactor) 'MPa_mm_Strain_' num2str(stPs) '.mat'];
    else
        files{2,2} = [newfile '\' num2str(StressIntensityFactor) 'MPa_mm_Strain.mat'];
    end
    save(files{2,2},'alldata')
    
    Operation{2,3} = 'Str';               unit{2,3}  = 'um';
    alldata = [x(:).*1e6 y(:).*1e6 exx(:) eyy(:) exy(:)]; % [m]
    if exist('stPs','var')
        files{2,3} = [newfile '\' num2str(StressIntensityFactor) 'MPa_um_Strain_' num2str(stPs) '.mat'];
    else
        files{2,3} = [newfile '\' num2str(StressIntensityFactor) 'MPa_um_Strain.mat'];
    end
    save(files{2,3},'alldata')

fprintf ('DONE\n\n');