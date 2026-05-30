%% 3-D Peptide Simulation with Helix/Coil Propensity and Residue Labels
clear; clc; close all;

% -------------------------
% PARAMETERS
% -------------------------
N            = 31; % gets changed based on sequence length
mass         = 1;
dt           = 0.002;
steps        = 3000;
restLen      = 1;        % backbone bond rest length
base_k       = 5;        % backbone spring stiffness
nonbond_k    = 1;        % non-bonded spring constant
damping      = 0.2;
cutoff       = 4;        % cutoff for non-bonded interactions
maxSpringLen = 5;        % maximum spring length

% -------------------------
% AMINO ACID SEQUENCE
% -------------------------
seq     = 'MGINTRELFLNFTIVLITVILMWLLVRSYQY';
AAorder = 'ACDEFGHIKLMNPQRSTVWY';
aaIndex = arrayfun(@(c) find(AAorder == c), seq);

% -------------------------
% HELIX AND COIL PROPENSITIES  (Chou-Fasman scale)
% -------------------------
helixProp = [1.45 0.77 1.01 1.51 1.13 0.53 1.00 1.08 1.16 1.34 ...
             1.20 0.73 0.59 1.17 0.79 0.79 0.82 1.06 1.08 0.69];

coilProp  = [0.66 0.86 1.01 0.26 0.35 1.56 0.87 0.47 1.01 0.58 ...
             0.60 1.56 1.52 0.99 1.14 1.43 0.96 0.50 0.41 0.97];

% -------------------------
% AFFINITY MATRIX (original)
% -------------------------
Affinity = [
   0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25;
   0.25 1.00 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25;
   0.25 0.25 -0.25 -0.25 0.25 0.25 0.25 0.25 1.00 0.25 0.25 0.25 0.25 0.25 1.00 0.25 0.25 0.25 0.25 0.25;
   0.25 0.25 -0.25 -0.25 0.25 0.25 0.25 0.25 1.00 0.25 0.25 0.25 0.25 0.25 1.00 0.25 0.25 0.25 0.25 0.25;
   0.25 0.25 0.25 0.25 0.50 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.50 0.50;
   0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25;
   0.25 0.25 0.25 0.25 0.50 0.25 0.25 0.25 0.75 0.25 0.25 0.75 0.25 0.75 0.75 0.25 0.25 0.25 0.50 0.50;
   0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25;
   0.25 0.25 1.00 1.00 0.75 0.25 0.25 0.25 0.10 0.25 0.25 0.75 0.25 0.75 0.10 0.25 0.25 0.25 0.75 0.75;
   0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25;
   0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25;
   0.25 0.25 0.25 0.25 0.25 0.25 0.75 0.25 0.75 0.25 0.25 0.25 0.75 0.75 0.75 0.75 0.75 0.25 0.25 0.75;
   0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25;
   0.25 0.25 0.25 0.25 0.25 0.25 0.75 0.25 0.75 0.25 0.25 0.25 0.25 0.75 0.75 0.75 0.75 0.25 0.25 0.75;
   0.25 0.25 1.00 1.00 0.75 0.25 0.25 0.25 0.10 0.25 0.25 0.75 0.25 0.75 -0.25 0.75 0.75 0.25 0.75 0.75;
   0.25 0.25 0.25 0.25 0.25 0.25 0.75 0.25 0.75 0.25 0.25 0.75 0.25 0.75 0.75 0.75 0.75 0.25 0.25 0.75;
   0.25 0.25 0.25 0.25 0.25 0.25 0.75 0.25 0.75 0.25 0.25 0.75 0.25 0.75 0.75 0.75 0.75 0.25 0.25 0.75;
   0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25;
   0.25 0.25 0.25 0.25 0.50 0.25 0.50 0.25 0.75 0.25 0.25 0.25 0.25 0.25 0.75 0.25 0.25 0.25 0.50 0.50;
   0.25 0.25 0.25 0.25 0.50 0.25 0.50 0.25 0.75 0.25 0.25 0.75 0.25 0.75 0.75 0.75 0.75 0.25 0.50 0.50];

% -------------------------
% INITIAL POSITIONS AND VELOCITIES
% -------------------------
pos = [(0:N-1)' * restLen, zeros(N,1), zeros(N,1)] + 0.05*randn(N,3);
vel = zeros(N,3);

% -------------------------
% COLOR MAP BY RESIDUE PHYSICOCHEMISTRY
%   Hydrophobic : VILMFYW  → blue
%   Polar       : STNQ     → green
%   Charged(+)  : KRH      → red
%   Charged(-)  : DE       → orange
%   Special     : CGP      → purple
% -------------------------
colorMap = zeros(N,3);
for i = 1:N
    aa = seq(i);
    if any(aa == 'VILMFYW'),     colorMap(i,:) = [0.1  0.4  0.8 ];
    elseif any(aa == 'STNQ'),    colorMap(i,:) = [0.1  0.7  0.4 ];
    elseif any(aa == 'KRH'),     colorMap(i,:) = [0.85 0.2  0.2 ];
    elseif any(aa == 'DE'),      colorMap(i,:) = [0.9  0.5  0.1 ];
    else,                        colorMap(i,:) = [0.55 0.2  0.75];
    end
end

% -------------------------
% FIGURE SETUP
% -------------------------
fig  = figure('Name','3-D Peptide Simulation','Color','w');
ax   = axes('Parent', fig);

hLine = plot3(ax, pos(:,1), pos(:,2), pos(:,3), '-', ...
              'Color',[0.6 0.6 0.6], 'LineWidth', 1.5);
hold(ax, 'on');

hScat = scatter3(ax, pos(:,1), pos(:,2), pos(:,3), 60, colorMap, 'filled');

% --- Residue labels: number on bead, AA letter just above ---
hText = gobjects(N,1);
for i = 1:N
    hText(i) = text(ax, pos(i,1), pos(i,2), pos(i,3), ...
                    sprintf('%d\n%s', i, seq(i)), ...
                    'FontSize', 6, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment',   'bottom', ...
                    'Color', 'k');
end

grid(ax,'on'); axis(ax,'equal');
xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
view(ax, 45, 20); rotate3d(ax, 'on');

% Legend
patch(ax,NaN,NaN,NaN,'FaceColor',[0.1  0.4  0.8 ],'DisplayName','Hydrophobic (VILMFYW)');
patch(ax,NaN,NaN,NaN,'FaceColor',[0.1  0.7  0.4 ],'DisplayName','Polar (STNQ)');
patch(ax,NaN,NaN,NaN,'FaceColor',[0.85 0.2  0.2 ],'DisplayName','Charged+ (KRH)');
patch(ax,NaN,NaN,NaN,'FaceColor',[0.9  0.5  0.1 ],'DisplayName','Charged- (DE)');
patch(ax,NaN,NaN,NaN,'FaceColor',[0.55 0.2  0.75],'DisplayName','Special (CGP)');
legend(ax,'show','Location','bestoutside','FontSize',8);

hTitle = title(ax, 'Step 0');

% -------------------------
% SIMULATION LOOP
% -------------------------
for t = 1:steps

    F = zeros(N,3);

    % --- BACKBONE HARMONIC SPRINGS ---
    for i = 1:N-1
        r_vec = pos(i+1,:) - pos(i,:);
        dist  = norm(r_vec);
        dist  = min(dist, maxSpringLen);
        dir   = r_vec / (dist + 1e-12);
        F_back = base_k * (dist - restLen) * dir;
        F(i,:)   = F(i,:)   + F_back;
        F(i+1,:) = F(i+1,:) - F_back;
    end

    % --- NON-BONDED INTERACTIONS (harmonic, affinity-modulated) ---
    % Replace this block with your LJ potential when ready.
    for i = 1:N-2
        for j = i+2:N
            r_vec = pos(j,:) - pos(i,:);
            dist  = norm(r_vec);
            if dist < cutoff && dist <= maxSpringLen
                dir  = r_vec / (dist + 1e-12);
                k_ij = nonbond_k * Affinity(aaIndex(i), aaIndex(j)) * ...
                       sqrt(helixProp(aaIndex(i)) * helixProp(aaIndex(j))) / ...
                       sqrt(coilProp(aaIndex(i))  * coilProp(aaIndex(j)));
                rest_nonbond = 1.5;
                F_non = k_ij * (dist - rest_nonbond) * dir;
                F(i,:) = F(i,:) + F_non;
                F(j,:) = F(j,:) - F_non;
            end
        end
    end

    % --- DAMPING ---
    F = F - damping * vel;

    % --- INTEGRATE ---
    vel = vel + (F / mass) * dt;
    pos = pos + vel * dt;

    % --- UPDATE PLOT every 5 steps ---
    if mod(t, 5) == 0
        set(hLine, 'XData', pos(:,1), 'YData', pos(:,2), 'ZData', pos(:,3));
        set(hScat, 'XData', pos(:,1), 'YData', pos(:,2), 'ZData', pos(:,3));
        for i = 1:N
            set(hText(i), 'Position', pos(i,:));
        end
        set(hTitle, 'String', sprintf('Step %d / %d', t, steps));
        drawnow limitrate;
    end

end

fprintf('Simulation complete.\n');
