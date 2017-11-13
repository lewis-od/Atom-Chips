%% Solve laplace's equation on a square region with a circle cut out
% Requires PDE Toolbox
clear all;

% Parameters
sigma = 5.96e7; % Conductivity of conductor [S*m^-1]
V0 = 1.0; % Voltage is V0 and -V0 at ends of conductor [V]
d = 0.1; % Thickness of conductor [m]
z = 0.5; % Distance from conductor [m]

params = [sigma, V0, d, z];

% Specify geometry
R1 = [3, 4, -1, 1, 1, -1, 1, 1, -1, -1];
C1 = [1, 0, 0, 0.5]';
C1 = [C1;zeros(length(R1) - length(C1),1)];
R1 = R1';
gm = [R1, C1];
sf = 'R1-C1';
ns = char('R1', 'C1')';

[phi, xq, yq, Bx, By] = calc_field(params, gm, ns, sf, [4, 2, 1, 3, 5, 6, 7, 8]);

% Plot results
figure(1);
surf(xq, yq, phi, 'EdgeColor', 'none', 'FaceColor', 'interp');
title("Electric Potential", 'FontSize', 18);
xlabel('x', 'FontSize', 16);
ylabel('y', 'FontSize', 16);
c = colorbar();
c.Label.String = "\phi";
c.Label.FontSize = 16;
colormap jet;
axis equal;
view(2);

figure(2);
B = sqrt(Bx.^2 + By.^2);
surf(xq, yq, B, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap('jet');
c = colorbar();
c.Label.String = '|B|';
c.Label.FontSize = 16;
title('B-Field for z=0.5 and d=0.1', 'FontSize', 18);
xlabel('x', 'FontSize', 16);
ylabel('y', 'FontSize', 16);
view(2);
