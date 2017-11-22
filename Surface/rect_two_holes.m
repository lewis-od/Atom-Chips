%% Solve laplace's equation on a rectanguar region with 2 circles cut out
% Requires PDE Toolbox
clear all;

%% Parameters
sigma = 5.96e7; % Conductivity of conductor [S*m^-1]
V0 = 1.0; % Voltage is V0 and -V0 at ends of conductor [V]
d = 0.1; % Thickness of conductor [m]
z = 0.5; % Distance from conductor [m]

%% Specify geometry
R1 = [3, 4, -1, 1, 1, -1, 2, 2, -2, -2]; % 1x2 Rectangle
C1 = [1, 0, 0.75, 0.5]'; % Circle of radius 0.5 centreed at (0, 1)
C1 = [C1;zeros(length(R1) - length(C1),1)];
C2 = [1, 0, -0.75, 0.5]'; % Circle of radius 0.5 centreed at (0, -1)
C2= [C2;zeros(length(R1) - length(C2),1)];
R1 = R1';
gm = [R1, C1, C2];
sf = 'R1-C1-C2';
ns = char('R1', 'C1', 'C2')';

%% Calculate the electric potential
[phi, xq, yq, model, result] = calc_potential(V0, gm, ns, sf);

%% Calculate B for varying z and plot how it varies in the x=0 and y=0 planes
dims = size(xq);
% Indexes of 0 entry in xq and qy arrays
x0 = ceil(dims(2)/2);
y0 = ceil(dims(1)/2);

zq = linspace(0, 1.5, 7); % z values to evaluate b at
zq = zq(2:end); % Don't evaluate at z=0
min_B = zeros(1, length(zq)); % Array to hold min value of B for each z
figure();
hold on;
for i = 1:length(zq)
    z = zq(i);
    [Bx, By] = calc_field(phi, sigma, d, z);
    B = sqrt(Bx.^2 + By.^2);
    min_B(i) = B(y0, x0); % B_min is at the origin in this case
    yq_slice = yq(:, x0); % y values corresponding to x=0
    xq_slice = xq(y0, :); % x values corresponding to y=0
    B_slice_y = B(:, x0); % B values corresponding to x=0
    B_slice_x = B(y0, :); % B values corresponding to y=0
    subplot(1, 2, 1);
    hold on;
    plot(yq_slice, B_slice_y);
    subplot(1, 2, 2);
    hold on;
    plot(xq_slice, B_slice_x);
end
subplot(1, 2, 1);
title('x=0');
xlabel('y');
ylabel('|B| [T]');
legend(cellstr(num2str(zq', 'z=%.2f')));
hold off;

subplot(1, 2, 2);
title('y=0');
xlabel('x');
ylabel('|B| [T]');
legend(cellstr(num2str(zq', 'z=%.2f')));
hold off;

figure();
plot(zq, min_B);
xlabel('z');
ylabel('B(0, 0, z)');

% figure();
% for i = 1:4
%     subplot(2, 2, i);
%     z = 0.1*i;
%     [Bx, By] = calc_field(phi, sigma, d, z);
%     B = sqrt(Bx.^2 + By.^2);
%     xq_slice = xq(y0, :);
%     B_slice = B(y0, :);
%     plot(xq_slice, B_slice);
%     title(['z=' num2str(z) '; y=0']);
%     xlabel('x');
%     ylabel('B');
% end
