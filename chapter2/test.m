%% TEST

x = [1+1j, -1-1j, 1-1j, -1+1j];
phi_0 = pi/6;
x_rotated = x * exp(1j*phi_0);
scatterplot(x_rotated);
