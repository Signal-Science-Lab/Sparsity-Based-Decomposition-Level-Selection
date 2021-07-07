% Script to generate model data
% -----------------------------
line1 = Lorenz(1e-7,1e-5,3180,3180);
line2 = Lorenz(1e-7,1e-5,3162,3180);
line3 = Lorenz(1e-7,1e-5,3198,3180);

SUPERP = [line1(:,1), line1(:,2)+ line2(:,2) + line3(:,2)];
x_y = 2*(10^8)*SUPERP(:,2);
x_y = x_y';

x = x_y; % x is model data (1 x 1024 double)
