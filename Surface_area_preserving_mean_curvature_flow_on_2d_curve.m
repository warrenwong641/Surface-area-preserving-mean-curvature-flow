theta = linspace(0, 2*pi, 128);

%curve 
n = 5;
epsilon = 0.1;
space_theta = theta(2) - theta(1);
r = 1 + epsilon*cos(n*theta);
X = r.*cos(theta);
Y = r.*sin(theta);

X_forward_one = [X(end), X(1:end-1)]; %shift element in X array to the right by one unit
Y_forward_one = [Y(end), Y(1:end-1)]; %shift element in Y array to the right by one unit

length_curve_initial = norm(sqrt((X_forward_one-X).^2+ (Y_forward_one-Y).^2),1) % length of initial curve
length_curve_list = []; %store the difference between length of curve at each time and length of initial curve 
t_interval = []; % record the corresponding time of each element in length_curve_list

for i=1:4000

r = sqrt(X.^2 + Y.^2); %radius at each moment

%unit tangent
T_x = gradient([X(end),X,X(1)], space_theta);
T_x = T_x(2:end-1);
T_y = gradient([Y(end),Y,Y(1)], space_theta);
T_y = T_y(2:end-1);
tangent_len = sqrt(T_x.^2 + T_y.^2);
T_x = T_x./tangent_len;
T_y = T_y./tangent_len;



%unit normal
N_x = T_y;
N_y = -T_x;

H = mean_curvature(theta, r);

%constant h
H_ds = H.*sqrt(gradient(r, theta).^2 + r.^2);
H_ds = [H_ds(end), H_ds(1:end), H_ds(1)];

H2_ds = (H.^2).*sqrt(gradient(r, theta).^2 + r.^2);
H2_ds = [H2_ds(end), H2_ds(1:end), H2_ds(1)];

h = simpson(theta(2)-theta(1), H_ds)/simpson(theta(2)-theta(1), H2_ds);

%X = (1-h*H).*N_x*0.01 + X; forward eular evolution, not used here.
%Y = (1-h*H).*N_y*0.01 + Y;

X_forward_one = [X(end), X(1:end-1)]; %shift element in X array to the right by one unit
Y_forward_one = [Y(end), Y(1:end-1)]; %shift element in Y array to the right by one unit
length_curve = norm(sqrt((X_forward_one-X).^2+ (Y_forward_one-Y).^2),1) % calculate length of curve at this time
length_curve_list = [length_curve_list, length_curve-length_curve_initial]; %append the difference between length of curve at this time and length of initial curve to the array
t_interval = [t_interval, 0.01*i]; %append the current time interval to time record array

%implementation of ode to find X and Y after dt = 0.01;
tspan = linspace(0,0.01, 500);
ic = transpose(X);
%opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t, x] = ode45(@(t,x) myode_1(t,x,h,H,N_x), tspan, ic);
X = x(end, 1:end);


ic = transpose(Y);
[t, y] = ode45(@(t,y) myode_2(t,y,h,H,N_y), tspan, ic);
Y = y(end, 1:end);

X_forward_one = [X(end), X(1:end-1)]; %shift element in X array to the right by one unit
Y_forward_one = [Y(end), Y(1:end-1)]; %shift element in Y array to the right by one unit
length_curve = norm(sqrt((X_forward_one-X).^2+ (Y_forward_one-Y).^2),1) % calculate length of curve at this time
length_curve_list = [length_curve_list, length_curve-length_curve_initial]; %append the difference between length of curve at this time and length of initial curve to the array
t_interval = [t_interval, 0.01*i]; %append the current time interval to time record array

end

plot(t_interval, length_curve_list, ' bx') % plot the graphe of the difference between length of curve at this time and length of initial curve against time
hold on;
plot(X, Y, 'rx') % plot the final curve


function dydt_x = myode_1(t,x,h,H,N_x) % function for ode45
dydt_x = transpose((1 - h*H).*N_x);
end

function dydt_y = myode_2(t,y,h,H,N_y) % function for ode45
dydt_y = transpose((1 - h*H).*N_y);
end

function [Output_arg] = simpson(h,integrand) % function for simpson rule
integral_value = 0;
i = 1;
integrand_forward = [integrand(end), integrand(1:end-1)];
integrand_backward = [integrand(2:end), integrand(1)]; 
while i <=length(integrand)
    integral_value = integral_value + h.*(integrand_forward(i) + 4.*integrand(i) + integrand_backward(i))./3;
    i = i + 1;
end
Output_arg = integral_value;
end

function H = mean_curvature(theta, r) %function for calculating mean curvature
space_theta = theta(2) - theta(1);
r_prime = gradient([r(end),r,r(1)], space_theta);
r_prime = r_prime(2:end-1);
r_double_prime = gradient([r_prime(end),r_prime,r_prime(1)], space_theta);
r_double_prime = r_double_prime(2:end-1);
H = (2*(r_prime.^2) + r.^2 - r.*r_double_prime)./((r_prime.^2 + r.^2).^1.5);
end



