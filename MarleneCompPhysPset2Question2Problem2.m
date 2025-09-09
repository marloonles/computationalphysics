%% MATLAB Primer Problem 2
t = [0:1000]/1000; %% 1000 element array from 0 to 1
f =  exp(5*t);
f(f >= 10) = 10;
plot(t, f)
