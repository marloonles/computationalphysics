%% MATLAB Primer Problem 1
t = [0:1000]/1000 %% 1000 element array from 0 to 1
f = sin(2*pi*t)
plot(t, f, 'blue');
hold on
plot(t,exp(-5*t),'red');
plot(t,(cos(t).*exp(-t)),'yellow')
xlabel('t')
ylabel('f(t)')
legend('sin(2\pit)', 'exp(-5t)', 'cos(t) \cdot exp(-t)')
hold off


