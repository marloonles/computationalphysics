t = [0:100]/100;
y = sin(2*pi*t);
I = (y(2:end) + y(1:end-1))/2;
I = I.*diff(t);
I = sum(I)

y2 = exp(-5*t);
I2 = (y2(2:end) + y2(1:end-1))/2 .* diff(t);
I2 = sum(I2)

y3 = cos(t).*exp(-t);
I3 = (y3(2:end) + y3(1:end-1))/2 .* diff(t);
I3 = sum(I3)


