clear all
close all

data = load('multi36p_cargo0.txt');

x = data(:,1);
y = data(:,2);

len = size(x,1);
nframe = 10;
msdx = zeros(nframe,1);
msdy = zeros(nframe,1);
msd = zeros(nframe,1);
for i = 1:nframe
    for j=1:len-i
        msdx(i) = msdx(i) + (x(j+i)-x(j))^2;
        msdy(i) = msdy(i) + (y(j+i)-y(j))^2;
    end
end

msd = msdx + msdy;

plot(msd)