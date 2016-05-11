clear all
close all

N = 20;

% first gen first line
% attraction is 2.158, dl is 2.2
a = 2.14/2.0;
tag = 'd=2.14'
size = N;
    line_y = zeros(size,1);
    line_x = 2*a*[1:size]';
    config = [line_x line_y];
    for i = 1:size
        line_x = 2*a*[1:size]' - a;
        line_y = line_y + sqrt(3)*a;
        if(mod(i,2) == 0)
            line_x = 2*a*[1:size]';
        end
        config = [config; line_x line_y];
    end
    
    config(:,1) = config(:,1) - mean(config(:,1));
    config(:,2) = config(:,2) - mean(config(:,2));
    rx = config(:,1);
    ry = config(:,2);
rz = zeros(length(rx),1);

dist = rx.^2 + ry.^2;
[B, index] = sort(dist);

phi = randn(length(rx),1);
theta = zeros(length(rx),1);
config = [];


config = [[1:length(rx)]' rx(index) ry(index) rz phi theta];



config(:,2)=config(:,2)-config(1,2);
config(:,3)=config(:,3)-config(1,3);
figure(1)
plot(config(:,2),config(:,3),'linestyle','none','marker','o');

np = length(config(:,2));

for i = 1:np
    figure(1)
    x1 = config(i,2);
    y1 = config(i,3);
    txt1 = num2str(i);
    text(x1,y1,txt1);
end

targetSubset = [3 4 7];
configSubset = config(targetSubset,:);
dlmwrite(['target_3p' tag '.txt'],configSubset,'delimiter','\t');
figure(3)
plot(configSubset(:,2),configSubset(:,3),'linestyle','none','marker','o');
xlim([-10 10])
ylim([-10 10])
axis off

targetSubset = [2:7];
configSubset = config(targetSubset,:);
dlmwrite(['target_6p' tag '.txt'],configSubset,'delimiter','\t');
figure(4)
plot(configSubset(:,2),configSubset(:,3),'linestyle','none','marker','o');
xlim([-10 10])
ylim([-10 10])
axis off


targetSubset = [2:19];
configSubset = config(targetSubset,:);
dlmwrite(['target_18p' tag '.txt'],configSubset,'delimiter','\t');
figure(5)
plot(configSubset(:,2),configSubset(:,3),'linestyle','none','marker','o');
xlim([-10 10])
ylim([-10 10])
axis off

targetSubset = [2:36 39];
configSubset = config(targetSubset,:);
dlmwrite(['target_36p' tag '.txt'],configSubset,'delimiter','\t');
figure(6)
plot(configSubset(:,2),configSubset(:,3),'linestyle','none','marker','o');
xlim([-10 10])
ylim([-10 10])
axis off


targetSubset = [2:65];
targetSubset = setdiff(targetSubset,[64 59 60 63]);
configSubset = config(targetSubset,:);
dlmwrite(['target_60p' tag '.txt'],configSubset,'delimiter','\t');
figure(7)
plot(configSubset(:,2),configSubset(:,3),'linestyle','none','marker','o');
xlim([-10 10])
ylim([-10 10])
axis off


targetSubset = [2:101];
targetSubset = setdiff(targetSubset,[98 100 93 89 92 91 90 94 97 99]);
configSubset = config(targetSubset,:);
dlmwrite(['target_90p' tag '.txt'],configSubset,'delimiter','\t');
figure(8)
plot(configSubset(:,2),configSubset(:,3),'linestyle','none','marker','o');
xlim([-15 15])
ylim([-15 15])
axis off




configSubset = config([100:290],:);
dlmwrite('config_forlocking.txt',configSubset,'delimiter','\t');