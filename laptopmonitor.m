monitor.viewDist = 58; %cm
monitor.center = screenRect(3:4)./2;
monitor.size(1) = 29; %Also cm
monitor.size(2)= monitor.size(1).*(monitor.center(2)./monitor.center(1));

save('laptopparameter.mat')