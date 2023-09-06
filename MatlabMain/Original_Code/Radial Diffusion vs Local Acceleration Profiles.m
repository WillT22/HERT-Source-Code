close all

subplot(1,2,1)
x = 1:1000;
y = zeros(1,length(x));
for i = 1:length(x)
    if i >50
        y(i) = 5*log((x(i)-50).^2);
    end
end
plot(x,y,'LineWidth',2)
set(gca,'XTick',[], 'YTick', [])
ylabel('Electron Density')
xlabel('Radius')
title('Inward Radial Diffusion Radial Profile')

subplot(1,2,2)
for i = 1:length(x)
    if i >50
        y(i) = 5*log((x(i)-50).^2).^2/x(i);
    end
end

plot(x,y,'LineWidth',2)
set(gca,'XTick',[], 'YTick', [])
ylabel('Electron Density')
xlabel('Radius')
title('Local Acceleration Radial Profile')