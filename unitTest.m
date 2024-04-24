
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% InternalGravityWaveSpectrumUnitTest
% 
% Leticia fabre de Lima
%
% April, 2024   Version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is Coriolis frequency correct?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if im.f0 >= -1 && im.f0<=1
    disp("Coriolis frequency has a resonable value")
else
    disp("Error: Coriolis frequency DOES NOT have a resonable value")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to check if N2 is right or "smooth enough"?
% Is N2 smooth enouth?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IDEIA 1
N2Diff = diff(im.N2zInitial);
threshold=10^-10;

if N2Diff <= threshold
    disp("N2 is smooth")
else
    disp("N2 IS NOT smooth enough")
end

% IDEA 2:

N2Compare = smoothdata(im.N2zInitial);
threshold=10^-10;

if abs(N2Diff - im.N2zInitial) <= threshold
    disp("N2 is smooth")
else
    disp("N2 IS NOT smooth enough")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visual check of modes on highest freq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Check on figure 1 if most part of variance is constrained between turning points")

figure(1)
subplot(1,3,1)
sgtitle("Check if most part of variance is btw turning points")

plot(sqrt(im.N2zInitial)*3600/(2*pi),im.zInitial ,'k',LineWidth=1.5)
hold on
xline(0.8*sqrt(im.N2max)*3600/(2*pi))
ylabel('depth')
xlabel('cph')
title('N(z)')
ylim([-im.Lz*0.25 0])

subplot(1,3,2)
plot(im.FInitial(:,end),im.zInitial ,'k',LineWidth=1.5)                      
title('Initial F - Highest Mode')   
ylim([-im.Lz*0.25 0])

subplot(1,3,3)
plot(im.GInitial(:,end),im.zInitial ,'k',LineWidth=1.5)                        
title('Initial G - Highest Mode')   
ylim([-im.Lz*0.25 0])       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is Kradial equally spaced?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diffLin = diff(im.KRadialLin);
if sum(abs(diff(diffLin))) <= 10^-13
    disp("KRadial Linear is equally spaced in linear scale")
else
    disp("Error: KRadial Linear IS NOT equally spaced in linear scale")
end

diffLog = diff(log(im.KRadialLog));
if sum(abs(diff(diffLog))) <= 10^-13  
    disp("KRadial Log is equally spaced in log scale")
else
    disp("Error: KRadial Log IS NOT equally spaced in log scale")
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visual check of Kradial distribution in linear scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(1:length(im.KRadialLin),(2*pi./im.KRadialLin),".")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is the size of main matrix correct after scattered interp?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(im.F) == [length(im.zNew) im.nModes length(im.KRadialLin)]
    disp("Size of F matriz is right")
else
    disp("Error: Size of F matriz IS NOT right")
end

if size(im.G) == [length(im.zNew) im.nModes length(im.KRadialLin)]
    disp("Size of G matriz is right")
else
    disp("Error: Size of G matriz IS NOT right")
end

if size(im.h) == [im.nModes length(im.KRadialLin)]
    disp("Size of h matriz is right")
else
    disp("Error: Size of h matriz IS NOT right")
end

if size(im.omega) == [im.nModes length(im.KRadialLin)]
    disp("Size of omega matriz is right")
else
    disp("Error: Size of omega matriz IS NOT right")
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visual check of F an G and h after interp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(2)

sgtitle("Check if most part of variance is btw turning points")

subplot(2,2,1)
plot(im.F(:,3,2),im.zNew ,'k',LineWidth=1.5) 
ylabel("Low mode")
title(["Long Wave (L=", num2str(2*pi/im.KRadialLin(2)),"m)"])

subplot(2,2,2)
plot(im.F(:,3,end),im.zNew ,'k',LineWidth=1.5)
title(["Short wave(L=", num2str(2*pi/im.KRadialLin(end)),"m)"])

subplot(2,2,3)
plot(im.F(:,end-5,2),im.zNew ,'k',LineWidth=1.5)
ylabel(["High mode(", num2str(im.nModes-5),")"])


subplot(2,2,4)
plot(im.F(:,end-5,end),im.zNew ,'k',LineWidth=1.5)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check value of coeficients (HKE, VKE and PE)
% Should they be >=0 and <=1???
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
