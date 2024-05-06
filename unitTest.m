
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
% check is not zero
if im.latitude >= -1 && im.latitude <=1
    disp("Coriolis frequency has a resonable value")
else
    disp("Error: Coriolis frequency DOES NOT have a resonable value")
end


import matlab.unittest.constraints.Throws
testCase = matlab.unittest.TestCase.forInteractiveUse;
testCase.verifyThat(@() InternalGravityWaveSpectrum(N2Func,Lz,"latitude",[-5:1:5]),Throws("Latitude:MustBeAwayEquator"))

testCase.verifyThat(@() InternalGravityWaveSpectrum(N2Func,Lz,"latitude",100),Throws("Latitude:WrongValue"))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to check if N2 is right or "smooth enough"?
% Is N2 smooth enouth?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%monotonic use Jeffrey's test internal modes

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


dOmegaVector = diff(omega);
if any(dOmegaVector<0)
    error('omega must be strictly monotonically increasing.')
end

warning('Mean stratification (N2) changes by %d orders of magnitude. This may lead to numerical instability.',round(dStrat));

validateattributes( yourVector, { 'numeric' }, { 'vector', 'increasing' } )

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

% Use internal mode exponential estrat toolbox and plot same mode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is Kradial equally spaced?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diffLin = diff(im.KRadialLin);
if sum(abs(diff(diffLin))) <= 10^-10
    disp("KRadial Linear is equally spaced in linear scale")
else
    disp("Error: KRadial Linear IS NOT equally spaced in linear scale")
end

diffLog = diff(log(im.KRadialLog));
if sum(abs(diff(diffLog))) <= 10^-10  
    disp("KRadial Log is equally spaced in log scale")
else
    disp("Error: KRadial Log IS NOT equally spaced in log scale")
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visual check of Horizontal wavelength distribution in linear scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
semilogy(1:length(im.KRadialLin),(2*pi./im.KRadialLin),".")
title("Check the resolution for longwaves")
xlim([0 0.1*length(im.KRadialLin)])
xlabel("Count")
ylabel("Wavelength [m]")

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
figure(3)

sgtitle("F and G; Vertical mode = 3; Vary wavelength")

subplot(1,2,1)
for ii = 1:5:20
 plot(im.F(:,3,ii),im.zNew ,'k',LineWidth=1.5) 
 hold on
% pause
end
title("F")
ylabel("Depth [m]")

subplot(1,2,2)
for ii = 1:5:20
 plot(im.G(:,3,ii),im.zNew ,'k',LineWidth=1.5) 
 hold on
% pause
end
title("G")

%%%%%%%%
figure(4)

sgtitle("F and G; Vertical mode = 12; Vary wavelength")

subplot(1,2,1)
for ii = 1:2:20
 plot(im.F(:,12,ii),im.zNew ,'k',LineWidth=1.5) 
 hold on
% pause
end
title("F")
ylabel("Depth [m]")

subplot(1,2,2)
for ii = 1:2:20
 plot(im.G(:,12,ii),im.zNew ,'k',LineWidth=1.5) 
 hold on
% pause
end
title("G")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check value of coeficients (HKE, VKE and PE)
% Should they be >=0 and <=1???
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
