

% Unit tests %
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

figure(1)
subplot(1,3,1)
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
% Check if Kradial is equally spaced
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check size of main matrix before/after interp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check value of coeficients (HKE, VKE and PE)
% Should they be >=0 and <=1???
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
