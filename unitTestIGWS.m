classdef unitTestIGWS < matlab.unittest.TestCase

    % InternalGravityWaveSpectrumUnitTest
    % 
    % Leticia fabre de Lima
    %
    % August, 2024   Version 2.0
    %
    % result = run(matlab.unittest.TestSuite.fromClass(?unitTestIGWS));
    % rt = table(result)
    % rt.Details{4,1}.DiagnosticRecord.Report
    %
    % import matlab.unittest.TestSuite
    % suite = TestSuite.fromMethod(?TestIGWSInitialization, 'testInitWithLatitude');
    % result = run(suite)


    properties
        im  % Property to store the InternalGravityWaveSpectrum instance
    end

    
    methods (Static)
        function arbitraryStratFunc = getArbitraryStratFunc()
            % This method returns an arbitrary stratification function based on atlas data
            lat0 = -50.0;
            lon0 = -25.0;
            atlas = VerticalModeAtlas('PhD/Reps/InternalGravityWaveSpectrum/modeAtlasFile.nc');
            [N2, z] = atlas.N2(lat0, lon0);
            arbitraryStratFunc = @(zin) interp1(z, N2, zin);  % Interpolated stratification function
        end
    end

    properties (ClassSetupParameter)
        Lz = struct('Lz', 4000);  % Vertical extent of the domain
        latitudeInit = {10, 33};  % Initialize latitude parameter
        stratification = {'exponential', 'constant', 'arbitrary'};  % Types of stratification
    end

    methods (TestClassSetup)
        function classSetup(testCase, Lz, latitudeInit, stratification)
            % Set up the InternalGravityWaveSpectrum based on stratification type
            switch stratification
                case 'exponential'
                    N2 = @(z) 3*2*pi/3600 * 3*2*pi/3600 * exp(2*z/1300);
                case 'constant'
                    N2 = @(z) 5.2e-3;
                case 'arbitrary'
                    N2 = unitTestIGWS.getArbitraryStratFunc();
            end
            testCase.im = InternalGravityWaveSpectrum(N2, Lz, 'latitude', latitudeInit);
        end
    end



methods (Test)
        function testVerticalVarianceWVM(testCase)
            % Test the vertical variance against interquartile ranges

            plot = 0;  % Set to 1 to enable plotting
            D = 4000;  % Depth
            zvect = linspace(-testCase.im.Lz, 0, 1000);

            % Initialize Wave Vortex Model (WVM) transform
            Lx = 200e3;
            Ly = Lx;  % Assuming square domain for simplicity
            Nx = 256; Ny = 256; Nz = 129;  % Grid points

            wvt = WVTransformBoussinesq([Lx, Ly, D], [Nx, Ny, Nz], 'latitude', testCase.im.latitude, 'N2', testCase.im.N2);
            wvt.initWithAlternativeSpectrum;

            % Define energy terms to evaluate
            energyTerms = {'HKE', 'VKE', 'PE', 'TE'};

            if plot
                figure(17);
            end

            % Loop over each energy term
            for i = 1:length(energyTerms)
                energyTerm = energyTerms{i};
                verticalVariance = testCase.im.verticalVariance(energyTerm, 'zVector', zvect, 'Plot', 0, 'Mask', 1);

                % Calculate energy based on term
                switch energyTerm
                    case 'PE'
                        Energy = 0.5 * (wvt.eta.^2 .* reshape(wvt.N2, [1, 1, size(wvt.eta, 3)]));
                    case 'HKE'
                        Energy = 0.5 * (wvt.u.^2 + wvt.v.^2);
                    case 'VKE'
                        Energy = 0.5 * wvt.w.^2;
                    otherwise
                        Energy = 0.5 * (wvt.u.^2 + wvt.v.^2 + wvt.w.^2 + ...
                                        (wvt.eta.^2 .* reshape(wvt.N2, [1, 1, size(wvt.eta, 3)])));
                end

                % Calculate mean and quartiles along z-dimension
                Energy_bar = squeeze(mean(Energy, [1, 2]));
                q25 = squeeze(prctile(Energy, 25, [1, 2]));
                q75 = squeeze(prctile(Energy, 75, [1, 2]));

                % Interpolate to zvect
                Energy_bar_interp = interp1(wvt.z, Energy_bar, zvect);
                q25_interp = interp1(wvt.z, q25, zvect);
                q75_interp = interp1(wvt.z, q75, zvect);

                % Unit test: Verify that vertical variance is within interquartile range
                testCase.verifyGreaterThan(verticalVariance, q25_interp, ...
                    'Vertical Variance is not above q25.');
                testCase.verifyLessThan(verticalVariance, q75_interp, ...
                    'Vertical Variance is not below q75.');

                % Plotting if enabled
                if plot
                    subplot(2, 2, i);
                    hold on;

                    % Plot quartiles
                    fill([q25_interp, fliplr(q75_interp)] * 100, [zvect, fliplr(zvect)], ...
                         'cyan', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                    plot(Energy_bar_interp * 100, zvect, 'b', 'LineWidth', 2);
                    plot(verticalVariance * 100, zvect, 'k--', 'LineWidth', 2);

                    % Label and title
                    xlabel([energyTerm, ' (cm^2s^{-2})']);
                    ylabel('Depth (m)');
                    title([energyTerm, ' with Interquartile Range']);
                    grid on;
                    hold off;
                end
            end

            if plot
                sgtitle('Energy Terms with Interquartile Range');
            end
        end
    end
end







        
 
 % 
 % 
 %        function plotModeHighestFrequency(im)
 %            % Visual check of modes on highest freq
 %            disp("Check on figure 1 if most part of variance is constrained between turning points")
 % 
 %            figure(1)
 %            subplot(1,3,1)
 %            sgtitle("Check if most part of variance is btw turning points")
 % 
 %            plot(sqrt(im.N2zInitial)*3600/(2*pi),im.zInitial ,'k',LineWidth=1.5)
 %            hold on
 %            xline(0.8*sqrt(im.N2max)*3600/(2*pi))
 %            ylabel('depth')
 %            xlabel('cph')
 %            title('N(z)')
 %            ylim([-im.Lz*0.25 0])
 % 
 %            subplot(1,3,2)
 %            plot(im.FInitial(:,end),im.zInitial ,'k',LineWidth=1.5)                      
 %            title('Initial F - Highest Mode')   
 %            ylim([-im.Lz*0.25 0])
 % 
 %            subplot(1,3,3)
 %            plot(im.GInitial(:,end),im.zInitial ,'k',LineWidth=1.5)                        
 %            title('Initial G - Highest Mode')   
 %            ylim([-im.Lz*0.25 0])   
 % 
 %        end
 % 
 % 
 %        function testEquallySpacedLinearLog(testCase)
 %            % Check if linear data is equally spaced
 %            linearDiff = diff(im.KRadialLin);
 %            linearSpacing = all(linearDiff == linearDiff(1));
 %            testCase.verifyTrue(linearSpacing, 'Linear data is not equally spaced.');
 % 
 %            % Check if log data is equally spaced
 %            logDiff = diff(log(im.KRadialLog));
 %            logSpacing = all(logDiff == logDiff(1));
 %            testCase.verifyTrue(logSpacing, 'Log data is not equally spaced.');
 % 
 %        end
 % 
 %        function plotHorizontalScaleResolution(im)
 %            %
 %            % Visual check of Horizontal wavelength distribution in linear scale
 % 
 %            figure(2)
 %            semilogy(1:length(im.KRadialLin),(2*pi./im.KRadialLin),".")
 %            title("Check the resolution for longwaves")
 %            xlim([0 0.1*length(im.KRadialLin)])
 %            xlabel("Count")
 %            ylabel("Wavelength [m]")
 %        end
 % 
 % 
 % 
 %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %        % Is the size of main matrix correct after scattered interp?
 %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % 
 %        if size(im.F) == [length(im.zNew) im.nModes length(im.KRadialLin)]
 %            disp("Size of F matriz is right")
 %        else
 %            disp("Error: Size of F matriz IS NOT right")
 %        end
 % 
 %        if size(im.G) == [length(im.zNew) im.nModes length(im.KRadialLin)]
 %            disp("Size of G matriz is right")
 %        else
 %            disp("Error: Size of G matriz IS NOT right")
 %        end
 % 
 %        if size(im.h) == [im.nModes length(im.KRadialLin)]
 %            disp("Size of h matriz is right")
 %        else
 %            disp("Error: Size of h matriz IS NOT right")
 %        end
 % 
 %        if size(im.omega) == [im.nModes length(im.KRadialLin)]
 %            disp("Size of omega matriz is right")
 %        else
 %            disp("Error: Size of omega matriz IS NOT right")
 %        end
 % 
 % 
 % 
 % 
 %        function plotFGAfterInterpolation(im)
 %            % Visual check of F an G and h after interp
 % 
 %            figure(3)
 % 
 %            sgtitle("F and G; Vertical mode = 3; Vary wavelength")
 % 
 %            subplot(1,2,1)
            % for ii = 1:5:20
            %  plot(im.F(:,3,ii),im.zNew ,'k',LineWidth=1.5) 
            %  hold on
            % % pause
            % end
            % 
            %  for ii = 1:3
            %  plot(FiK(:,3,ii),zPerModeLog(:,ii) ,'k',LineWidth=1.5) 
            %  hold on
            % % pause
            % end
 %            title("F")
 %            ylabel("Depth [m]")
 % 
 %            subplot(1,2,2)
 %            for ii = 1:5:20
 %             plot(im.G(:,3,ii),im.zNew ,'k',LineWidth=1.5) 
 %             hold on
 %            % pause
 %            end
 %            title("G")
 % 
 %            %%%%%%%%
 %            figure(4)
 % 
 %            sgtitle("F and G; Vertical mode = 12; Vary wavelength")
 % 
 %            subplot(1,2,1)
 %            for ii = 1:2:20
 %             plot(im.F(:,12,ii),im.zNew ,'k',LineWidth=1.5) 
 %             hold on
 %            % pause
 %            end
 %            title("F")
 %            ylabel("Depth [m]")
 % 
 %            subplot(1,2,2)
 %            for ii = 1:2:20
 %             plot(im.G(:,12,ii),im.zNew ,'k',LineWidth=1.5) 
 %             hold on
 %            % pause
 %            end
 %            title("G")
 % 
 %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %        % Check value of coeficients (HKE, VKE and PE)
 %        % Should they be >=0 and <=1???
 %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % end