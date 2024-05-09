classdef unitTestIGWS < matlab.unittest.TestCase

    % InternalGravityWaveSpectrumUnitTest
    % 
    % Leticia fabre de Lima
    %
    % April, 2024   Version 1.0


    properties
        im
    end

    
    methods (Static)
        function arbitraryStratFunc = getArbitraryStratFunc()
            lat0 = -50.0;
            lon0 = -25.0;
            atlas = VerticalModeAtlas('PhD/Reps/vertical-mode-atlas/modeAtlasFile.nc');
            rho = atlas.rho(lat0, lon0);
            [N2, z] = atlas.N2(lat0, lon0);
            arbitraryStratFunc = @(zin) interp1(z, N2, zin);
        end
    end


    properties (ClassSetupParameter)
        Lz = struct('Lz',4000);
        latitudeInit= (-90:10:100); % Define latitude here if needed
        stratification = {'exponential','constant','arbitrary'};        
    end

    methods (TestClassSetup)
        function classSetup(testCase, Lz, latitudeInit,stratification)
            switch stratification
                case 'exponential'
                   N2=@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300);
                case 'constant'
                   N2=@(z)5.2e-3;      
                case 'arbitrary'                    
                   N2=MyUnitTest.getArbitraryStratFunc();         
            end
           testCase.im = InternalGravityWaveSpectrum(N2, Lz, 'latitude', latitudeInit);

        end
    end

    methods (TestParameterDefinition, Static)
        function latitude = initializeLatitudeParameter()
            latitude = (-90:10:100);
        end
    end

   properties (TestParameter)
        latidute
    end

methods (Test)
    function isLatitudeCorrect(testCase, latitude)
        % Test if the correct error is thrown based on latitude value

        % Assuming InternalGravityWaveSpectrum has a method that sets latitude
        testCase.im.setLatitude(latitude);

        % Depending on the latitude, the class should throw different errors
        if latitude < -10
            % Expecting 'MustBeAwayEquator' error for latitude less than -10
            testCase.verifyThat(@() testCase.im.someMethod(), Throws("Latitude:MustBeAwayEquator"));
        else
            % Expecting 'WrongValue' error for latitude greater than or equal to -10
            testCase.verifyThat(@() testCase.im.someMethod(), Throws("Latitude:WrongValue"));
        end
    end
end
end



        
   %%%%%%%%%%%%%%%%%% TO BE UPDATED %%%%%%%%%%%%%%%%%%%%%%%%%%%%     

 %         function testMonotonicDensityValidation(testCase)       
 %            % Mock the options structure
 %            options.shouldForceMonotonicDensity = 1;
 % 
 %            % Test with increasing vector
 %            N2zInitial = [1, 2, 3, 4];
 %            verifyError(testCase, @()validateattributes(N2zInitial, {'numeric'}, {'vector', 'increasing'}), '');
 % 
 %            % Test with decreasing vector
 %            N2zInitial = [4, 3, 2, 1];
 %            verifyError(testCase, @()validateattributes(N2zInitial, {'numeric'}, {'vector', 'increasing'}), '');
 % 
 %            % Test with non-increasing vector
 %            N2zInitial = [1, 3, 2, 4];
 %            verifyError(testCase, @()validateattributes(N2zInitial, {'numeric'}, {'vector', 'increasing'}), '');
 %        end
 % 
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
 %            for ii = 1:5:20
 %             plot(im.F(:,3,ii),im.zNew ,'k',LineWidth=1.5) 
 %             hold on
 %            % pause
 %            end
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