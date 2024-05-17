classdef TestIGWSInitialization < matlab.unittest.TestCase

    % InternalGravityWaveSpectrumUnitTest
    % 
    % Leticia fabre de Lima
    %
    % April, 2024   Version 1.0
    %
    % result = run(matlab.unittest.TestSuite.fromClass(?TestIGWSInitialization));
    % rt = table(result)
    %rt.Details{4,1}.DiagnosticRecord.Report

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



    properties (TestParameter)
        latitude = {0,5,10,90,-5,-10,-90}
        %shouldForceMonotonicDensity={0,1};
        %N2zInitial={ones(1000)*5.2e-3, unitTestIGWS.getArbitraryStratFunc(1:1000),rand(1000,1)}
    end

     methods (Test)
        function testInitWithLatitude(testCase,latitude)
            if latitude >= -5 && latitude <= 5
            testCase.verifyError(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,latitude=latitude),"Latitude:MustBeAwayEquator");
            elseif latitude <= -85 || latitude > 85
                testCase.verifyError(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,latitude=latitude),"Latitude:WrongValue");
            else
                testCase.verifyWarningFree(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,latitude=latitude));
            end
        end

        
        function testValidationWhenShouldForceMonotonicDensityEnabled(testCase)            
            testCase.verifyWarningFree(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,shouldForceMonotonicDensity=1));
        end

        function testValidationWhenShouldForceMonotonicDensityDisabled(testCase)
            testCase.verifyWarningFree(@()InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,shouldForceMonotonicDensity=0));
        end

        function testValidationWithInvalidInput(testCase)
            testCase.verifyError(@()InternalGravityWaveSpectrum(@(z) randn(size(z))*5.2e-3,4000,shouldForceMonotonicDensity=1), 'MATLAB:expectedIncreasing');
        end
    end



    end

   
