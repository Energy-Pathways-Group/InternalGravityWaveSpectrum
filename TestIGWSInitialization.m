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


    properties (TestParameter)
        latitude = {0,5,10,90}
    end

    methods (Test)
        function testInitWithLatitude(testCase,latitude)
            if latitude <= 5
            testCase.verifyError(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,latitude=latitude),"Latitude:MustBeAwayEquator");
            elseif latitude > 85
                testCase.verifyError(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,latitude=latitude),"Latitude:WrongValue");
            else
                testCase.verifyWarningFree(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,latitude=latitude));
            end
        end


    end

end

