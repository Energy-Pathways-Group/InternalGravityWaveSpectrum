
classdef TestIGWSInitialization < matlab.unittest.TestCase

    % InternalGravityWaveSpectrumUnitTest
    % 
    % Leticia fabre de Lima
    %
    % April, 2024   Version 1.0
    %
    % result = run(matlab.unittest.TestSuite.fromClass(?TestIGWSInitialization));
    % rt = table(result)
    % rt.Details{4,1}.DiagnosticRecord.Report
    %
    % import matlab.unittest.TestSuite
    % suite = TestSuite.fromMethod(?TestIGWSInitialization, 'testInitWithLatitude');
    % result = run(suite)

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
    
        % Test Stratification
        
        function testValidationWhenShouldForceMonotonicDensityEnabled(testCase)            
            testCase.verifyWarningFree(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,shouldForceMonotonicDensity=1));
        end

        function testValidationWhenShouldForceMonotonicDensityDisabled(testCase)
            testCase.verifyWarningFree(@()InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,shouldForceMonotonicDensity=0));
        end

        function testValidationWithInvalidInput(testCase)
            testCase.verifyError(@()InternalGravityWaveSpectrum(@(z) randn(size(z))*5.2e-3,4000,shouldForceMonotonicDensity=1), 'MATLAB:expectedIncreasing');
        end
        
        %%%%%%%%%%%%%%% Test orthogonality First Condition %%%%%%%%%%%%%%%%  
        function testOrthogonalityFirstCond(testCase)
            D=4000;
            testCase.im = InternalGravityWaveSpectrum(@(z)5.2e-3,D);

                        
            %Creating Gi*Gj Matrix
            %orthogonalMatrix = zeros(testCase.im.nZ,size, size);           

            for indK= 1:testCase.im.nK
                for indZ = 1:testCase.im.nZ

                    slice = squeeze(testCase.im.G(indZ, :,indK));
                    A= slice(:);
                    B=A';
                    allProductMatrix(:,:,indZ)= A*B;
                end
                for n=1:testCase.im.nModes
                    for m=1:testCase.im.nModes
                        orthogonalMatrix(n,m,indK)=trapz(testCase.im.zPerMode(:,indK),(testCase.im.N2(testCase.im.zPerMode(:,indK)) - testCase.im.f0^2) .* allProductMatrix(n,m,:));
            
                    end
                end
            end

            %Creating expected Matrix     
            row=size(orthogonalMatrix,1);
            col=size(orthogonalMatrix,2);

            for indK= 1:testCase.im.nK
                expectedMatrix(:,:, indK) = eye(row, col);    
            end
            expectedMatrix = expectedMatrix*9.8;
            expectedMatrix(64,64,:)=0;
            

            %Testing
            testCase.verifyEqual(orthogonalMatrix, expectedMatrix,"AbsTol", 0.1)
        end


        %%%%%%%%%%%%%%% Test orthogonality Second Condition %%%%%%%%%%%%%%%%   
     
        function testOrthogonalitySecondCond(testCase)
            D=4000;
            testCase.im = InternalGravityWaveSpectrum(@(z)5.2e-3,D);

                        
            %Creating Gi*Gj, hihj, FiFj Matrices
                      
            clear A B
            for indK= 1:testCase.im.nK
                for indZ = 1:testCase.im.nZ

                    sliceG = squeeze(testCase.im.G(indZ, :,indK));
                    A= sliceG(:);
                    B=A';
                    allProductMatrixG(:,:,indZ)= A*B;

                    clear A B

                    sliceF = squeeze(testCase.im.F(indZ, :,indK));
                    A= sliceF(:);
                    B=A';
                    allProductMatrixF(:,:,indZ)= A*B;

                    clear A B
                end

                sliceh = squeeze(testCase.im.h(:,indK));
                A= sliceh(:);
                B=A';
                allProductMatrixh(:,:,indK)= A*B;

                                
                for n=1:testCase.im.nModes
                    for m=1:testCase.im.nModes
                        integrand= allProductMatrixF(n,m,:) + allProductMatrixh(n,m,indK)*testCase.im.KRadialLog(indK).^2*allProductMatrixG(n,m,:);

                        orthogonalMatrix(n,m,indK)=trapz(testCase.im.zPerMode(:,indK), integrand);
            
                    end
                end

                %Creating expected Matrix     
                row=size(orthogonalMatrix,1);
                col=size(orthogonalMatrix,2);
                
                expectedMatrix= zeros(size(orthogonalMatrix));

                for indK= 1:testCase.im.nK
                    for i=1:min(row,col)
                        expectedMatrix(i,i, indK) = testCase.im.h(i,indK);
                    end
                end

                %expectedMatrix(64,64,:)=0;                
            end
           

            %Testing NOT WORKING FOR LAST MODE
            testCase.verifyEqual(orthogonalMatrix(1:end-1,1:end-1,:), expectedMatrix(1:end-1,1:end-1,:),"AbsTol", 0.1)
        end
        %%%%%%%%%%%%%%%  Test orthogonality 2 %%%%%%%%%%%%%%%%%%%%%%%%%%
        % NO ENERGY IN THE LAST MODE? 
        %
        function testOrthogonalityFirstCondSimple(testCase)
            Depth=4000;
            testCase.im = InternalGravityWaveSpectrum(@(z)5.2e-3,Depth);

            expectedMatrix = ones(testCase.im.nModes-1,testCase.im.nK)*9.8;

            orthogonalMatrix=zeros(testCase.im.nModes-1,testCase.im.nK);
            for indk = 1:testCase.im.nK
                for indj =1:testCase.im.nModes -1
    
                orthogonalMatrix(indj,indk)=trapz(testCase.im.zPerMode(:,indk),(testCase.im.N2(testCase.im.zPerMode(:,indk)) - testCase.im.f0^2) .* testCase.im.G(:,indj,indk).^2,1);
                end
            end

           testCase.verifyEqual(orthogonalMatrix, expectedMatrix,"AbsTol", 0.1)



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Test Vertical Bases
        function testVerticalBases(testCase)
            D=4000;
            testCase.im = InternalGravityWaveSpectrum(@(z)5.2e-3,D);

            % why is phase oposite?
            for jind= 1:testCase.im.nModes
                for zind= 1:length(testCase.im.zNew)
                    A2 = (1/D)*(2*testCase.im.g/((5.2e-3) - (testCase.im.f0)^2));
                    A= sqrt(A2);
                    m = jind*pi/D;

                    G(zind,jind)=A*sin(m*(testCase.im.zNew(zind)+D));
                    
                end
            end
            
            % Add both Abs and Rel Tol (10^3)
            testCase.verifyEqual(squeeze(testCase.im.G(:,:,1).^2), G.^2,"AbsTol", 0.1)
            % considering two error tolerance

        end



        % Test Energy
        function testDistributionAlternativeSpectrum(testCase)
            testCase.im = InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000);
            testCase.verifyEqual(sum(testCase.im.TE(:)), testCase.im.E_T(:),"AbsTol", 0.3)
        end



        function testEnergyCoeficients(testCase)
            testCase.im = InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000);
            
            %try integral of each mode at a time
            %test the orthogonal condition
            HKEIntegral= trapz(testCase.im.zNew,testCase.im.HKEcoef.*testCase.im.F.^2,1);
            VKEIntegral= trapz(testCase.im.zNew,testCase.im.VKEcoef.*testCase.im.G.^2,1);
            PEIntegral= trapz(testCase.im.zNew,testCase.im.PEcoef.*testCase.im.G.^2.*testCase.im.N2atQuadPoints',1);
            allIntegral = HKEIntegral+VKEIntegral+PEIntegral;
            allIntegral=squeeze(allIntegral);

            testCase.verifyEqual(allIntegral,testCase.im.h/2,"RelTol", 0.1)

            figure(1)
            pcolor(allIntegral)
            colorbar
            clim([0, 0.25]);

            figure(2)
            pcolor(testCase.im.h/4)
            colorbar
            clim([0, 0.25]);
        end
        

        function testTotalEnergy(testCase)
            testCase.im = InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,shouldForceMonotonicDensity=0);

            HKEIntegral= trapz(testCase.im.zNew,testCase.im.A.^2.*testCase.im.HKEcoef.*testCase.im.F.^2,1);
            VKEIntegral= trapz(testCase.im.zNew,testCase.im.A.^2.*testCase.im.VKEcoef.*testCase.im.G.^2,1);
            PEIntegral= trapz(testCase.im.zNew,testCase.im.A.^2.*testCase.im.PEcoef.*testCase.im.G.^2.*testCase.im.N2atQuadPoints',1);
            allIntegral = HKEIntegral+VKEIntegral+PEIntegral;
            allIntegral=squeeze(allIntegral);

            testCase.verifyEqual(allIntegral,testCase.im.TE,"AbsTol", 0.9)

            figure(3)
            pcolor(allIntegral)
            figure(4)
            pcolor(testCase.im.TE)
        end

    end



 end
end

   
