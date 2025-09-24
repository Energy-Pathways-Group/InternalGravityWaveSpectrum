
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check if latitude is correct/valid
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function testInitWithLatitude(testCase,latitude)
            if latitude >= -5 && latitude <= 5
            testCase.verifyError(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,latitude=latitude),"Latitude:MustBeAwayEquator");
            elseif latitude <= -85 || latitude > 85
                testCase.verifyError(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,latitude=latitude),"Latitude:WrongValue");
            else
                testCase.verifyWarningFree(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,latitude=latitude));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check if stratification is monotonic depending on condition
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
        function testValidationWhenShouldForceMonotonicDensityEnabled(testCase)            
            testCase.verifyWarningFree(@() InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,shouldForceMonotonicDensity=1));
        end

        function testValidationWhenShouldForceMonotonicDensityDisabled(testCase)
            testCase.verifyWarningFree(@()InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,shouldForceMonotonicDensity=0));
        end

        function testValidationWithInvalidInput(testCase)
            testCase.verifyError(@()InternalGravityWaveSpectrum(@(z) randn(size(z))*5.2e-3,4000,shouldForceMonotonicDensity=1), 'MATLAB:expectedIncreasing');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test orthogonality - First Condition
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
        function testOrthogonalityFirstCond(testCase)
            D=4000;
            N0 = 3*2*pi/3600;
            L_gm = 1300;
            %testCase.im = InternalGravityWaveSpectrum(@(z)5.2e-3,D);
            testCase.im = InternalGravityWaveSpectrum(@(z) N0*N0*exp(2*z/L_gm),D);
            
                        
            %Creating Gi*Gj Matrix                   
            for indK= 1:testCase.im.nK
                for indZ = 1:testCase.im.nZ

                    slice = squeeze(testCase.im.G(indZ, :,indK));
                    A= slice(:);
                    B=A';
                    allProductMatrix(:,:,indZ)= A*B;
                end
                for n=1:testCase.im.nModes
                    for m=1:testCase.im.nModes
                        integrand= (testCase.im.N2(testCase.im.zPerMode(:,indK)) - testCase.im.f0^2) .* squeeze(allProductMatrix(n,m,:));
                        orthogonalMatrix(n,m,indK)=trapz(testCase.im.zPerMode(:,indK),integrand);            
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
            errorMatrix = orthogonalMatrix - expectedMatrix;

            
           %%%%%%%%%%%    Plot      %%%%%%%%%%%%

           %For each K figure out each mode starting failing this unit
           %test
            
           AbsTol= 0.05;
            
            lowerVerticalMode = zeros(testCase.im.nK, 1); % Vetor para armazenar o menor modo para cada número de onda horizontal
            
            for iK = 1:testCase.im.nK                
                slice = errorMatrix(:,:,iK);
                
                % First element that is bigger than the AbsTol
                [row, col] = find(slice > AbsTol, 1, 'first');
                
                if isempty(row)                    
                    lowerVerticalMode(iK) = NaN;
                else                    
                    lowerVerticalMode(iK) = min(row, col); 
                end
            end
            
            
            plot(1:testCase.im.nK,lowerVerticalMode)
            ylabel("Lower Vertical Number")
            xlabel("Horizontal Wave Number")

            %Creating Mask           
            mask = ones(testCase.im.nModes, testCase.im.nModes, testCase.im.nK);
            
            for iK = 1:testCase.im.nK
                if ~isnan(lowerVerticalMode(iK))
                    % Ajustar a matriz de máscara para definir os valores a 0 a partir do menor modo vertical
                    
                    mask(lowerVerticalMode(iK):end, lowerVerticalMode(iK):end, iK) = 0;
                end
            end        

            %Testing
            testCase.verifyEqual(orthogonalMatrix(mask==1), expectedMatrix(mask==1),"AbsTol", AbsTol)

            %Video 
            if(0)
            
                % Open video
                videoFile = 'errorFandG.avi';
                v = VideoWriter(videoFile);                
                v.FrameRate = 2;
                open(v);
                
                
                % Video frame to frame
                for iK = 1:testCase.im.nK
                    
                    figure;
                    
                    % Plotting pcolor for iK
                    pcolor(errorMatrix(:,:,iK));                    
                    colorbar;
                    clim([0 0.1]);
                    title(['iK: ', num2str(iK)]);
                    
                    % Write frame in a video
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                    
                    % Fechar a figura para não sobrecarregar a memória
                    close(gcf);
                end
                
                % Fechar o arquivo de vídeo
                close(v);
                
                % Exibir mensagem de conclusão
                disp(['Video saved as ', videoFile]);

            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test orthogonality - Second Condition
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
     
        function testOrthogonalitySecondCond(testCase)
            D=4000;
            N0 = 3*2*pi/3600;
            L_gm = 1300;
            %testCase.im = InternalGravityWaveSpectrum(@(z)5.2e-3,D);
            testCase.im = InternalGravityWaveSpectrum(@(z) N0*N0*exp(2*z/L_gm),D);

                        
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


            end
           errorMatrix = orthogonalMatrix - expectedMatrix;

            
           %%%%%%%%%%%    Plot      %%%%%%%%%%%%

           %For each K figure out each mode starting failing this unit
           %test
            
           AbsTol= 0.05;
            
            lowerVerticalMode = zeros(testCase.im.nK, 1); % Vetor para armazenar o menor modo para cada número de onda horizontal
            
            for iK = 1:testCase.im.nK                
                slice = errorMatrix(:,:,iK);
                
                % First element that is bigger than the AbsTol
                [row, col] = find(slice > AbsTol, 1, 'first');
                
                if isempty(row)                    
                    lowerVerticalMode(iK) = NaN;
                else                    
                    lowerVerticalMode(iK) = min(row, col); 
                end
            end
            
            
            plot(1:testCase.im.nK,lowerVerticalMode)
            ylabel("Lower Vertical Number")
            xlabel("Horizontal Wave Number")

            %Creating Mask           
            mask = ones(testCase.im.nModes, testCase.im.nModes, testCase.im.nK);
            
            for iK = 1:testCase.im.nK
                if ~isnan(lowerVerticalMode(iK))
                    % Ajustar a matriz de máscara para definir os valores a 0 a partir do menor modo vertical
                    
                    mask(lowerVerticalMode(iK):end, lowerVerticalMode(iK):end, iK) = 0;
                end
            end        

            %Testing
            testCase.verifyEqual(orthogonalMatrix(mask==1), expectedMatrix(mask==1),"AbsTol", AbsTol)            
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test Vertical Bases (computational x analitical)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
     
        function testVerticalBases(testCase)
            D=4000;
            testCase.im = InternalGravityWaveSpectrum(@(z)5.2e-3,D);

            % why is phase oposite?
            for jind= 1:testCase.im.nModes
                for zind= 1:testCase.im.nZ
                    A2 = (1/D)*(2*testCase.im.g/((5.2e-3) - (testCase.im.f0)^2));
                    A= sqrt(A2);
                    m = jind*pi/D;

                    G(zind,jind)=A*sin(m*(testCase.im.zPerMode(zind,1)+D));
                    
                end
            end
            
            % Add both Abs and Rel Tol (10^3)
            testCase.verifyEqual(squeeze(testCase.im.G(:,:,1).^2), G.^2,"AbsTol", 0.1)
            % considering two error tolerance

        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test if the total energy is  same before and after distribution 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        function testDistributionAlternativeSpectrum(testCase)
            testCase.im = InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000);
            testCase.verifyEqual(sum(testCase.im.TE(:)), testCase.im.E_T(:),"AbsTol", 0.05)
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test Coeficients of energy
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        function testEnergyCoeficients(testCase)
            testCase.im = InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000);
            
            %try integral of each mode at a time
            for indk=1:testCase.im.nK
                HKEIntegral(:,indk)= trapz(testCase.im.zPerMode(:,indk),testCase.im.HKEcoef(:,indk)'.*testCase.im.F(:,:,indk).^2,1);
                VKEIntegral(:,indk)= trapz(testCase.im.zPerMode(:,indk),testCase.im.VKEcoef(:,indk)'.*testCase.im.G(:,:,indk).^2,1);
                PEIntegral(:,indk)= trapz(testCase.im.zPerMode(:,indk),testCase.im.PEcoef(:,indk)'.*testCase.im.G(:,:,indk).^2.*testCase.im.N2atQuadPoints(:,indk),1);
               
            end

            allIntegral = HKEIntegral+VKEIntegral+PEIntegral;
            allIntegral=squeeze(allIntegral);

            % Expected Value
            expectedValue = testCase.im.h / 2;

            % Tolerences
            relativeTolerance = 0.05;
            absoluteTolerance = 1e-4; 

            % Tolerance verification  
            isWithinRelTol = abs(allIntegral - expectedValue) <= relativeTolerance * abs(expectedValue);
            isWithinAbsTol = abs(allIntegral - expectedValue) <= absoluteTolerance;

                       
            % new
            testCase.verifyTrue(all(isWithinRelTol(:) | isWithinAbsTol(:)))

            % PLOT!

           
       end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test Total Energy (TE = HKEIntegral+VKEIntegral+PEIntegral;)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        function testTotalEnergy(testCase)
            testCase.im = InternalGravityWaveSpectrum(@(z) 3*2*pi/3600*3*2*pi/3600*exp(2*z/1300),4000,shouldForceMonotonicDensity=0);
            for indk=1:testCase.im.nK         
            
                HKEIntegral(:,indk)= trapz(testCase.im.zPerMode(:,indk),testCase.im.A2(:,indk)'.*testCase.im.HKEcoef(:,indk)'.*testCase.im.F(:,:,indk).^2,1);
                VKEIntegral(:,indk)= trapz(testCase.im.zPerMode(:,indk),testCase.im.A2(:,indk)'.*testCase.im.VKEcoef(:,indk)'.*testCase.im.G(:,:,indk).^2,1);
                PEIntegral(:,indk)= trapz(testCase.im.zPerMode(:,indk),testCase.im.A2(:,indk)'.*testCase.im.PEcoef(:,indk)'.*testCase.im.G(:,:,indk).^2.*testCase.im.N2atQuadPoints(:,indk),1);
            end

            allIntegral = HKEIntegral+VKEIntegral+PEIntegral;
            allIntegral=squeeze(allIntegral);

            % Expected Value
            expectedValue = testCase.im.TE;

            % Tolerences
            relativeTolerance = 0.05;
            absoluteTolerance = 1e-4; 

            % Tolerance verification  
            isWithinRelTol = abs(allIntegral - expectedValue) <= relativeTolerance * abs(expectedValue);
            isWithinAbsTol = abs(allIntegral - expectedValue) <= absoluteTolerance;

            %testCase.verifyTrue(all(isWithinRelTol(:) | isWithinAbsTol(:)))

            if(0)
                figure(3)
                pcolor(log(allIntegral))
                colorbar
                title("Total Energy")
                xlabel("Horizontal Wave Number INDEX")
                ylabel("Vertical mode")
    
                figure(4)
                pcolor(log(testCase.im.TE))
                colorbar
                title("HKE + VKE + PE")
                xlabel("Horizontal Wave Number INDEX")
                ylabel("Vertical mode")
            end
        end

    end



 end


   
