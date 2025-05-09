classdef InternalGravityWaveSpectrum < handle
    properties (Access = public)
        latitude % Latitude for which the modes are being computed.
        f0 % Coriolis parameter at the above latitude.
        Lz % Depth of the ocean.
        N2 %function_handle
        g
        shouldForceMonotonicDensity
        N2max

        nModes, nK, nZ
        
        KRadialLog  % size(k) = nK
        j

        F  % size(F_k) = [nZ,nModes,nK]
        G  % size(G_k) = [nZ,nModes,nK]
        h  % size(h_k) = [nModes,nK]
        omega % size(omega_k) = [nK,nModes]       
        zPerMode % [nZ,nModes]        
        N2atQuadPoints
        Lr2
                  
        E_T 

        HKEcoef, VKEcoef, PEcoef

        Am, Ap, A2

        HKE, VKE, PE, TE            


        zNew

        cutoff_modes, cutoff_k

        M, B

        
    end
        properties (Access = private, Hidden)
            %Proprieties for test:
            FInitial % [nZ,nModes]
            GInitial % [nZ,nModes]
            zInitial % [nZ]
            N2zInitial
                   
        end

   

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalGravityWaveSpectrum(N2,Lz,options)
        
            arguments
                N2 %handle function
                Lz (1,1) {mustBePositive}
                options.latitude (1,1) double = 33 %set condition. How to modify erro mesage? costume validator
                options.nModes (1,1) double = 64
                options.nK (1,1) double = 64           
                options.shouldForceMonotonicDensity {mustBeNumericOrLogical} = 0
            end

           
            if options.latitude >= -5 && options.latitude <= 5
                error("Latitude:MustBeAwayEquator","This toolbox does not provide a good solution near the Equator (5°S to 5°N)")
            
            elseif options.latitude >= 85 || options.latitude <= -85
                error("Latitude:WrongValue","Latitude not valid")
            end
            

            
            self.N2=N2;  
            self.latitude=options.latitude;          
            self.nModes=options.nModes;
            self.nK=options.nK;
            self.nZ=options.nModes + 1;
            self.g=9.80665;
            self.Lz=Lz;
            self.shouldForceMonotonicDensity=options.shouldForceMonotonicDensity;
            self.j=1:self.nModes;

            % Calculate the cutoff for the last one-third of the second and third dimensions
            self.cutoff_modes = ceil(self.nModes * 2/3);
            self.cutoff_k = ceil(self.nK * 2/3);
                

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 1: Computation of min and max Kh based on the
        % stratification and latitude
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Step 1.1: Compute the Coriolis frequency based on the latitude
        f0 = 2* 7.2921*10^-5 * sind(self.latitude);
        self.f0=f0;

        % Step 1.2: Get the max bouyancy frequency based on the
        % stratification function

        % Step 1.2.1: define vertical vector (z) based on Lz
        % and nZ

        zInitial = linspace(-Lz,0,10001);
        N2zInitial= self.N2(zInitial);
        N2max = max(N2zInitial);


       if self.shouldForceMonotonicDensity == 1 
           validateattributes( N2zInitial, { 'numeric' }, { 'vector', 'increasing' } )
       end


        self.zInitial=zInitial;
        self.N2zInitial=N2zInitial;
        self.N2max=N2max;
       
                    
        % Step 1.3: Compute the K associated with max(N2)   

        im = InternalModesWKBSpectral(N2=self.N2,zIn=[-Lz 0],zOut=zInitial,latitude=self.latitude,nModes=self.nModes);       


        [FInitial,GInitial,h,k] = im.ModesAtFrequency(0.95*sqrt(N2max));
        Kmax= max(k);

        self.FInitial = FInitial;
        self.GInitial = GInitial;

        % Step 1.4: Define KRadial based on Kmin=0, Kmax and nK        
        minOrder = 0; 
        % minOrder = floor(log10(2*pi/Kmax));
        % if minOrder<=0
        %      minOrder=0;
        % end

        % KRadial equally spaced in log scale
        wavelengthLog=logspace(minOrder,6,self.nK);     
        KRadialLog=fliplr((2*pi)./wavelengthLog);
        self.KRadialLog = KRadialLog;
        self.KRadialLog(1)=0.5*self.KRadialLog(1);

        % % KRadial from frequency
        % omegaVector = linspace(self.f0*1.005,0.8*sqrt(self.N2max),32);
        % 
        % KVector=[];
        % for i= 1:length(omegaVector)
        %     [~,~,~,K] = im.ModesAtFrequency(omegaVector(i));
        %     KVector=[KVector;K];        
        % end
        % %KRadialLog=sort(KVector(:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        % Step 2: Computation of F and G Matriz [nK, nZ, nModes]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        upperBoundary = UpperBoundary.rigidLid;
        normalization = Normalization.kConstant;

        for iK=1:length(KRadialLog)    

            im = InternalModesSpectral(N2=self.N2,zIn=[-Lz 0],zOut=zInitial,latitude=self.latitude,nModes=self.nModes);


            im.normalization = normalization;
            im.upperBoundary = upperBoundary;  


            zPerModeLog(:,iK) = im.GaussQuadraturePointsForModesAtWavenumber(self.nModes+1,KRadialLog(iK));

            im = InternalModesSpectral(N2=self.N2,zIn=[-Lz 0],zOut=zPerModeLog(:,iK),latitude=self.latitude,nModes=self.nModes);
            [FiK(:,:,iK),GiK(:,:,iK),hiK(:,iK),omegaiK(:,iK)] = im.ModesAtWavenumber(KRadialLog(iK)); %modes at quadrature points and not equally spaced

        end

            %%%%%%%% IMPORTANT!!!!! %%%%%%%%%%%%%

            % the format for FiK and GiK is [depth (nZ), vertical modes (nModes), horiz wavenumber (nK)]
            % the format for hiK and omegaiK is [vertical modes (nModes), horiz wavenumber (nK)]
            % the format for hiK and omegaiK is [vertical modes (nModes), horiz wavenumber (nK)]


 
        self.zPerMode =zPerModeLog;
        self.F = FiK;
        self.G = GiK;
        self.h = hiK;
        self.omega = omegaiK;        
        self.zPerMode =zPerModeLog;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 3: Computation of the energy coeficients based on 
        % the squared equations (Jeffrey's paper)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % 
        % % Step 1: Compute coeficients in 2D [nModes,nKLin]

        self.HKEcoef = (1/4)*(1+ (f0^2./(self.omega.^2))); %must be at least 0.25
        self.VKEcoef= (1/4)* ((self.KRadialLog).^2 .* self.h.^2);        
        self.PEcoef =  (1/4)* (((self.KRadialLog).^2.*self.h.^2)./self.omega.^2);    


        % Make a separate function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 4: Computation of energy distribution according 
        % with the alternative Internal Wave Spectrum
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %self = HorizontalWavenumberGMSpectrum(self);
        %self = NonSeparableSpectrum(self); 
        %self = SeparableSpectrum(self);

      end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Spectral functions
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function selfUpdated = HorizontalWavenumberGMSpectrum(self,parameters)
        
            arguments
                self
                parameters.j_star = 3;
                parameters.slope = 3;
            end
        
            % -------------------------
            % Parameters
            % -------------------------
            j_star = 3;
            slope = 3/2;
            GMAmplitude = 1;
        
            % GM Parameters
            L_gm = 1.3e3;         % Thermocline exponential scale [m]
            invT_gm = 5.2e-3;     % Reference buoyancy frequency [rad/s]
            E_gm = 6.3e-5;        % Non-dimensional energy parameter
        
            % Total GM energy
            E_T = (L_gm^3) * (invT_gm^2) * E_gm * GMAmplitude;
            self.E_T = E_T;
        
            % -------------------------
            % Vertical mode weighting function M(j)
            % -------------------------
            M_raw = @(j) (j_star^2 + j.^2).^(-5/4);
            M_norm = sum(M_raw(1:1024));
            M = @(j) M_raw(j) / M_norm;
            self.M = M;
        
            % Sanity check: ensure M sums to 1
            % sum(M(1:1024))
        
            % -------------------------
            % Compute Rossby radius of deformation
            % -------------------------
            Lr2 = (self.g .* self.h) ./ (self.f0 ^ 2);
            self.Lr2 = Lr2;
        
            % -------------------------
            % Define unnormalized B(k,j)
            % -------------------------
            B_unnormalized = @(k, jInd) ...
                (1 ./ (k.^2 .* (interp1(self.KRadialLog, Lr2(jInd,:), k) + 1)).^slope) .* ...
                 sqrt(interp1(self.KRadialLog, Lr2(jInd,:), k));
        
            % Integration bounds
            kmin = min(self.KRadialLog(1));
            kmax = max(self.KRadialLog(end));
            %kmin=10^(-13);
            %kmax=1;
        
            % -------------------------
            % Normalize B(k,j)
            % -------------------------
            B_norm = ones(self.nModes, 1);
            for jIdx = 1:length(self.j)
                B_norm(jIdx) = integral(@(k) B_unnormalized(k, self.j(jIdx)), kmin, kmax);
            end
        
            B = @(k, jInd) ...
                ((1 ./ (k.^2 .* (interp1(self.KRadialLog, Lr2(jInd,:), k) + 1)).^slope) .* ...
                 sqrt(interp1(self.KRadialLog, Lr2(jInd,:), k)))./B_norm(jInd);
                self.B = B;
        
            % -------------------------
            % Sanity check: confirm integrals are normalized
            test_integrals = zeros(self.nModes , 1);
            for jIdx = 1:self.nModes 
                test_integrals(jIdx) = integral(@(k) B(k, self.j(jIdx)), kmin, kmax);
            end
            % Uncomment to view:
             disp('Sanity check integrals:'), disp(test_integrals)
        
            % -------------------------
            % Define model spectrum
            % -------------------------
            model_spectrum = @(k,jInd) E_T * B(k, jInd) * M(jInd);
        
            % -------------------------
            % Compute total energy (A2) using spectrum
            % -------------------------
            TE = amplitudesWithSpectrum(self, model_spectrum);
            self.A2 = 2 * TE ./ self.h;
        
            % -------------------------
            % Compute energy components
            % -------------------------
            N2atQuadPoints = self.N2(self.zPerMode);
        
            HKE = shiftdim(self.A2 .* self.HKEcoef, -1) .* self.F.^2;
            VKE = shiftdim(self.A2 .* self.VKEcoef, -1) .* self.G.^2;
        
            if isscalar(N2atQuadPoints)
                PE = shiftdim(self.A2 .* self.PEcoef, -1) .* self.G.^2 .* N2atQuadPoints;
            else
                PE = shiftdim(self.A2 .* self.PEcoef, -1) .* self.G.^2 .* ...
                     reshape(N2atQuadPoints, [self.nZ, 1, self.nK]);
            end
        
            % -------------------------
            % Store results in object
            % -------------------------
            self.HKE = HKE;
            self.VKE = VKE;
            self.PE = PE;
            self.TE = TE;
            self.N2atQuadPoints = N2atQuadPoints;
        
            selfUpdated = self;
        end

   

        function selfUpdated = NonSeparableSpectrum(self)
            arguments
                self
            end
        
            % -------------------------
            % Parameters
            % -------------------------
            j_star = 3;
            slope = 1;
            GMAmplitude = 1;
        
            % GM Parameters
            L_gm = 1.3e3;           % Thermocline exponential scale [m]
            invT_gm = 5.2e-3;       % Reference buoyancy frequency [rad/s]
            E_gm = 6.3e-5;          % Non-dimensional energy parameter
        
            % Total GM energy
            E_T = (L_gm^3) * (invT_gm^2) * E_gm * GMAmplitude;
            self.E_T = E_T;
        
            % -------------------------
            % Compute Rossby radius of deformation
            % -------------------------
            Lr2 = (self.g .* self.h) ./ (self.f0 ^ 2);
            self.Lr2 = Lr2;
            k_star2 = 1 ./ Lr2(3,:); % Using mode-3 deformation radius (Leticia 04/02/2025)
        
            % -------------------------
            % Define unnormalized spectral shape B(k,j)
            % -------------------------
            B_unnormalized = @(k, jInd) 1 ./ ...
                (k.^2 + ...
                (interp1(self.KRadialLog, Lr2(jInd,:), k)).^(-1) + ...
                interp1(self.KRadialLog, k_star2, k)).^slope;
        
            % Integration bounds
            kmin = min(self.KRadialLog(1));
            kmax = max(self.KRadialLog(end));
        
            % -------------------------
            % Normalize B(k,j)
            % -------------------------
            B_norm = ones(self.nModes, 1);
            for jIdx = 1:length(self.j)
                B_norm(jIdx) = integral(@(k) B_unnormalized(k, self.j(jIdx)), kmin, kmax);
            end
        
            % Normalized B(k,j)
            B = @(k, jInd) B_unnormalized(k, jInd) ./ B_norm(jInd);
            self.B = B;

            % -------------------------
            % Sanity check: Ensure normalization
            % -------------------------
            test_integrals = zeros(self.nModes, 1);
            for jIdx = 1:self.nModes
                test_integrals(jIdx) = integral(@(k) B(k, self.j(jIdx)), kmin, kmax);
            end
            % Uncomment the next line if you want to see the test output
            % disp('Sanity check integrals:'), disp(test_integrals)
        
            % -------------------------
            % Define model spectrum
            % -------------------------
            model_spectrum = @(k, j) (E_T / self.nModes) * B(k, j);
        
            % -------------------------
            % Compute total energy (A2) using spectrum
            % -------------------------
            TE = amplitudesWithSpectrum(self, model_spectrum);
            self.A2 = 2 * TE ./ self.h;
        
            % -------------------------
            % Compute energy components
            % -------------------------
            N2atQuadPoints = self.N2(self.zPerMode);
        
            HKE = shiftdim(self.A2 .* self.HKEcoef, -1) .* self.F.^2;
            VKE = shiftdim(self.A2 .* self.VKEcoef, -1) .* self.G.^2;
        
            if isscalar(N2atQuadPoints)
                PE = shiftdim(self.A2 .* self.PEcoef, -1) .* self.G.^2 .* N2atQuadPoints;
            else
                PE = shiftdim(self.A2 .* self.PEcoef, -1) .* self.G.^2 .* ...
                     reshape(N2atQuadPoints, [self.nZ, 1, self.nK]);
            end
        
            % -------------------------
            % Store results in object
            % -------------------------
            self.HKE = HKE;
            self.VKE = VKE;
            self.PE = PE;
            self.TE = TE;
            self.N2atQuadPoints = N2atQuadPoints;
        
            selfUpdated = self;
        end




        function selfUpdated = SeparableSpectrum(self)
            arguments
                self
            end
        
            % -------------------------
            % Parameters
            % -------------------------
            j_star = 3;
            slopeK = 1;
            slopej=1;
            GMAmplitude = 1;
        
            % GM Parameters
            L_gm = 1.3e3;           % Thermocline exponential scale [m]
            invT_gm = 5.2e-3;       % Reference buoyancy frequency [rad/s]
            E_gm = 6.3e-5;          % Non-dimensional energy parameter
        
            % Total GM energy
            E_T = (L_gm^3) * (invT_gm^2) * E_gm * GMAmplitude;
            self.E_T = E_T;
        
            % -------------------------
            % Compute Rossby radius of deformation
            % -------------------------
            Lr2 = (self.g .* self.h) ./ (self.f0 ^ 2);
            self.Lr2 = Lr2;
            k_star2 = 1 ./ Lr2(3,:); % Using mode-3 deformation radius (Leticia 04/02/2025)
        
            % -------------------------
            % Define unnormalized spectral shape B(k,j)
            % -------------------------
            B_unnormalized = @(k, jInd) 1 ./ ...
               ( (k.^2 + interp1(self.KRadialLog, k_star2, k)).^(slopeK).*...
                    ((interp1(self.KRadialLog, Lr2(jInd,:), k)).^(-2) + ...
                    (interp1(self.KRadialLog, Lr2(3,:), k)).^(-2)).^(slopej) );
                
               
        
            % Integration bounds
            kmin = min(self.KRadialLog(1));
            kmax = max(self.KRadialLog(end));
        
            % -------------------------
            % Normalize B(k,j)
            % -------------------------
            B_norm = ones(self.nModes, 1);
            for jIdx = 1:length(self.j)
                B_norm(jIdx) = integral(@(k) B_unnormalized(k, self.j(jIdx)), kmin, kmax);
            end
        
            % Normalized B(k,j)
            B = @(k, jInd) B_unnormalized(k, jInd) ./ B_norm(jInd);
            self.B = B;

            % -------------------------
            % Sanity check: Ensure normalization
            % -------------------------
            test_integrals = zeros(self.nModes, 1);
            for jIdx = 1:self.nModes
                test_integrals(jIdx) = integral(@(k) B(k, self.j(jIdx)), kmin, kmax);
            end
            % Uncomment the next line if you want to see the test output
            % disp('Sanity check integrals:'), disp(test_integrals)
        
            % -------------------------
            % Define model spectrum
            % -------------------------
            model_spectrum = @(k, j) (E_T / self.nModes) * B(k, j);
        
            % -------------------------
            % Compute total energy (A2) using spectrum
            % -------------------------
            TE = amplitudesWithSpectrum(self, model_spectrum);
            self.A2 = 2 * TE ./ self.h;
        
            % -------------------------
            % Compute energy components
            % -------------------------
            N2atQuadPoints = self.N2(self.zPerMode);
        
            HKE = shiftdim(self.A2 .* self.HKEcoef, -1) .* self.F.^2;
            VKE = shiftdim(self.A2 .* self.VKEcoef, -1) .* self.G.^2;
        
            if isscalar(N2atQuadPoints)
                PE = shiftdim(self.A2 .* self.PEcoef, -1) .* self.G.^2 .* N2atQuadPoints;
            else
                PE = shiftdim(self.A2 .* self.PEcoef, -1) .* self.G.^2 .* ...
                     reshape(N2atQuadPoints, [self.nZ, 1, self.nK]);
            end
        
            % -------------------------
            % Store results in object
            % -------------------------
            self.HKE = HKE;
            self.VKE = VKE;
            self.PE = PE;
            self.TE = TE;
            self.N2atQuadPoints = N2atQuadPoints;
        
            selfUpdated = self;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Fitting Spectral functions 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [energySpectrumModel,self] = gmSpectrumFunctionHandle(self,params)
        
            arguments
                self
                params 
                %params should be [j_star,slope,GMAmplitude]
                %params.j_star = 3;
                %params.slope = 3/2;
                %GMAmplitude = 1;
            end
        
            % -------------------------
            % Parameters
            % -------------------------
            j_star = params(1);
            slope = params(2);
            GMAmplitude = params(3);
        
            % GM Parameters
            L_gm = 1.3e3;         % Thermocline exponential scale [m]
            invT_gm = 5.2e-3;     % Reference buoyancy frequency [rad/s]
            E_gm = 6.3e-5;        % Non-dimensional energy parameter
        
            % Total GM energy
            E_T = (L_gm^3) * (invT_gm^2) * E_gm * GMAmplitude;
            self.E_T = E_T;
        
            % -------------------------
            % Vertical mode weighting function M(j)
            % -------------------------
            M_raw = @(j) (j_star^2 + j.^2).^(-5/4);
            M_norm = sum(M_raw(1:1024));
            M = @(j) M_raw(j) / M_norm;
            self.M = M;
        
            % Sanity check: ensure M sums to 1
            % sum(M(1:1024))
        
            % -------------------------
            % Compute Rossby radius of deformation
            % -------------------------
            Lr2 = (self.g .* self.h) ./ (self.f0 ^ 2);
            self.Lr2 = Lr2;
        
            % -------------------------
            % Define unnormalized B(k,j)
            % -------------------------
            B_unnormalized = @(k, jInd) ...
                (1 ./ (k.^2 .* (interp1(self.KRadialLog, Lr2(jInd,:), k) + 1)).^slope) .* ...
                 sqrt(interp1(self.KRadialLog, Lr2(jInd,:), k));
        
            % Integration bounds
            kmin = min(self.KRadialLog(1));
            %kmin=1.2566e-05;
            kmax = max(self.KRadialLog(end));
            %kmax=3.5186e-04;
        
            % -------------------------
            % Normalize B(k,j)
            % -------------------------
            B_norm = ones(self.nModes, 1);
            for jIdx = 1:length(self.j)
                B_norm(jIdx) = integral(@(k) B_unnormalized(k, self.j(jIdx)), kmin, kmax);
            end
        
            % B = @(k, jInd) ...
            %     (1 ./ (k.^2 .* (interp1(self.KRadialLog, Lr2(jInd,:), k) + 1)).^slope) .* ...
            %      sqrt(interp1(self.KRadialLog, Lr2(jInd,:), k)) / ...
            %     B_norm(jInd);
            B = @(k, jInd) B_unnormalized(k,jInd)/B_norm(jInd);
            self.B = B;
        
            % -------------------------
            % Sanity check: confirm integrals are normalized
            test_integrals = zeros(self.nModes - 1, 1);
            for jIdx = 1:self.nModes - 1
                test_integrals(jIdx) = integral(@(k) B(k, self.j(jIdx)), kmin, kmax);
            end
            % Uncomment to view:
            % disp('Sanity check integrals:'), disp(test_integrals)
            % disp('Sanity check total integrals:'), disp(sum(test_integrals))
        
            % -------------------------
            % Define model spectrum
            % -------------------------
            energySpectrumModel = @(k,jInd) self.E_T * self.B(k, jInd) * self.M(jInd);      

            % Sanity check: confirm integrals are normalized
            test_integrals = zeros(self.nModes - 1, 1);
            for jIdx = 1:self.nModes - 1
                test_integrals(jIdx) = integral(@(k) energySpectrumModel(k, self.j(jIdx)), kmin, kmax);
            end
            % Uncomment to view:
            % disp('Sanity check integrals:'), disp(test_integrals)
            disp('Sanity check total integrals:'), disp(sum(test_integrals))
        
            selfUpdated = self;
        end


        function val = objectiveFunction(self,params,k,j,E_model)
            arguments
                self
                params
                k
                j
                E_model
            end
            energySpectrum = gmSpectrumFunctionHandle(self,params);
            limk=round((length(k) * (2/3)));
            limj=round((length(j) * (2/3)));
            dk = (k(2)-k(1));

            for indj = 1:limj-1
                for indk = 1:limk
                    Ekj(indj,indk) = energySpectrum(k(indk)+(dk/1.9),j(indj+1));
                end
            end
            E_model=E_model(2:limj,1:limk);
            %Ekj = energySpectrum(k,j);
            % +1/2 dk

            if(1)
                % First pcolor plot
                subplot(1, 2, 1);
                jpcolor(log(E_model));                     
                colorbar;
                %clim([-18 -2])
                title('Model');                
                
                % Second pcolor plot
                subplot(1, 2, 2);
                jpcolor(log(Ekj));
                colorbar;
                %clim([-18 -2])
                title('Toolbox');
            end



            val = mean((E_model(:) - Ekj(:)).^2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Displaying energy distribution
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function verticalVariance =  verticalVariance(self,energyTerm,options)
            %
            arguments
                self
                energyTerm %options are: 'TE','HKE','VKE' and 
                options.zVector = linspace(-self.Lz,0,1000)   
                options.mask logical = false
                options.plot logical = true
                
            end
            
              
            if strcmp(energyTerm, 'TE')
                energy = self.HKE + self.VKE +self.PE;
            elseif strcmp(energyTerm, 'HKE')
                energy = self.HKE;
            elseif strcmp(energyTerm, 'VKE')
                energy = self.VKE;
            else
                energy = self.PE;
            end

            if options.mask == 1
                
                energyMask=energy(:, 1:self.cutoff_modes, 1:self.cutoff_k);
                KRadialLogMask=self.KRadialLog(1:self.cutoff_k);
                verticalModeMask=1:self.cutoff_modes;

                energyInterp = scatteredInterpolation(self, energyMask, options.zVector, KRadialLogMask, verticalModeMask,'mask',options.mask); 
            else
                energyInterp = scatteredInterpolation(self, energy, options.zVector, self.KRadialLog, 1:self.nModes); 
            end                       


            %Sum over modes and k
            verticalVariance=sum(sum(energyInterp,3),2);

            %%% plot %%%%
            if options.plot ==1

                fig = figure(10);

                plot(verticalVariance*1000,options.zVector)
                title("HKE")
                ylabel("Depth [m]")
                xlabel("Variance [cm^2/s^2]")
                grid on
            else
            end
        end


    
        %%%%%%%%%%%%%%%
        % I need to fix the part of mask
        function energyAtHorizontalWavenumber = energyAtHorizontalWavenumber(self, z,energyTerm, options)
                    % 
            arguments
                self
                z (1,1) double %or vector
                energyTerm %options are: 'TE','HKE','VKE' and 'PE'
                options.KRadial double = []
                options.mask logical = false
                options.plot logical = true

            end


            if strcmp(energyTerm, 'TE')
                energy = self.HKE + self.VKE +self.PE;
            elseif strcmp(energyTerm, 'HKE')
                energy = self.HKE;
            elseif strcmp(energyTerm, 'VKE')
                energy = self.VKE;
            else
                energy = self.PE;
            end
           
            %%%            

            % interp the matriz [nz, nK] for the same position on the
            % vertical (z)
            if options.mask ==1
                energyMask=energy(:, 1:self.cutoff_modes, 1:self.cutoff_k);
                energyMaskAtzK=squeeze(sum(energyMask,2));

                for i= 1:length(self.KRadialLog(1:self.cutoff_k))                                    
                    HKEatk(i)= interp1(self.zPerMode(:,i),energyMaskAtzK(:,i),z);  
                end

            else
                % energyinterp = scatteredInterpolation(self, energy, -self.Lz:0, self.KRadialLog, 1:self.nModes);
                % energyAtzK=squeeze(sum(energyinterp,2));
                % wvtHKEatZ=interp1(-self.Lz:0,energyAtzK,z);
                energyAtzK=squeeze(sum(energy,2));

                for i= 1:length(self.KRadialLog)                                  
                    HKEatk(i)= interp1(self.zPerMode(:,i),energyAtzK(:,i),z);   
                end                  
            end

            if options.mask==1
                dKLog = self.dKLog;
                HKEatk=HKEatk./dKLog(1:self.cutoff_k);
            else
                HKEatk=HKEatk./self.dKLog;
            end

            %interp on the KRadial vector specified by the user
            if  ~isempty(options.KRadial)
                if options.mask == 1
                   energyAtHorizontalWavenumber= interp1(self.KRadialLog(1:self.cutoff_k),HKEatk,options.KRadial);
                else
                   energyAtHorizontalWavenumber= interp1(self.KRadialLog,HKEatk,options.KRadial); 
                end
            else
                energyAtHorizontalWavenumber=squeeze(HKEatk);
            end
                       

        end



        %%%%%%%%%%%%%%%

        function energyAtVerticalMode = energyAtVerticalMode(self, z, options)
                        % 
            arguments
                self
                z (1,1) double % define double or vector
                options.modeVector double = (1:self.nModes)  
                %options.zVector = linspace(-self.Lz,0,4000) 
                options.plot logical = true
            end

            %Befome summing over K I need to interpolate to the same z
            %Vector. This loop is definitely not the best way of doing it            

            for i= 1:length(self.KRadialLinear)
                for j=1:self.nModes
                    HKESamez(:,j,i)=interp1(self.zPerMode(:,i),self.HKE(:,j,i),z);  
                end
            end

            %Sum over K
            HKEatzMode = squeeze(sum(HKESamez,3));

            %interp at desired depth

            %for i= 1:self.nModes
            %    HKEatMode(i)=interp1(options.zVector,HKEatzMode(:,i),z);             
            %end
            
            energyAtVerticalMode = HKEatzMode;

            if options.plot ==1

                figure(40)

                semilogy(options.modeVector, HKEAtVerticalMode*100)    
                title("HKE")
                ylabel("Variance [cm^2/s^2]")
                xlabel("Vertical Mode")
                grid on
            else
            end     
       
        end
      



        %%%%%%%%%%%%%%%

        function [S,EnergyFrequency] = energyAtFrequencies(self,z,energyTerm,options)
            arguments
                self
                z
                energyTerm %options are: 'TE','HKE','VKE' and 'PE'                
                options.omegaVector = linspace(self.f0,0.8*sqrt(self.N2max),self.nK);
                options.spectrumType
                options.plot = false
            end

            if strcmp(energyTerm, 'TE')
                data = self.HKE + self.VKE +self.PE;
                energy = squeeze(scatteredInterpolation(self, data, z, self.KRadialLog, 1:self.nModes)); 
                
            elseif strcmp(energyTerm, 'HKE')
                data = self.HKE;
                energy = squeeze(scatteredInterpolation(self, data, z, self.KRadialLog, 1:self.nModes)); 
            
            elseif strcmp(energyTerm, 'VKE')
                data = self.VKE;
                energy = squeeze(scatteredInterpolation(self, data, z, self.KRadialLog, 1:self.nModes)); 
            elseif strcmp(energyTerm, 'PE')
                data = self.PE;
                energy = squeeze(scatteredInterpolation(self, data, z, self.KRadialLog, 1:self.nModes)); 
            else
                disp('This option does not exist')
            end
                        

            %Defining omegaVector
            omegaj1=self.omega(1,:);
            dOmega=max(diff(sort(omegaj1(:))));            
            omegaVector=min(self.omega(:)):2*dOmega:max(self.omega(:));

            EnergyFrequency=zeros(length(self.j),length(omegaVector));

            % Redistributing energy over frequency
            for indj=(1:self.nModes) 
                
                for i=(1:length(omegaVector)-1)    
                                  
                % find all the kl point btw the two values of Kh
                    indForOmega = self.omega(indj,:)>=omegaVector(i) & self.omega(indj,:)<omegaVector(i+1);

                    %isRepresented(indj,i)= sum(indForOmega);
                    EnergyFrequency(indj,i) = EnergyFrequency(indj,i) + sum(squeeze(energy(indj,indForOmega)));
                    
                              
                end  
            end



            S=sum(EnergyFrequency,1);          


            
            if options.plot ==1

                figure(60)

                loglog((omegaVector)*(24*3600)/(2*pi), S*100)    
                title(energyTerm)
                ylabel("Variance [cm^2/s^2]")
                xlabel("Frequency [cycle/day]")
                grid on
            else
            end     
              
            
        end


       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Interpolation        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function SIModes = initScatteredInterpolant(self,data,KRadialLog,nModes,options)
             arguments                
                 self 
                 data
                 KRadialLog 
                 nModes
                 options.mask logical = false
             end

                if options.mask 
                     % ZLog
                     ZVectorLog=reshape(self.zPerMode(:,1:self.cutoff_k),[],1);
                else
                     % ZLog
                     ZVectorLog=reshape(self.zPerMode,[],1);
                end

                % KLog
                KVectorRepLog = reshape(repmat(KRadialLog,[self.nZ 1]),[],1);
                lambdaVectorLog= (2*pi)./KVectorRepLog;

             if size(data,3)>1
                 for n=nModes
                     dataVector = reshape(data(:,n,:),[],1);
    
                     SIModes{n}=scatteredInterpolant(ZVectorLog,lambdaVectorLog,dataVector);                          
                 end
             else
                 dataVector = reshape(data,[],1);
    
                 SIModes=scatteredInterpolant(ZVectorLog,lambdaVectorLog,dataVector); 
             end

         end           

         
   

        %%%%%%%%%
        function DataInterpMat = scatteredInterpolation(self, data, zVectorNew, KVectorNew, verticalMode, options)
             arguments                
                 self 
                 data
                 zVectorNew
                 KVectorNew
                 verticalMode
                 options.mask logical = false
             end

             SIModes = initScatteredInterpolant(self,data,KVectorNew, verticalMode,'mask',options.mask);

            
             %zLin
             ZVectorLin=reshape(repmat(zVectorNew,[1, length(KVectorNew)]),[],1);
            
             %Klin    
             lengthZ=length(zVectorNew);
             %VectorRepNew =  reshape(repmat(self.KRadialLog,[self.nZ 1]),[],1);
             KVectorRepNew = reshape(repmat(KVectorNew,[lengthZ 1]),[],1);
             lambdaVectorNew= (2*pi)./KVectorRepNew;

             if size(data,3)>1

                 for i=1:length(verticalMode)
                    
                     DataInterp = SIModes{verticalMode(i)}(ZVectorLin, lambdaVectorNew);
    
                     DataInterpMat1(:,:,i) =reshape(DataInterp, length(zVectorNew),length(KVectorNew));
                 end

             else

                 for i=1:length(verticalMode)

                     DataInterp = SIModes(ZVectorLin, lambdaVectorNew);

                     DataInterpMat1(:,:,i) =reshape(DataInterp, length(zVectorNew),length(KVectorNew));
                 end

             end
             DataInterpMat = permute(DataInterpMat1,[1,3,2]);

          end

       %%%%%%%%%%%%
       function DataInterp2D = interp2D(self, data, KVectorNew, verticalMode)
             arguments                
                 self 
                 data
                 KVectorNew
                 verticalMode
             end
             [X,Y] = ndgrid((1:self.nModes),self.KRadialLog);
             [Xq,Yq]= ndgrid(verticalMode,KVectorNew);

             DataInterp2D = interpn(X,Y,data,Xq,Yq);
           end
             


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Other Usefull Tools        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            function plotStratifcationHighMode(self)
                        % 
            arguments
                self                            
            end

            figure(30)
            subplot(1,3,1)
            plot(sqrt(self.N2zInitial)*3600/(2*pi),self.zInitial ,'k',LineWidth=1.5)
            hold on
            xline(0.8*sqrt(self.N2max)*3600/(2*pi))
            ylabel('depth')
            xlabel('cph')
            title('N(z)')

            subplot(1,3,2)
            plot(self.FInitial(:,end),self.zInitial ,'k',LineWidth=1.5)                      
            title('Initial F - Highest Mode')   

            subplot(1,3,3)
            plot(self.GInitial(:,end),self.zInitial ,'k',LineWidth=1.5)                        
            title('Initial G - Highest Mode')   
            end


        

        %%%%%%%%%%%%%%%%
        function plotQuadraturePoints(self,Mode)
            arguments
                self 
                Mode (1,1) integral
            end
        end


        %%%%%%%%%%%%%%%%
        function dKLog = dKLog(self)
            for iK = 1:self.nK
                    if iK == 1
                        % First Interval
                        lowerBound = 0;  % Use um valor pequeno para evitar zero
                        upperBound = self.KRadialLog(iK)/2;
                    elseif iK == length(self.KRadialLog)
                        % Last Interval
                        lowerBound = upperBound;
                        upperBound = self.KRadialLog(iK);
                    else
                        % Others Intervals
                        lowerBound = upperBound;
                        upperBound = self.KRadialLog(iK) + (self.KRadialLog(iK + 1) - self.KRadialLog(iK) )/2;
                    end
                 dKLog(iK) = (upperBound-lowerBound);        
            end   
        end

        function omegaVector = omegaVector(self)
            omegaj1=self.omega(1,:);
            dOmega=max(diff(sort(omegaj1(:))));            
            omegaVector=min(self.omega(:)):2*dOmega:max(self.omega(:));
        end

    end
end

    

