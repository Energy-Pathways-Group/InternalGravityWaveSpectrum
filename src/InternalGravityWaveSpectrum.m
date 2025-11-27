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

        A2

        HKE, VKE, PE, TE            
        TEPSD

        zNew

        cutoff_modes, cutoff_k

        new_KRadialLog      

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
        %N2max = max(N2zInitial);
        N2max = 1.2474e-05;

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
        minOrder = 1; 
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
        
    end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 4: Computation of energy distribution according 
        % with spectrum S
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function selfUpdated = assignEnergySpectrum(self,S) 
            arguments
                self
                S = self.generalSpectrum(j_star=3,slope_j=1,slope_k=1,A=1)                
            end
            
            
   
            % -------------------------
            % Compute total energy (A2) using spectrum
            % -------------------------
            clear TE
                    
            TE = self.amplitudesWithSpectrum1(S);
                       
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computes total energy per mode (K,j) using integral  % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


        function totalEnergyPerComponent = amplitudesWithSpectrum1(self, spectrum)
            arguments
                self {mustBeNonempty}
                spectrum {mustBeNonempty}                
            end
        
            % Initialize output
            totalEnergyPerComponent = zeros(self.nModes, self.nK);
        
            % Ensure KRadialLog is column vector
            K = self.KRadialLog(:);
            nK = length(K);
        
            % Compute geometric bin edges for log-spaced grid
            edges = zeros(nK+1,1);
            edges(1) = K(1);
            for i = 1:nK-1
                edges(i+1) = sqrt(K(i) * K(i+1));
            end
            edges(end) = K(end);
        
            % Safe spectrum function to avoid out-of-range evaluation
            Ssafe = @(k,jInd) spectrum(min(max(k, K(1)), K(end)), jInd);
        
            % Integrate spectrum over each bin
            for iJ = 1:self.nModes
                jval = self.j(iJ);
                for iK = 1:nK
                    lb = edges(iK);
                    ub = edges(iK+1);
        
                    % Total energy in the bin (TP)
                    TP = integral(@(kk) Ssafe(kk,jval), lb, ub, ...
                                  'RelTol',1e-8, 'AbsTol',1e-12);
                    totalEnergyPerComponent(iJ,iK) = TP;
        
                    % Convert to PSD and                    
                    binWidth = ub - lb;
                    self.TEPSD(iJ,iK) = TP / binWidth;                   
                    totalEnergyPerComponent(iJ,iK) = TP;
                    
                end
            end

            % Diagnostic printout (only makes sense for TP)        
            for iJ = 1:self.nModes
                jval = self.j(iJ);
                fullIntegral = integral(@(k) Ssafe(k,jval), K(1), K(end));
                sumBins = sum(totalEnergyPerComponent(iJ,:));
                fprintf('Mode %d: fullIntegral=%.6g, sumBins=%.6g, rel diff=%.2e\n', ...
                        iJ, fullIntegral, sumBins, ...
                        abs(fullIntegral - sumBins)/fullIntegral);
            end
            
        end
                    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Fitting Spectral functions 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function S_normalized = normalizeSpectrum(self,S)
            % Takes a function_handle S with arguments (k,j) and returns
            % a function_handle with the same arguments normalized to GM
            % energy level 1.
            kmin = min(self.KRadialLog(1));
            kmax = max(self.KRadialLog(end));

            S_norm = ones(self.nModes, 1);
            for jIdx = 1:length(self.j)
                S_norm(jIdx) = integral(@(k) S(k, self.j(jIdx)), kmin, kmax);
            end

            L_gm = 1.3e3;         % Thermocline exponential scale [m]
            invT_gm = 5.2e-3;     % Reference buoyancy frequency [rad/s]
            E_gm = 6.3e-5;        % Non-dimensional energy parameter
            E = (L_gm^3) * (invT_gm^2) * E_gm; % Total GM energy

            norm = E/sum(S_norm);
            S_normalized = @(k, jInd) norm*S(k,jInd);
        end

        function S = gmSpectrum(self,p)
            %igw.gmSpectrum(j_star=params(1),slope_j=params(2),slope_k=params(3),A=params(4));
            
            arguments
                self 
                p.j_star = 3;
                
            end

            A = 1;
            slope_j=5/4;
            slope_k=1;    
            % -------------------------
            % Compute Rossby radius of deformation
            % -------------------------
            Lr2_ = (self.g .* self.h) ./ (self.f0 ^ 2);
            Lr2_func = @(k,jInd) interp1(self.KRadialLog, Lr2_(jInd,:), k);
            

            S_unnorm = @(k,jInd) sqrt(Lr2_func(k, jInd)) ./ ( ((k.^2 .* Lr2_func(k, jInd) + 1).^slope_k) .* (p.j_star^2 + self.j(jInd).^2).^(slope_j) );
            S_normalized = self.normalizeSpectrum(S_unnorm);
            S = @(k,jInd) A * S_normalized(k,jInd);
        end


       
        function S = generalSpectrum(self,p)
            arguments
                self 
                p.j_star = 3;
                p.slope_j = 1;
                p.slope_k = 1;
                p.A = 1;
            end                
            
            % -------------------------
            % Compute Rossby radius of deformation
            % -------------------------
            Lr2_ = (self.g .* self.h) ./ (self.f0 ^ 2);
            Lr2_func = @(jInd, k) interp2(self.KRadialLog, self.j, Lr2_, k, jInd);

            % -------------------------
            % Compute k*^2 using mode-3 Rossby radius of deformations 
            % -------------------------
            %jVec = p.j_star * ones(size(self.KRadialLog));  % vetor same size than K
            %jVec = 3 * ones(size(self.KRadialLog));

            %kstar2 = squeeze(1 ./ Lr2_func(jVec, self.KRadialLog));
            %kstar2 = squeeze(1 ./ Lr2_func(self.j, self.KRadialLog));
            %kstar2_func = @(k) interp1(self.KRadialLog, kstar2, k);
            kstar2_func= @(jInd,k) 1./Lr2_func(jInd,k);
        
            % -------------------------
            % Define spectral function S(k,j)
            % -------------------------
            S_unnorm = @(k, jInd) sqrt(Lr2_func(jInd,k)) ./ ...
               ( ((k.^2./kstar2_func(jInd,k) + 1).^(p.slope_k)).*...
               ((Lr2_func(p.j_star,k)./Lr2_func(jInd,k)+1).^(p.slope_j)) );

            S_normalized = self.normalizeSpectrum(S_unnorm);
            S = @(k,jInd) p.A * S_normalized(k,jInd);    
               
        end


        function M = tidalSpectrum(self,p)
            
            arguments
                self 
                p.tidalOmega = 2*pi/(12.42*3600);
                p.A = 100
                p.c = 7*10^-6
                p.j0 = 3.5
                p.d=0
                
            end
            omegaFunc= @(k,jInd) interp1(self.KRadialLog, self.omega(jInd,:), k);
            
            
            M = @(k,jInd) (p.A^2 * p.c^2) ./ ((omegaFunc(k,jInd) - p.tidalOmega).^2 + p.c^2) .* ...
                 (1 ./ ((jInd - p.j0).^2 + p.d^2));
        end


        function SfM2 = fM2Spectrum(self,p)
            
            arguments
                self 
                p.fM2Omega = 2*pi/(12.42*3600) + self.f0;
                %p.j_center = 1
                p.A = 100
                p.c = 10
                p.j0 = 2
                p.d=0.5
                
                
                
            end
            omegaFunc= @(k,jInd) interp1(self.KRadialLog, self.omega(jInd,:), k);
            
            
            SfM2 = @(k,jInd) (p.A^2 * p.c^2) ./ ((omegaFunc(k,jInd) - p.fM2Omega).^2 + p.c^2).* ...
                 (1 ./ ((jInd - p.j0).^2 + p.d^2));
            
            
        end

         function Sf = fSpectrum(self,p)
            
            arguments
                self 
                p.A = 100
                p.c = 10
                p.j0 =0
                p.d=3        
                
                
            end
            omegaFunc= @(k,jInd) interp1(self.KRadialLog, self.omega(jInd,:), k);
            
            
            Sf = @(k,jInd) (p.A^2 * p.c^2) ./ ((omegaFunc(k,jInd) - self.f0).^2 + p.c^2).* ...
                 (1 ./ ((jInd - p.j0).^2 + p.d^2));
            

         end

        function val = objectiveFunction(self,S,k,j,S_wvt)
            arguments
                self
                S
                k
                j
                S_wvt               
            end
            energySpectrum = S;
            dk = (k(2)-k(1));

            Ejk = zeros(length(j),length(k));
            for indj = 1:length(j)
                for indk = 1:length(k)
                    Ejk(indj,indk) = energySpectrum(k(indk)+(dk/1.9),j(indj));
                end
            end

           
            if(1)
                % % First pcolor plot
                % subplot(1, 2, 1);
                % jpcolor(log(S_wvt));                     
                % colorbar;
                % clim([2 10])
                % title('Model');                
                
                % Second pcolor plot
                %subplot(1, 2, 2);
                jpcolor(2*pi*k, j, log10(Ejk))
                cb=colorbar;
                clim([1 3.5])
                cb.FontSize = 10;
                
                ax = gca;
                ax.FontSize = 12;
                ax.XLabel.FontSize = 12;
                ax.YLabel.FontSize = 12;
                ax.Title.FontSize = 12;
                %cmocean('haline')
                %title('Toolbox');
            end

            val = mean((log(S_wvt(:)) - log(Ejk(:))).^2);
        end


        function val = SGMTide(self,params, k, j, Ejk)
            model_GM = self.gmSpectrum(j_star=params(1));            
            model_M2 = self.tidalSpectrum(A=params(2),c=params(3));
            model_fM2 = self.fM2Spectrum(A=params(4),c=params(5));
            S = @(k,jInd) model_GM(k,jInd) + model_M2(k,jInd) + model_fM2(k,jInd);
            val = self.objectiveFunction(S,k,j,Ejk);
        end

        function val = SgenTide(self,params, k, j, Ejk)
            model_GM = self.generalSpectrum(j_star=params(1),slope_j=params(2),slope_k=params(3),A=params(4));
            model_M2 = self.tidalSpectrum(A=params(5),c=params(6));
            model_fM2 = self.fM2Spectrum(A=params(7),c=params(8));
            S = @(k,jInd) model_GM(k,jInd) + model_M2(k,jInd) + model_fM2(k,jInd);
            val = self.objectiveFunction(S,k,j,Ejk);
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
                %options.omega2D =self.omega
                options.omegaVector %= linspace(self.f0,0.8*sqrt(self.N2max),self.nK);
                %options.new_KRadialLog = self.KRadialLog
                options.spectrumType
                options.plot = false
                options.mask = true
            end
    
            % new_wavelengthLog=logspace(2,6,1000);     
            % new_KRadialLog=fliplr((2*pi)./new_wavelengthLog);
            % self.new_KRadialLog = new_KRadialLog;
            % self.new_KRadialLog(1)=0.5*self.new_KRadialLog(1);
           
            
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
            
            % omegaj1=self.omega(1,:);
            % dOmega=max(diff(sort(omegaj1(:))));            
            % omegaVector=min(self.omega(:)):2*dOmega:max(self.omega(:));

            EnergyFrequency=zeros(length(self.j),length(options.omegaVector));

            % Redistributing energy over frequency
            for indj=(1:self.nModes) 
                
                for i=(1:length(options.omegaVector)-1)    
                                  
                % find all the kl point btw the two values of Kh
                    indForOmega = self.omega(indj,:)>=options.omegaVector(i) & self.omega(indj,:)<options.omegaVector(i+1);

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

    

