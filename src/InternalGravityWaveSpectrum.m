classdef InternalGravityWaveSpectrum < handle
% INTERNALGRAVITYWAVESPECTRUM  Spectral model for internal gravity waves
%
%   This class constructs a spectral representation of internal gravity
%   waves in a stratified ocean using vertical modes and a logarithmically
%   spaced horizontal wavenumber grid.
%
%   The model:
%     - Computes vertical modes using InternalModesSpectral
%     - Evaluates modal structure functions (F, G), dispersion relation (ω),
%       and modal amplitudes (h)
%     - Constructs energy coefficients for horizontal kinetic energy (HKE),
%       vertical kinetic energy (VKE), and potential energy (PE)
%
%   The resulting fields are defined on a nonuniform vertical grid
%   (Gauss quadrature points per mode and wavenumber), which motivates the
%   use of scattered interpolation methods throughout the toolbox.
%
%   COORDINATE CONVENTIONS
%     - z < 0 below the surface, z = 0 at the surface
%     - Total horizontal wavenumber K is radial and log-spaced
%     - Frequencies ω satisfy f0 < ω < sqrt(N2max)
%
%   LIMITATIONS
%     - Not valid near the equator (|latitude| < 5°)
%     - Currently assumes scalar vertical interpolation queries
%
%   See also:
%     InternalModesSpectral
%     InternalModesWKBSpectral



    properties (Access = public)
        latitude 
        f0 
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
        omega % size(omega_k) = [nModes,nK]       
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

    properties %(Access = private, Hidden)
        %Stored for testing / debugging 
        FInitial % [nZ,nModes]
        GInitial % [nZ,nModes]
        zInitial % [nZ]
        N2zInitial
        Lr2Interpolant   % griddedInterpolant object
               
    end
   

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization/ Constructor
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function self = InternalGravityWaveSpectrum(N2,Lz,options)                    
        % INTERNALGRAVITYWAVESPECTRUM
        %  
        %   Constructs an internal gravity wave spectral model for a given
        %   stratification profile and ocean depth. It's a initialization.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %
        %   N2      : function handle, buoyancy frequency squared N²(z)
        %   Lz      : scalar, total ocean depth [helpful: positive number]
        %
        %   options : name-value pairs
        %       .latitude                    : latitude in degrees (default = 33)
        %       .nModes                      : number of vertical modes (default = 64)
        %       .nK                          : number of horizontal wavenumbers (default = 64)
        %       .shouldForceMonotonicDensity : enforce monotonic N²(z) (default = false)
        %
        % -------------------------------------------------------------------------
        % OUTPUT
        %
        %   self : initialized InternalGravityWaveSpectrum object
        %
        % -------------------------------------------------------------------------
        % NOTES
        %
        %   - The model is not valid near the equator (|latitude| < 5°)
        %   - Vertical modes are computed using Gauss quadrature points that
        %     vary with horizontal wavenumber
        %   - Energy coefficients follow the squared linear wave equations
        %     (see Jeffrey et al. 2021)
        %
        %--------------------------------------------------------------------------

            arguments
                N2 function_handle
                Lz (1,1) double {mustBePositive}
                options.latitude (1,1) double = 33
                options.nModes (1,1) double = 64
                options.nK (1,1) double = 64
                options.shouldForceMonotonicDensity (1,1) logical = false
            end

                    
            % ==================================================
            % Latitude validity checks
            % ==================================================
            if abs(options.latitude) <= 5
                error("Latitude:MustBeAwayEquator", ...
                    "Model is not valid within 5° of the equator.")
            elseif abs(options.latitude) >= 85
                error("Latitude:InvalidValue", ...
                    "Latitude must be between ±85°.")
            end          
           
            % ==================================================
            % Store basic parameters
            % ==================================================
            self.N2      = N2;
            self.Lz      = Lz;
            self.latitude = options.latitude;
            self.nModes  = options.nModes;
            self.nK      = options.nK;
            self.nZ      = options.nModes + 1;
            self.g       = 9.80665;
            self.shouldForceMonotonicDensity = options.shouldForceMonotonicDensity;
            self.j       = 1:self.nModes;

            self.cutoff_modes = ceil(2/3 * self.nModes);
            self.cutoff_k     = ceil(2/3 * self.nK);
                

            % ==================================================
            % Stratification and Coriolis frequency
            % ==================================================
            omegaEarth = 7.2921e-5;
            self.f0 = 2 * omegaEarth * sind(self.latitude);

            zInitial = linspace(-Lz, 0, 10001);
            N2zInitial = self.N2(zInitial);
            self.N2max = max(N2zInitial);

            %Value for tests
            %N2max = 1.2474e-05;
       
            if self.shouldForceMonotonicDensity
                validateattributes(N2zInitial, {'numeric'}, {'increasing'})
            end

            self.zInitial  = zInitial;
            self.N2zInitial = N2zInitial;

            % ==================================================
            % Estimate maximum horizontal wavenumber
            % ==================================================
            imWKB = InternalModesWKBSpectral( ...
                N2=self.N2, zIn=[-Lz 0], zOut=zInitial, ...
                latitude=self.latitude, nModes=self.nModes);

            [F0, G0, ~, k] = imWKB.ModesAtFrequency(0.95 * sqrt(self.N2max));
            Kmax = max(k);

            self.FInitial = F0;
            self.GInitial = G0;
            

            % ==================================================
            % Log-spaced horizontal wavenumber grid
            % ==================================================
            minOrder = max(2, floor(log10(2*pi / Kmax)));

            %Value for testing
            %minOrder=2;
        
            wavelengthLog = logspace(minOrder, 6, self.nK);
            self.KRadialLog = fliplr((2*pi) ./ wavelengthLog);
            self.KRadialLog(1) = 0.5 * self.KRadialLog(1);           
           

            % ==================================================
            % Vertical modes, eigendepths, and dispersion relation
            % ==================================================        
            upperBoundary = UpperBoundary.rigidLid;
            normalization = Normalization.kConstant;    
            

            for iK = 1:self.nK
                % First pass: obtain quadrature points
                im = InternalModesSpectral( ...
                    N2=self.N2, zIn=[-Lz 0], zOut=zInitial, ...
                    latitude=self.latitude, nModes=self.nModes);

                im.upperBoundary = upperBoundary;
                im.normalization = normalization;  
    
                zPerModeLog(:,iK) = im.GaussQuadraturePointsForModesAtWavenumber(...
                    self.nModes+1,self.KRadialLog(iK));
    
                % Second call: recompute modes on quadrature grid for efficiency
                im = InternalModesSpectral(...
                    N2=self.N2,zIn=[-Lz 0],zOut=zPerModeLog(:,iK),...
                    latitude=self.latitude,nModes=self.nModes);

                [F(:,:,iK),G(:,:,iK),h(:,iK),omega(:,iK)] = ...
                    im.ModesAtWavenumber(self.KRadialLog(iK));
    
            end
 
            self.zPerMode =zPerModeLog;
            self.F = F;
            self.G = G;
            self.h = h;
            self.omega = omega;  

            % ========================================================
            % Rossby radius squared on native grid
            % Size: [nModes x nK]
            % ========================================================        
           
            Lr2_ = (self.g .* self.h) ./ (self.f0^2);

            % IMPORTANT:
            % griddedInterpolant expects NDGRID ordering
            % Here: (j, k)
            self.Lr2Interpolant = griddedInterpolant( ...
                {self.j, self.KRadialLog}, ...
                Lr2_, ...
                'linear', ...   % interpolation
                'nearest');     % extrapolation

            self.Lr2Interpolant = @(j,k) reshape( interp1(self.KRadialLog,(self.g .* self.h(j,:)) ./ (self.f0^2),k,"linear"), size(k));
            

            % ========================================================
            % Energy coefficients from linear wave theory
            % ========================================================
        
            self.HKEcoef = 0.25 * (1+ (self.f0^2./(self.omega.^2))); 
            self.VKEcoef = 0.25 * ((self.KRadialLog).^2 .* self.h.^2);        
            self.PEcoef  = 0.25 * (((self.KRadialLog).^2.*self.h.^2)./self.omega.^2);          
        end

        function y = interp_at(self,d,j,k)
            FF = d{j};   % griddedInterpolant
            y = FF(k);   % k can be a vector
        end

        function selfUpdated = assignEnergySpectrum(self,S) 

        % ASSIGNENERGYSPECTRUM  Compute and assign wave–energy components from a spectrum
        %
        %   self = ASSIGNENERGYSPECTRUM(self, S) applies a user-defined spectrum function
        %   handle S to compute the total wave energy and its decomposition into:
        %       - Horizontal Kinetic Energy (HKE)
        %       - Vertical   Kinetic Energy (VKE)
        %       - Potential  Energy (PE)
        %
        %   The method stores all results back into the object.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %
        %   self  : Model object
        %
        %   S     : Function handle defining the energy spectrum.
        %           It must have the signature:   A = S(j, k)
        %
        %           If omitted, the default spectrum is:
        %               S = self.generalSpectrum(j_star=3, slope_j=1, slope_k=1, A=1)
        %
        % -------------------------------------------------------------------------
        % OUTPUT
        %
        %   selfUpdated : Updated model object containing:
        %       self.A2
        %       self.TE
        %       self.HKE
        %       self.VKE
        %       self.PE
        %       self.N2atQuadPoints
        %
        % -------------------------------------------------------------------------
        % EXAMPLE
        %
        %   S = myModel.generalSpectrum(j_star=4, slope_j=1.2, slope_k=1.1, A=2);
        %   myModel = myModel.assignEnergySpectrum(S);
        %
        % -------------------------------------------------------------------------

            arguments
                self
                S (1,1) function_handle = self.generalSpectrum( ...
                        j_star=3, slope_j=1, slope_k=1, A=1)
            end
            

            % ================================================================
            % Validate spectrum function handle with a test call
            % ================================================================
            try
                % Test using minimal valid input shapes 
                S(1e-4,1);
            catch ME
                error("Invalid spectrum function handle S. " + ...
                      "It must be callable as S(j, k).\nOriginal error:\n%s", ME.message);
            end
   
            % ================================================================
            % Compute total energy from the spectrum
            % ================================================================
            self.TE = self.amplitudesWithSpectrum(S,false);

            % Normalize by eigendepth
            self.A2 = 2 * self.TE ./ self.h;

            % ================================================================
            % Precompute stratification at modal quadrature points
            % ================================================================
            self.N2atQuadPoints = self.N2(self.zPerMode);

            % ================================================================
            % Energy components
            % ================================================================
            self.HKE = shiftdim(self.A2 .* self.HKEcoef, -1) .* self.F.^2;
            self.VKE = shiftdim(self.A2 .* self.VKEcoef, -1) .* self.G.^2;
        
            % PE requires correct handling of N² shape
            if isscalar(self.N2atQuadPoints)
                % Constant stratification
                N2local = self.N2atQuadPoints;
            else
                % Reshape to match dimensions [nZ 1 nK]
                N2local = reshape(self.N2atQuadPoints, [self.nZ 1 self.nK]);
            end
        
            self.PE = shiftdim(self.A2 .* self.PEcoef, -1) .* self.G.^2 .* N2local;
            selfUpdated = self;
    
        end



        function totalEnergyPerComponent = amplitudesWithSpectrum(self, spectrum, verbose)

        % AMPLITUDESWITHSPECTRUM  Compute total energy from a spectrum
        %
        %   totalEnergyPerComponent = AMPLITUDESWITHSPECTRUM(self, spectrum)
        %   computes the total energy for each vertical mode and horizontal wavenumber
        %   bin given a user-defined spectrum function handle.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %
        %   self     : Model object containing:
        %              - nModes : number of vertical modes
        %              - nK     : number of horizontal wavenumber bins
        %              - KRadialLog : radial wavenumber vector (log-spaced)
        %              - j      : vertical mode numbers
        %
        %   spectrum : Function handle defining the energy spectrum.
        %              Must have the signature:  S = spectrum(k, j)
        %
        %   verbose  : Logical flag (true/false). If true, prints diagnostic info.
        %              Default: false
        %
        % -------------------------------------------------------------------------
        % OUTPUTS
        %
        %   totalEnergyPerComponent : [nModes x nK] array containing the total energy
        %                             in each vertical mode and horizontal bin.
        %
        % -------------------------------------------------------------------------
        % EXAMPLES
        %
        %   S = @(k,j) exp(-k.^2) .* j; 
        %   energy = myModel.amplitudesWithSpectrum(S, true);
        %
        % -------------------------------------------------------------------------    

            arguments
                self {mustBeNonempty}
                spectrum {mustBeNonempty, mustBeA(spectrum,'function_handle')}
                verbose logical = false
            end
    
       
            % ========================================================
            % Initialize output
            % ========================================================
            totalEnergyPerComponent = zeros(self.nModes, self.nK);
        
            % Ensure KRadialLog is column vector
            K = self.KRadialLog(:);
            localnK = self.nK;
        
            % ========================================================
            % Geometric bin edges for log-spaced K
            % ========================================================
            edges = [K(1); sqrt(K(1:end-1).*K(2:end)); K(end)];
        
            % ========================================================
            % Safe spectrum function to avoid out-of-range evaluation
            % ========================================================
            Ssafe = @(k,jInd) spectrum(min(max(k, K(1)), K(end)), jInd);            
        
            % ========================================================
            % Integrate spectrum over each bin
            % ========================================================
            for iJ = 1:self.nModes
                jval = self.j(iJ);
                for iK = 1:localnK
                    lb = edges(iK);
                    ub = edges(iK+1);
        
                    % Total energy in the bin (TP)
                    TP = integral(@(kk) Ssafe(kk,jval), lb, ub, ...
                                  'RelTol',1e-8, 'AbsTol',1e-12);
                    totalEnergyPerComponent(iJ,iK) = TP;
        
                    % Store PSD per bin
                    self.TEPSD(iJ,iK) = TP / (ub - lb);          
                                      
                end
            end

            % ========================================================
            % Optional diagnostic printout
            % ========================================================    
            if verbose
                for iJ = 1:self.nModes
                    jval = self.j(iJ);
                    fullIntegral = integral(@(k) Ssafe(k,jval), K(1), K(end));
                    sumBins = sum(totalEnergyPerComponent(iJ,:));
                    fprintf('Mode %d: fullIntegral=%.6g, sumBins=%.6g, rel diff=%.2e\n', ...
                            iJ, fullIntegral, sumBins, abs(fullIntegral - sumBins)/fullIntegral);
                end
            end            
            
        end

                    
        function testTEequalsSumOfComponents(self)
        % TESTTEEQUALSSUMOFCOMPONENTS  Verify TE consistency via error diagnostics
        %
        %   This test checks whether the total energy (TE) computed directly from
        %   the spectrum is consistent with the sum of its components:
        %
        %       TE ≈ HKE + VKE + PE
        %
        %   Procedure:
        %     1. Interpolate HKE, VKE, and PE onto a common vertical grid
        %     2. Integrate each component over depth
        %     3. Apply the standard high-mode / high-k cutoff
        %     4. Compute absolute and relative error fields
        %
        %   Output:
        %     - Two-panel figure:
        %         (left)  log10 absolute error
        %         (right) relative error in percent (%)
        %     - Printed scalar absolute and relative errors (L2 norms)
        %
        
            arguments
                self {mustBeNonempty}
            end
        
            % ========================================================
            % Define common vertical grid
            % ========================================================
            zTest = linspace(-self.Lz, 0, 1000);
        
            % ========================================================
            % Interpolate energy components onto common z-grid
            % Output size: [nz x nModes x nK]
            % ========================================================
            HKEz = self.scatteredInterpolation(self.HKE, zTest, self.KRadialLog, self.j);
            VKEz = self.scatteredInterpolation(self.VKE, zTest, self.KRadialLog, self.j);
            PEz  = self.scatteredInterpolation(self.PE,  zTest, self.KRadialLog, self.j);
        
            % ========================================================
            % Integrate over depth (dimension 1 = z)
            % Result size: [nModes x nK]
            % ========================================================
            HKEint = squeeze(trapz(zTest, HKEz, 1));
            VKEint = squeeze(trapz(zTest, VKEz, 1));
            PEint  = squeeze(trapz(zTest, PEz,  1));
        
            % ========================================================
            % Apply standard cutoff (remove highest 1/3 modes and k)
            % ========================================================
            self.cutoff_modes = ceil(self.nModes * 2/3);
            self.cutoff_k     = ceil(self.nK     * 2/3);
        
            HKEint = HKEint(1:self.cutoff_modes, 1:self.cutoff_k);
            VKEint = VKEint(1:self.cutoff_modes, 1:self.cutoff_k);
            PEint  = PEint( 1:self.cutoff_modes, 1:self.cutoff_k);
        
            TE_components = HKEint + VKEint + PEint;
            TE_direct     = self.TE(1:self.cutoff_modes, 1:self.cutoff_k);
        
            % ========================================================
            % Error fields
            % ========================================================
            absErrField = abs(TE_components - TE_direct);
        
            % Relative error in %
            relErrField = 100 * absErrField ./ max(abs(TE_direct), eps);
        
            % Log-scaled absolute error (avoid log(0))
            logAbsErrField = log10(max(absErrField, eps));
        
            % ========================================================
            % Scalar error metrics (single values)
            % ========================================================
            absErr = norm(absErrField(:), 2);
            relErr = absErr / norm(TE_direct(:), 2);
        
            fprintf('Energy consistency check:\n');
            fprintf('  Absolute L2 error : %.3e\n', absErr);
            fprintf('  Relative L2 error : %.2f %%\n', 100*relErr);
        
            % ========================================================
            % Visualization: error diagnostics
            % ========================================================
            figure('Color','w');
        
            % --- Panel 1: log10 absolute error ---
            subplot(1,2,1)
            pcolor(self.KRadialLog(1:self.cutoff_k), ...
                   1:self.cutoff_modes, ...
                   AbsErrField);
            shading flat
            colorbar
            xlabel('k [rad/m]')
            ylabel('Vertical mode j')
            title('log_{10} |TE_{components} - TE_{direct}|')
        
            % --- Panel 2: relative error (%) ---
            subplot(1,2,2)
            pcolor(self.KRadialLog(1:self.cutoff_k), ...
                   1:self.cutoff_modes, ...
                   relErrField);
            shading flat
            colorbar
            xlabel('k [rad/m]')
            ylabel('Vertical mode j')
            title('Relative error (%)')
        
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Fitting Spectral functions 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function S_normalized = normalizeSpectrum(self,S_in)

        % NORMALIZESPECTRUM  Normalize a spectrum function handle to GM energy level 1
        %
        %   S_normalized = NORMALIZESPECTRUM(self, S) takes a function handle S(k,j)
        %   representing the energy spectrum and returns a new function handle
        %   normalized such that its total energy matches the Garrett–Munk (GM)
        %   reference energy level.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %
        %   self : Model object containing:
        %          - nModes     : number of vertical modes
        %          - KRadialLog : radial wavenumber vector (log-spaced)
        %          - j          : vector of vertical mode numbers
        %
        %   S    : Function handle of the energy spectrum, signature S(k,j)
        %
        % -------------------------------------------------------------------------
        % OUTPUTS
        %
        %   S_normalized : Function handle of the normalized spectrum
        %
        % -------------------------------------------------------------------------
        % EXAMPLE
        %
        %   S = @(k,j) exp(-k.^2).*j;
        %   S_norm = myModel.normalizeSpectrum(S);
        %
        % -------------------------------------------------------------------------

            arguments
                self {mustBeNonempty}
                S_in {mustBeNonempty, mustBeA(S_in,'function_handle')}
            end

            % ========================================================
            % Determine integration bounds
            % ========================================================
            kmin = self.KRadialLog(1);
            kmax = self.KRadialLog(end);

            % ========================================================
            % Integrate spectrum over whole wavenumber range for each vertical mode
            % ========================================================
            S_norm = ones(self.nModes,1);
    
            for jIdx = 1:self.nModes
                S_norm(jIdx) = integral(@(k) S_in(k, jIdx), kmin, kmax, ...
                              ...
                              'RelTol',1e-8, 'AbsTol',1e-12);
            end

            % ========================================================
            % Garrett–Munk reference energy
            % ========================================================
            L_gm = 1.3e3;         % Thermocline exponential scale [m]
            invT_gm = 5.2e-3;     % Reference buoyancy frequency [rad/s]
            E_gm = 6.3e-5;        % Non-dimensional energy parameter
            E = (L_gm^3) * (invT_gm^2) * E_gm; % Total GM energy

            % ========================================================
            % Compute normalization factor and return normalized spectrum
            % ========================================================
            normFactor = E / sum(S_norm);
            S_normalized = @(k, jInd) normFactor * S_in(k, jInd);
        end

        function S = gmSpectrum(self, p)
        % GMSPECTRUM  Generate a Garrett–Munk type internal wave spectrum
        %
        %   S = GMSPECTRUM(self, p) returns a function handle S(k,j) representing
        %   a Garrett–Munk (GM) type internal wave spectrum for vertical modes and
        %   horizontal wavenumbers.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %   self : Model object containing relevant physical parameters
        %   p    : Structure with optional parameters:
        %          - j_star : GM reference vertical mode (default: 3)
        %
        % OUTPUTS
        %   S : Function handle of the GM spectrum normalized to GM energy level 1.
        %       Signature: S(k,j)
        %
        % -------------------------------------------------------------------------
        
            arguments
                self
                p.j_star = 3;
            end
        
            % ========================================================
            % Fixed GM parameters
            % ========================================================
            A = 1;
            slope_j = 5/4;
            slope_k = 1;
        
            % ========================================================
            % Compute Rossby radius of deformation squared per mode
            % ========================================================
            Lr2_ = (self.g .* self.h) ./ (self.f0 ^ 2);
            Lr2_func = @(k,jInd) interp1(self.KRadialLog, Lr2_(jInd,:), k, 'linear');
             
            % ========================================================
            % Unnormalized GM spectrum
            % ========================================================
            S_unnorm = @(k,jInd) sqrt(Lr2_func(k, jInd)) ./ ...
                ( ((k.^2 .* Lr2_func(k, jInd) + 1).^slope_k) .* ...
                  (p.j_star^2 + self.j(jInd).^2).^slope_j );

            % ========================================================
            % Normalize spectrum to GM energy level 1
            % ========================================================
            S_normalized = self.normalizeSpectrum(S_unnorm);
        
            % ========================================================
            % Apply amplitude scaling
            % ========================================================
            S = @(k,jInd) A * S_normalized(k,jInd);
        
        end

       

        function S = generalSpectrum(self, p)
        % GENERALSPECTRUM  Generate a tunable internal wave energy spectrum
        %
        %   S = GENERALSPECTRUM(self, p) returns a function handle S(k,j) defining
        %   a generalized internal wave energy spectrum with user-specified
        %   slopes, reference mode, and amplitude.
        %
        %   The spectrum is normalized to GM energy level 1 before applying
        %   the amplitude scaling factor p.A.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %
        %   self : Model object containing:
        %          - g          : gravitational acceleration
        %          - h          : eigendepths or mode depths
        %          - f0         : Coriolis parameter
        %          - j          : vertical mode numbers
        %          - KRadialLog : radial wavenumber vector (log-spaced)
        %
        %   p    : Structure with parameters:
        %          - j_star  : reference vertical mode (default: 3)
        %          - slope_j : vertical mode slope      (default: 1)
        %          - slope_k : horizontal wavenumber slope (default: 1)
        %          - A       : amplitude scaling factor (default: 1)
        %
        % -------------------------------------------------------------------------
        % OUTPUTS
        %
        %   S : Function handle of the generalized spectrum.
        %       Signature: S(k,j)
        %
        % -------------------------------------------------------------------------
        
            arguments
                self
                p.j_star  = 3;
                p.slope_j = 1;
                p.slope_k = 1;
                p.A       = 1;
            end

            assert(~isempty(self.Lr2Interpolant), ...
            'Lr2Interpolant not initialized.');

        
            % ========================================================
            % Rossby radius of deformation squared
            % ========================================================

            % Old, simple and performance killer method
            % Lr2_ = (self.g .* self.h) ./ (self.f0 ^ 2);
            % 
            % Lr2_func = @(k,jInd) interp2(self.KRadialLog, self.j, Lr2_, ...
            %                              k, jInd, 'linear');

            % New method (Jan 6th, 2026) 
            % ========================================================
            % Rossby radius interpolant handles
            % ========================================================
            % Safe wrapper (vectorized in k)
            
            % self.Lr2 = @(k,j) self.Lr2Interpolant(j .* ones(size(k)), k );
            % self.Lr2  = @(k,j) self.Lr2Interpolant(j, k);
            kstar2 = @(k,j) 1 ./ self.Lr2(k,j);

            % ========================================================
            % Unnormalized generalized spectrum
            % ========================================================
            S_unnorm = @(k,jInd) sqrt(self.Lr2(k,jInd)) ./ ...
                ( (k.^2 ./ kstar2(k,jInd) + 1).^p.slope_k .* ...
                  (self.Lr2(k,p.j_star) ./ self.Lr2(k,jInd) + 1).^p.slope_j );
        
            % ========================================================
            % Normalize and apply amplitude scaling
            % ========================================================
            S_normalized = self.normalizeSpectrum(S_unnorm);
            S = @(k,jInd) p.A * S_normalized(k,jInd);
        
        end

        function S = generalSpectrumKLr(self, p)
        % GENERALSPECTRUM  Generate a tunable internal wave energy spectrum
        %
        %   S = GENERALSPECTRUM(self, p) returns a function handle S(k,j) defining
        %   a generalized internal wave energy spectrum with user-specified
        %   slopes, reference mode, and amplitude.
        %
        %   The spectrum is normalized to GM energy level 1 before applying
        %   the amplitude scaling factor p.A.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %
        %   self : Model object containing:
        %          - g          : gravitational acceleration
        %          - h          : eigendepths or mode depths
        %          - f0         : Coriolis parameter
        %          - j          : vertical mode numbers
        %          - KRadialLog : radial wavenumber vector (log-spaced)
        %
        %   p    : Structure with parameters:
        %          - j_star  : reference vertical mode (default: 3)
        %          - slope_j : vertical mode slope      (default: 1)
        %          - slope_k : horizontal wavenumber slope (default: 1)
        %          - A       : amplitude scaling factor (default: 1)
        %
        % -------------------------------------------------------------------------
        % OUTPUTS
        %
        %   S : Function handle of the generalized spectrum.
        %       Signature: S(k,j)
        %
        % -------------------------------------------------------------------------
        
            arguments
                self
                p.j_star  = 3;
                p.slope_j = 1;
                p.slope_k = 1;
                p.A       = 1;
            end

            assert(~isempty(self.Lr2Interpolant), ...
            'Lr2Interpolant not initialized.');

        
            % ========================================================
            % Rossby radius of deformation squared
            % ========================================================

            % Old, simple and performance killer method
            % Lr2_ = (self.g .* self.h) ./ (self.f0 ^ 2);
            % 
            % Lr2_func = @(k,jInd) interp2(self.KRadialLog, self.j, Lr2_, ...
            %                              k, jInd, 'linear');

            % New method (Jan 6th, 2026) 
            % ========================================================
            % Rossby radius interpolant handles
            % ========================================================
            % Safe wrapper (vectorized in k)
            
            % self.Lr2 = @(k,j) self.Lr2Interpolant(j .* ones(size(k)), k );
            % self.Lr2  = @(k,j) self.Lr2Interpolant(j, k);
            kstar2 = @(k,Lr) 1 ./ (Lr*Lr);

            % ========================================================
            % Unnormalized generalized spectrum
            % ========================================================
            S_unnorm = @(k,Lr) Lr ./ ...
                ( (k.^2 ./ kstar2(k,Lr) + 1).^p.slope_k .* ...
                  (Lr*Lr ./ self.Lr2(k,jInd) + 1).^p.slope_j );
        
            % ========================================================
            % Normalize and apply amplitude scaling
            % ========================================================
            S_normalized = self.normalizeSpectrum(S_unnorm);
            S = @(k,jInd) p.A * S_normalized(k,jInd);
        
        end


        function S = frequencyLocalizedSpectrum(self, p)
        % frequencyLocalizedSpectrum  Narrowband frequency–mode localized spectrum
        %
        %   S = frequencyLocalizedSpectrum(self, p)
        %
        %   Returns a function handle S(k,j) representing a spectrum that is
        %   localized around a target frequency omega0 and a target vertical
        %   mode j0. The localization is implemented using Lorentzian profiles
        %   in both frequency and vertical mode number.
        %
        %   This function is intended as a *generic building block* for constructing
        %   spectra associated with specific physical processes (e.g. tides,
        %   inertial oscillations, near-inertial waves). Specific cases should be
        %   implemented as wrapper functions that set omega0 appropriately.
        %
        %   Inputs
        %   ------
        %   self : igw object
        %       Toolbox object containing dispersion relation and grids.
        %
        %   p : struct with fields (all optional)
        %       omega0 : central frequency [rad/s]
        %           Default corresponds to the M2 tidal frequency.
        %       A      : amplitude scaling factor
        %       c      : frequency bandwidth (half-width of Lorentzian) [rad/s]
        %       j0     : central vertical mode number
        %       d      : mode-number bandwidth (half-width of Lorentzian)
        %
        %   Output
        %   ------
        %   S : function_handle
        %       Spectrum S(k,j) defined for continuous horizontal wavenumber k
        %       and discrete vertical mode index j.
        %
        %   Example
        %   -------
        %   % Generic frequency-localized spectrum
        %   S = igw.frequencyLocalizedSpectrum();
        %
        %   % Inertial spectrum (wrapper example)
        %   p.omega0 = igw.f0;
        %   S = igw.frequencyLocalizedSpectrum(p);
        
            arguments
                self
                p.omega0 = 2*pi/(12.42*3600);   % Default: M2 tidal frequency
                p.A      = 100;
                p.c      = 7e-6;
                p.j0     = 3.5;
                p.d      = 0;
            end
        
            % ==================================================
            % Interpolated dispersion relation omega(k,j)
            % ==================================================
            omegaFunc = @(k, jInd) interp1(self.KRadialLog, ...
                                          self.omega(jInd, :), ...
                                          k);
        
            % ==================================================
            % Define frequency–mode localized spectrum
            % ==================================================
            S = @(k, jInd) ...
                (p.A^2 * p.c^2) ./ ...
                ( (omegaFunc(k, jInd) - p.omega0).^2 + p.c^2 ) .* ...
                ( 1 ./ ( (jInd - p.j0).^2 + p.d^2 ) );
        end
       

        function S = tidalSpectrum(self, p)
            arguments
                self
                p.A = 100;
                p.c = 7e-6;
                p.j0 = 3.5;
                p.d  = 0;
            end
        
            p.omega0 = 2*pi/(12.42*3600); % M2
            S = self.frequencyLocalizedSpectrum(p);
        end


        function S = fM2Spectrum(self,p)
            arguments
                self
                p.A = 100;
                p.c = 10;
                p.j0 = 2;
                p.d  = 0.5;
            end
        
            p.omega0 = 2*pi/(12.42*3600) + self.f0;
            S = self.frequencyLocalizedSpectrum(p);
        end


        function S = fSpectrum(self, p)
            arguments
                self
                p.A = 100;
                p.c = 10;
                p.j0 = 0;
                p.d  = 3;
            end
        
            p.omega0 = self.f0;
            S = self.frequencyLocalizedSpectrum(p);
        end


        function spectra = listAvailableSpectra(self)
        % LISTAVAILABLESPECTRA  List spectrum constructors available in the toolbox
        %
        %   spectra = LISTAVAILABLESPECTRA(self) returns a struct array describing
        %   the spectral models implemented in the toolbox.
        %
        %   If called with no output, the list is printed to the command window.
        %
        % -------------------------------------------------------------------------
        % OUTPUT
        %
        %   spectra : struct array with fields
        %       - Name        : Function name
        %       - Type        : 'generic' or 'wrapper'
        %       - Description : Short physical description
        %       - Parameters  : Tunable parameters exposed to the user
        %
        % -------------------------------------------------------------------------
        % EXAMPLE
        %
        %   igw.listAvailableSpectra();
        %
        %   specs = igw.listAvailableSpectra();
        %   {specs.Name}
        %
        % -------------------------------------------------------------------------
        
            spectra = struct( ...
                'Name', {}, ...
                'Type', {}, ...
                'Description', {}, ...
                'Parameters', {} );
        
            % ==================================================
            % Generic building blocks
            % ==================================================
            spectra(end+1) = struct( ...
                'Name', 'generalSpectrum', ...
                'Type', 'generic', ...
                'Description', 'General separable k–j spectrum with tunable slopes and reference mode.', ...
                'Parameters', '{{j_star, slope_j, slope_k, A}}' );
        
            spectra(end+1) = struct( ...
                'Name', 'frequencyLocalizedSpectrum', ...
                'Type', 'generic', ...
                'Description', 'Spectrum localized in frequency and vertical mode using Lorentzian profiles.', ...
                'Parameters', '{{omega0, A, c, j0, d}}' );
        
            % ==================================================
            % Physical wrappers
            % ==================================================
            spectra(end+1) = struct( ...
                'Name', 'tidalSpectrum', ...
                'Type', 'wrapper', ...
                'Description', 'M2 tidal spectrum (wrapper around frequencyLocalizedSpectrum).', ...
                'Parameters', '{{A, c, j0, d}}' );
        
            spectra(end+1) = struct( ...
                'Name', 'fSpectrum', ...
                'Type', 'wrapper', ...
                'Description', 'Near-inertial spectrum centered at the Coriolis frequency f0.', ...
                'Parameters', '{{A, c, j0, d}}' );
        
            spectra(end+1) = struct( ...
                'Name', 'fM2Spectrum', ...
                'Type', 'wrapper', ...
                'Description', 'Frequency-shifted M2 spectrum centered at f0 + omega_M2.', ...
                'Parameters', '{{A, c, j0, d}}' );
        
            % ==================================================
            % Print nicely if no output is requested
            % ==================================================
            if nargout == 0
                fprintf('\nAvailable spectra:\n');
                fprintf('------------------\n');
                for i = 1:numel(spectra)
                    fprintf('• %s (%s)\n', spectra(i).Name, spectra(i).Type);
                    fprintf('    %s\n', spectra(i).Description);
                    fprintf('    Parameters: %s\n\n', spectra(i).Parameters);
                end
                clear spectra
            end
        end


        function self = removeAllEnergy(self)
        % Need to organize function
            self.TE    = [];
            self.HKE   = [];
            self.VKE   = [];
            self.PE    = [];
            self.A2    = [];
            self.TEPSD = [];
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Post-processing / Diagnostic
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function vertVar = verticalVariance(self, energyTerm, options)
        % VERTICALVARIANCE  Compute vertical profile of integrated energy variance
        %
        %   vertVar = VERTICALVARIANCE(self, energyTerm) computes the vertical
        %   variance profile of the selected energy component by interpolating
        %   energy onto a specified vertical grid and summing over modes and
        %   horizontal wavenumbers.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %
        %   self : Model object containing energy fields:
        %          HKE, VKE, PE, and spectral grids.
        %
        %   energyTerm : string
        %       Energy component to process. One of:
        %           'TE'  - Total energy (HKE + VKE + PE)
        %           'HKE' - Horizontal kinetic energy
        %           'VKE' - Vertical kinetic energy
        %           'PE'  - Potential energy
        %
        %   options : name-value arguments
        %       zVector : vertical grid for interpolation (default: linspace(-Lz,0,1000))
        %
        % -------------------------------------------------------------------------
        % OUTPUT
        %
        %   vertVar : column vector
        %       Vertical profile of energy variance evaluated on options.zVector.
        %
        % -------------------------------------------------------------------------
        
            arguments
                self
                energyTerm (1,:) char {mustBeMember(energyTerm, {'TE','HKE','VKE','PE'})} = 'TE'
                options.zVector = linspace(-self.Lz, 0, 1000)
            end
        
            % ==================================================
            % Select energy component
            % ==================================================
            switch energyTerm
                case 'TE'
                    energy = self.HKE + self.VKE + self.PE;
                case 'HKE'
                    energy = self.HKE;
                case 'VKE'
                    energy = self.VKE;
                case 'PE'
                    energy = self.PE;
            end
        
            % ==================================================
            % Interpolate energy onto requested vertical grid
            % ==================================================
            energyInterp = scatteredInterpolation( ...
                    self, energy, options.zVector, self.KRadialLog, 1:self.nModes);

        
            % ==================================================
            % Sum over modes and horizontal wavenumbers
            % ==================================================
            vertVar = sum(sum(energyInterp, 3), 2);
        
        end

        function E_k = energyAtHorizontalWavenumber(self, z, energyTerm, options)
        % ENERGYATHORIZONTALWAVENUMBER  Horizontal wavenumber spectrum at specified depth
        %
        %   E_k = ENERGYATHORIZONTALWAVENUMBER(self, z, energyTerm) returns the
        %   energy spectrum as a function of horizontal wavenumber at a given vertical
        %   position z.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %
        %   self       : igw object containing energy fields and grids
        %   z          : scalar vertical position [m]
        %   energyTerm : string specifying energy component, one of:
        %                'TE', 'HKE', 'VKE', 'PE'
        %   options    : name-value arguments
        %       KRadial : vector of wavenumbers to interpolate onto (default: [])
        %      
        % -------------------------------------------------------------------------
        % OUTPUT
        %
        %   E_k : row vector
        %       Horizontal wavenumber spectrum at depth z.
        %       Size: 1 × nK (or 1 × length(KRadial) if provided)
        %
        % -------------------------------------------------------------------------
        
            arguments
                self
                z (1,1) double
                energyTerm (1,:) char {mustBeMember(energyTerm, {'TE','HKE','VKE','PE'})} = 'TE'
                options.KRadial double = []
            end
        
            % ==================================================
            % Select energy component
            % ==================================================

            switch energyTerm
                case 'TE'
                    energy = self.HKE + self.VKE + self.PE;
                case 'HKE'
                    energy = self.HKE;
                case 'VKE'
                    energy = self.VKE;
                case 'PE'
                    energy = self.PE;
            end
        
            % ==================================================
            % Sum over vertical modes 
            % ==================================================

            E_zK = squeeze(sum(energy, 2));     % [nZ x nK]
            KVec = self.KRadialLog;            
        
            % ==================================================
            % Interpolate to requested z
            % ==================================================
                       
            E_at_z = zeros(1, self.nK); 
        
            for iK = 1:self.nK
                E_at_z(iK) = interp1(self.zPerMode(:, iK), E_zK(:, iK), z);
            end
        
            % ==================================================
            % Normalize by dKLog
            % ==================================================
            
            E_at_z = E_at_z ./ self.dKLog;
            
        
            % ==================================================
            % Interpolate onto user-specified KRadial if provided
            % ==================================================

            if ~isempty(options.KRadial)
                E_k = interp1(KVec, E_at_z, options.KRadial);
            else
                E_k = E_at_z;
            end

        end
            

        function E_j = energyAtVerticalMode(self, z, options)
        % ENERGYATVERTICALMODE  Energy per vertical mode at specified depth
        %
        %   energyAtVerticalMode = ENERGYATVERTICALMODE(self, z, options)
        %   interpolates HKE energy at the scalar vertical position z, then sums
        %   over horizontal wavenumbers. Returns energy for each vertical mode.
        %
        % -------------------------------------------------------------------------
        % NOTE
        %   Currently, this function only works for scalar z.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %
        %   self         : igw object containing HKE and zPerMode
        %   z            : scalar vertical position [m]
        %   options.modeVector : vector of mode indices to include (default: 1:nModes)
        %
        % -------------------------------------------------------------------------
        % OUTPUT
        %
        %   energyAtVerticalMode : array [1 x length(modeVector)]
        %       Energy interpolated at z for each vertical mode
        %
        % -------------------------------------------------------------------------
        
            arguments
                self
                z (1,1) double   % scalar only
                options.modeVector double = 1:self.nModes
            end
            % ==================================================
            % Initialize
            % ==================================================
            nModesLocal = length(options.modeVector);
            HKESamez = zeros(1, nModesLocal, self.nK);
        
            % ==================================================
            % Interpolate HKE along vertical for each mode and horizontal wavenumber
            % ==================================================
            for iK = 1:self.nK
                for jIdx = 1:nModesLocal
                    jLocal = options.modeVector(jIdx);
                    HKESamez(1, jIdx, iK) = interp1(self.zPerMode(:, iK), self.HKE(:, jLocal, iK), z);
                end
            end
        
            % ==================================================
            % Sum over horizontal wavenumber
            % ==================================================
            E_j = squeeze(sum(HKESamez, 3));
        
        end

         
        function [S, EnergyFrequency] = energyAtFrequencies(self, z, energyTerm, options)
        % ENERGYATFREQUENCIES  Compute energy distribution across frequencies
        %
        %   [S, EnergyFrequency] = ENERGYATFREQUENCIES(self, z, energyTerm, options)
        %   interpolates the selected energy term at a given vertical position z,
        %   then redistributes it over the frequency vector omegaVector according
        %   to the object's dispersion relation.
        %
        %------------------------------------------------------------------
        % NOTE
        %   - Currently works for scalar z only.
        %   - The redistribution assumes a fixed mapping from horizontal
        %     wavenumber to frequency using the existing dispersion relation.
        %   - The option 'spectrumType' is included for future development         
        %     but is not yet active.
        %   - The frequency-bin masking logic may be refined in future updates.
        %
        %------------------------------------------------------------------             
        % INPUTS
        %
        %   self        : igw object containing HKE, VKE, PE, zPerMode, KRadialLog, omega
        %   z           : scalar vertical position [m]
        %   energyTerm  : string, one of 'TE', 'HKE', 'VKE', 'PE'
        %   options     : struct with fields
        %       .omegaVector   : vector of target frequencies [rad/s]
        %       .spectrumType  : optional, currently not used
        %
        %------------------------------------------------------------------
        % OUTPUTS
        %
        %   S               : [1 x length(omegaVector)-1] array, total energy per frequency bin
        %   EnergyFrequency : [nModes x length(omegaVector)] array, energy per mode and frequency bin
        %
        %------------------------------------------------------------------
        
            arguments
                self
                z (1,1) double      % scalar only
                energyTerm char
                options.omegaVector double = linspace(self.f0, 0.8*sqrt(self.N2max), self.nK)
                options.spectrumType                
            end
        
            % ==================================================
            % Select energy term and interpolate to z
            % ==================================================
            switch energyTerm
                case 'TE'
                    data = self.HKE + self.VKE + self.PE;
                case 'HKE'
                    data = self.HKE;
                case 'VKE'
                    data = self.VKE;
                case 'PE'
                    data = self.PE;
                otherwise
                    error('Invalid energyTerm. Must be ''TE'', ''HKE'', ''VKE'' or ''PE''.');
            end
        
            % Interpolate to z for all modes and horizontal wavenumbers
            energy = squeeze(scatteredInterpolation(self, data, z, self.KRadialLog, 1:self.nModes));
        
            % ==================================================
            % Initialize frequency energy matrix
            % ==================================================
            
            nOmega = length(options.omegaVector);
            EnergyFrequency = zeros(self.nModes, nOmega);
        
            % ==================================================
            % Redistribute energy over frequency bins
            % ==================================================
            for indj = 1:self.nModes
                for i = 1:(nOmega-1)
                    % Find horizontal wavenumbers corresponding to frequency bin
                    indForOmega = self.omega(indj,:) >= options.omegaVector(i) & ...
                                  self.omega(indj,:) < options.omegaVector(i+1);
                    % Sum energy for this mode and frequency bin
                    EnergyFrequency(indj, i) =  EnergyFrequency(indj,i) + sum(squeeze(energy(indj, indForOmega)));
                end
            end
            
        
            % ==================================================
            % Total energy per frequency (sum over modes)
            % ==================================================
            S = sum(EnergyFrequency, 1);
        
        end


       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Interpolation        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function SIModes = initScatteredInterpolant(self, data, KRadialLog, nModes)
        % INITSCATTEREDINTERPOLANT  Initialize scattered interpolants for energy fields
        
            arguments
                self
                data
                KRadialLog
                nModes
            end
        
            SIModes = cell(1, self.nModes);

            %---------------------------------------------------
            % Vertical coordinate
            %---------------------------------------------------            
            ZVectorLog=reshape(self.zPerMode,[],1);
        
            % --------------------------------------------------
            % Horizontal wavelength coordinate
            % --------------------------------------------------
            KVectorRepLog   = reshape(repmat(KRadialLog, [self.nZ 1]), [], 1);
            lambdaVectorLog = (2*pi) ./ KVectorRepLog;
        
            % --------------------------------------------------
            % Mode-dependent data: [nZ x nModes x nK]
            % --------------------------------------------------
            if size(data,3) > 1
        
                for n = nModes        
                    dataVector = reshape(data(:, n, :), [], 1);
        
                    SIModes{n} = scatteredInterpolant( ...
                        ZVectorLog, lambdaVectorLog, dataVector);
                end
        
            % --------------------------------------------------
            % Mode-independent data
            % --------------------------------------------------
            else
                dataVector = reshape(data, [], 1);
        
                SIModes = scatteredInterpolant( ...
                    ZVectorLog, lambdaVectorLog, dataVector);                   
            end
        end
   

        function DataInterpMat = scatteredInterpolation(self, data, zVectorNew, KVectorNew, verticalMode)
        % SCATTEREDINTERPOLATION  Interpolate energy fields onto a new (z, k) grid
        %
        %   DataInterpMat = SCATTEREDINTERPOLATION(self, data, zVectorNew, KVectorNew, verticalMode)
        %   interpolates the input energy field defined on the native (zPerMode, KRadialLog)
        %   grid onto a new vertical and horizontal wavenumber grid using scattered
        %   interpolation in (z, lambda) space, where lambda = 2*pi/k.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %
        %   self          : igw object containing grid information (zPerMode, nZ, etc.)
        %
        %   data          : energy array defined on the native grid.
        %                   Expected size:
        %                     - [nZ x nModes x nK] for mode-dependent fields
        %                     - [nZ x nK]          for mode-independent fields
        %
        %   zVectorNew    : vector of target vertical positions [m]
        %
        %   KVectorNew    : vector of target horizontal wavenumbers [rad/m]
        %
        %   verticalMode  : vector of vertical modes to interpolate
        %
        % -------------------------------------------------------------------------
        % OUTPUT
        %
        %   DataInterpMat : interpolated energy array of size
        %                   [length(zVectorNew) x length(verticalMode) x length(KVectorNew)]
        %
        %                   Dimension order:
        %                       1st : vertical coordinate (z)
        %                       2nd : vertical mode
        %                       3rd : horizontal wavenumber
        %
        % -------------------------------------------------------------------------
        % NOTES
        %
        %   - Interpolation is performed in (z, lambda) space, where lambda = 2*pi/k.
        %   - Linear interpolation is used with no extrapolation.
        %   - The native vertical grid may vary with mode through zPerMode.
        %
        % -------------------------------------------------------------------------
        
            arguments
                self
                data
                zVectorNew
                KVectorNew
                verticalMode
            end
            
            % ==================================================
            % Initialize interpolants for requested modes
            % ==================================================
            SIModes = initScatteredInterpolant(self, data, KVectorNew, verticalMode);
        
            % ==================================================
            % New (z, k) grid
            % ==================================================
            ZVectorLin = reshape(repmat(zVectorNew,[1, length(KVectorNew)]),[],1);
            lengthZ    = length(zVectorNew);
            KVectorRepNew = reshape(repmat(KVectorNew,[lengthZ 1]),[],1);
            lambdaVectorNew = (2*pi)./KVectorRepNew;
        
            % ==================================================
            % Preallocate
            % ==================================================
            DataInterpMat1 = zeros(lengthZ, length(KVectorNew), length(verticalMode));
        
            % ==================================================
            % Interpolate
            % ==================================================
            if size(data,3) > 1
                for i = 1:length(verticalMode)
                    DataInterp = SIModes{verticalMode(i)}(ZVectorLin, lambdaVectorNew);
                    DataInterpMat1(:,:,i) = reshape(DataInterp, lengthZ, length(KVectorNew));
                end
            else
                for i = 1:length(verticalMode)
                    DataInterp = SIModes(ZVectorLin, lambdaVectorNew);
                    DataInterpMat1(:,:,i) = reshape(DataInterp, lengthZ, length(KVectorNew));
                end
            end
        
            % ==================================================
            % Output as [z x mode x k]
            % ==================================================
            DataInterpMat = permute(DataInterpMat1,[1,3,2]);
        
        end


       
        function DataInterp2D = interp2D(self, data, KVectorNew, verticalMode)
        % INTERP2D  Interpolate mode–wavenumber data onto a new grid
        %
        %   DataInterp2D = INTERP2D(self, data, KVectorNew, verticalMode)
        %   interpolates data defined on the native (mode, KRadialLog) grid onto
        %   a new grid specified by verticalMode and KVectorNew.
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %
        %   self         : igw object containing KRadialLog and nModes
        %   data         : [nModes x nK] array defined on native grid
        %   KVectorNew   : vector of horizontal wavenumbers
        %   verticalMode : vector of vertical mode indices
        %
        % -------------------------------------------------------------------------
        % OUTPUT
        %
        %   DataInterp2D : [length(verticalMode) x length(KVectorNew)] array
        %
        % -------------------------------------------------------------------------
        
            arguments                
                self 
                data
                KVectorNew
                verticalMode
            end
        
            % ==================================================
            % Native grid
            % ==================================================
            [X, Y] = ndgrid(1:self.nModes, self.KRadialLog);
        
            % ==================================================
            % Query grid
            % ==================================================
            [Xq, Yq] = ndgrid(verticalMode, KVectorNew);
       
            % ==================================================
            % 2D interpolation in (mode, k) space
            % ==================================================
            DataInterp2D = interpn(X, Y, data, Xq, Yq, 'linear');
        
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


        
        function binWidth = dKLog(self)
            arguments
                self 
            end

            % ==================================================
            % Ensure KRadialLog is column vector
            % ==================================================
            K = self.KRadialLog(:);
            
            % ==================================================
            % Compute geometric bin edges for log-spaced grid
            % ==================================================
            edges = zeros(self.nK+1,1);
            edges(1) = K(1);
            for i = 1:self.nK-1
                edges(i+1) = sqrt(K(i) * K(i+1));
            end
            edges(end) = K(end);
            
            binWidth = diff(edges)';
            
        end



        function omegaVector = omegaVector(self)
            omegaj1=self.omega(1,:);
            dOmega=max(diff(sort(omegaj1(:))));            
            omegaVector=min(self.omega(:)):2*dOmega:max(self.omega(:));
        end

    end
end

    

