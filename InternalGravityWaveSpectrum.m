classdef InternalGravityWaveSpectrum < handle
    properties (Access = public)
        latitude % Latitude for which the modes are being computed.
        f0 % Coriolis parameter at the above latitude.
        Lz % Depth of the ocean.
        N2 %function_handle

        KRadial    % size(k) = nK
        F  % size(F_k) = [nZ,nModes,nK]
        G  % size(G_k) = [nZ,nModes,nK]
        h  % size(h_k) = [nModes,nK]
        omega % size(omega_k) = [nK,nModes]       
        zPerMode % [nZ,nModes]
        N2max
        FInitial % [nZ,nModes]
        GInitial % [nZ,nModes]
        zInitial % [nZ]
        N2zInitial

        nModes, nK, nZ
        
        

    end
    
    %properties (GetAccess = private) % I am tring to keep here the variables I want to use for Unit test proposes
        
    %end


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
                options.latitude (1,1) double = 33
                options.nModes (1,1) double = 64
                options.nK (1,1) double = 4
                options.nZ (1,1) double =  5              
            end

            
            self.N2=N2;  
            self.latitude=options.latitude;          
            self.nModes=options.nModes;
            self.nK=options.nK;
            self.nZ=options.nZ;
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 1: Computation of min and max Kh based on the
        % stratification and latitude
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Step 1.1: Compute the Coriolis frequency based on the latitude
        f0 = 2* 7.2921*10^-5 * sind(options.latitude);
        self.f0=f0;

        % Step 1.2: Get the max bouyancy frequency based on the
        % stratification function

        % Step 1.2.1: define vertical vector (z) based on Lz
        % and nZ
        
             
        %Lz = length, positive

        zInitial = linspace(-Lz,0,10001);
        N2zInitial= N2(zInitial);
        N2max = max(N2zInitial);

        self.zInitial=zInitial;
        self.N2zInitial=N2zInitial;
        self.N2max=N2max;
       
                    
        % Step 1.3: Compute the K associated with max(N2)   
       
        im = InternalModesWKBSpectral(N2=N2,zIn=[-Lz 0],zOut=zInitial,latitude=options.latitude,nModes=options.nModes);       
        
        
        %Unit test: how to test if this modes where computed rigth?
        %Unit test: plot FiK,GiK related with this mode

        [FInitial,GInitial,h,k] = im.ModesAtFrequency(0.8*sqrt(N2max));
        %Ks= k.*conj(k); %if its not real, there ia a problem! check h!!!
        Kmax= max(k);

        self.FInitial = FInitial;
        self.GInitial = GInitial;

        % Step 1.4: Define KRadial based on Kmin=0, Kmax and nK
        KRadial = linspace(0,Kmax,options.nK);   
        self.KRadial = KRadial;
                              
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        % Step 2: Computation of F and G Matriz [nK, nZ, nModes]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        upperBoundary = UpperBoundary.rigidLid;
        normalization = Normalization.kConstant;

        for iK=1:length(KRadial)    
    
            im = InternalModesSpectral(N2=N2,zIn=[-Lz 0],zOut=zInitial,latitude=options.latitude,nModes=options.nModes);
            %[FThis,GThis,hThis,omegaThis] = im.ModesAtWavenumber(KRadial(iK));         
        
            im.normalization = normalization;
            im.upperBoundary = upperBoundary;  
            
            %zPerMode(:,iK) = im.GaussQuadraturePointsForModesAtWavenumber(options.nModes+1,KRadial(iK));
            zPerMode(:,iK) = im.GaussQuadraturePointsForModesAtWavenumber(options.nModes+1,KRadial(iK));
            
            im = InternalModesSpectral(N2=N2,zIn=[-Lz 0],zOut=zPerMode(:,iK),latitude=options.latitude,nModes=options.nModes);
            [FiK(:,:,iK),GiK(:,:,iK),hiK(:,iK),omegaiK(:,iK)] = im.ModesAtWavenumber(KRadial(iK)); %modes at quadrature points and not equally spaced
        
        end
        self.zPerMode =zPerMode;
        self.F = FiK;
        self.G = GiK;
        self.h = hiK;
        self.omega = omegaiK;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 3: Computation of the energy coeficients based on 
        % the squared equations (Jeffrey's paper)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % 
        % % Step 1: Compute coeficients in 2D [nModes,nK]
        % 
        % HKEcoef = (1/4)*(1+ (f0^2./(omegaiK.^2)));
        % HVEcoef= (1/4)* (KRadial'.^2 .* h);
        % 
        % % Will the above multiplication work or do I need to do over a
        % % loop?? Or I can also transform KRadial [nK] in a 2D matrix [nModes, nK] 
        % 
        % PEcoef =  (1/4)* ((KRadial.^2.*h)./omegaiK);
        % 
        % self.HKEcoef=HKEcoef;
        % self.HVEcoef=HVEcoef;
        % self.PEcoef=PEcoef;


        % Step 2: Make the coeficients in 3D [nModes,nK,nZ]
        % (not necessary to create the 3D matriz if using write index on .* )

        % for i = 1:length(z)
        %     HKEcoef3D(:,:,i)=HKEcoef;
        %     HVEcoef3D(:,:,i)=HVEcoef;
        %     PEcoef3D(:,:,i)=PEcoef;
        % end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 4: Computation of energy distribution according 
        % with the alternative Internal Wave Spectrum
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %   % Should I define the function here or do something apart?
      %   % for now I will past here and just adapt:
      % 
      %   j_star=3;
      %   slope=1;
      % 
      %   % GM Parameters. We will use the same??           
      %   L_gm = 1.3e3; % thermocline exponential scale, meters
      %   invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
      %   E_gm = 6.3e-5; % non-dimensional energy parameter
      %   E_T = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm*GMAmplitude;
      % 
      %   % Compute the proper vertical function normalization
      %   M = (j_star^2 +(1:1024).^2).^((-5/4));
      %   M_norm = sum(M);
      % 
      %   %Create the energy matrix 3D 
      %   totalEnergy = zeros(length(nK),length(nModes));        
      % 
      % 
      % 
      %   % Step 4.1: Distributing the energy %%%
      %   for j=(1:length(nModes)-1)    %I need to think better about the inds here!!!
      %       %Kh_2D = Kh(:,:,j+1);  
      % 
      %       for i=(1:length(kRadial)-1)            
      % 
      %           %Defining LR                
      %           LR= sqrt(g*h(i,j))/f0;
      % 
      %           %Defining Bfunc and B_norm
      %           fun = @(k) (1./(k.^2*LR^2 + 1).^(1*slope))*LR;
      %           B_norm = integral(fun,kRadial(1),kRadial(end));
      % 
      %           % Integrate the energy btw 2 Kh
      %           E = E_T*(integral(fun,kRadial(i),kRadial(i+1))/B_norm)*(((j^2 + j_star^2).^((-5/4)))/M_norm);
      %           totalEnergy(i,j) = E;     
      %           clear E
      % 
      %       end
      % 
      % 
      % 
      %   end 
      %       % Step 4.2: Get APlus AMinus Matrix ???
      %       % The matrix will have the size [nk, nl, nModes]
      %       % After here I will need to transform to [nK, nModes, nZ]???
      %       % I don't think I need this, there must be other way!    
      % 
      % 
      %       A = sqrt((TotalEnergy./self.h)/2);
      % 
      end
       


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Kinetic Energy 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = HKEVariance(self, z, options)
            % 
            arguments
                self
                z
                options
            end
            
            %Also can call HKEAtHorizontalWavenumber 

            % Step 1: Multiplay G and F by Energy Coeficients and A matrix
            % Step 2: Sum(?) over K and Nmodes to get HKE by z only
            % Step 3: Interpolate in the z vector that the user inputed
            % Step 4: Return the data ou return an error. Like if the user
            % asked for energy in a depth deeper than the local depth.


        end


        
        function E = HKEAtHorizontalWavenumber(self,  z, KRadial, options)
                    % 
            arguments
                self
                z
                KRadial
                options               
            end
           
            % Step 1: Multiplay G and F by Energy Coeficients and A matrix
            % Step 2: Sum(?) over Nmodes and select on depth (DO NOT INTEGRATE OVER DEPTH)
            % Step 3: Interpolate in the kRadial vector that the user inputed
            % Step 4: Return the data ou return an error. 
        end



        function E = HKEAtVerticalMode(self, modeVector, options)
                        % 
            arguments
                self
                modeVector
                options
            end

        % Step 1: Multiplay G and F by Energy Coeficients and A matrix
        % Step 2: Sum(?) over K and integrate over depth
        % Step 3: Interpolate in the modeVector vector that the user inputed
        % Step 4: Return the data ou return an error. 
        end
      

        function S = HKEAtFrequencies(self,omega,spectrumType)
            arguments
                self
                omega
                spectrumType
            end


        % Step 1: Multiplay G and F by Energy Coeficients and A matrix
        % Step 2: Sum(?) over nModes and integrate over depth
        % Step 3: Convert from wave number to frequency. That is the most
        % complicated part, because even that I coded this function already
        % I will need to adapt this. Or would be just to get the correspondent
        % frequency for each K based on dispertion relation??
        % Step 4: Interpolate in the omega vector that the user inputed
        % Step 5: Return the data ou return an error. 
        
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Vertical Kinetic Energy 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for the next section the step-by-step would be the same than for
        % Horizontal Kinetic Energy.

        function E = VKEVariance(self, z, options)
            % 
            arguments
                self
                z
                options
            end

        %Step 1: 

        end


        
        function E = VKEAtHorizontalWavenumber(self, KRadial, options)
                    % 
            arguments
                self
                KRadial
                options
            end
        
        end



        function E = VKEAtVerticalMode(self, modeVector, options)
                        % 
            arguments
                self
                modeVector
                options
            end
    
        end



        function S = VKEAtFrequencies(self,z,omega,spectrumType)
            arguments
                self
                z
                omega
                spectrumType
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Potential Energy 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = PEVariance(self, z, options)
            % 
            arguments
                self
                z
                options
            end

        end


        
        function E = PEAtHorizontalWavenumber(self, KRadial, options)
                    % 
            arguments
                self
                KRadial
                options
            end
        
        end



        function E = PEAtVerticalMode(self, modeVector, options)
                        % 
            arguments
                self
                modeVector
                options
            end
            
        end



        function S = PEAtFrequencies(self,z,omega,spectrumType)
            arguments
                self
                z
                omega
                spectrumType
                
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Tests
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            function plotStratifcationHighMode(self)
                        % 
            arguments
                self                            
            end

            
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


        function checkOrthogonality(self)
            arguments
                self               
            end

            figure()
            subplot(1,3,1)
            plot(self.G(:,2,2),self.zPerMode(:,2),'k',LineWidth=1.5)            
            ylabel('depth')            
            title('G - Mode=2; indK=1')

            subplot(1,3,2)
            plot(self.G(:,2,2),self.zPerMode(:,2),'k',LineWidth=1.5)            
            ylabel('depth')            
            title('G - Mode=3; indK=1')

            subplot(1,3,3)
            plot((self.G(:,3,2).*self.G(:,2,2)),self.zPerMode(:,2),'k',LineWidth=1.5)            
            ylabel('depth')            
            title('GMode=2 times GMode=3; indK=1')
            
            delFunc = trapz(self.zPerMode(:,2),self.G(:,2,2).*self.G(:,2,2),1);

            disp(['delFunc= ', num2str(delFunc)])


        % indK=1;
        % 

        end

        function plotQuadraturePoints(self,Mode)
            arguments
                self 
                Mode (1,1) integral
            end
        end

    end
end

    

