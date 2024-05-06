classdef InternalGravityWaveSpectrum < handle
    properties (Access = public)
        latitude % Latitude for which the modes are being computed.
        f0 % Coriolis parameter at the above latitude.
        Lz % Depth of the ocean.
        N2 %function_handle
        g

        KRadialLog  % size(k) = nK
        KRadialLin
        wavelengthLinear

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
        test
        delFuncAll

        nModes, nK, nZ

        HKEcoef, VKEcoef, PEcoef

        HKE, VKE, PE
        
        HKEatK, VKEatK, PEatK, TEatK

        A

        zPerModeLog
        zNew


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
                options.latitude (1,1) double = 33 %set condition. How to modify erro mesage? costume validator
                options.nModes (1,1) double = 64
                options.nK (1,1) double = 32
                options.nZ (1,1) double =  65          
            end

            
            if options.latitude >= -5 && options.latitude < 5
                error("Latitude:MustBeAwayEquator","This toolbox does not provide a good solution near the Equator (5°S to 5°N)")
            
            elseif options.latitude > 90 || options.latitude < -90
                error("Latitude:WrongValue","Latitude not valid")
            end


            
            self.N2=N2;  
            self.latitude=options.latitude;          
            self.nModes=options.nModes;
            self.nK=options.nK;
            self.nZ=options.nZ;
            self.g=9.80665;
            self.Lz=Lz;
            % 
            % self.N2=N2Func;  
            % self.latitude=latitude;          
            % self.nModes=64;
            % self.nK=10;
            % self.nZ=65;
            % self.g=9.80665;
            % self.Lz=Lz;
            
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
        
             
        %Lz = length, positive

        zInitial = linspace(-Lz,0,10001);
        N2zInitial= self.N2(zInitial);
        N2max = max(N2zInitial);


       % This function checks if the vector is monotonic, but if the test
       % fails it throws an error. I didn't find a way to change for a
       % warning. 

        validateattributes( N2zInitial, { 'numeric' }, { 'vector', 'increasing' } )


        %The N2Func is used in InternalModesSpectral. So if I correct the
        %vector I would need to create a new function from that?

        self.zInitial=zInitial;
        self.N2zInitial=N2zInitial;
        self.N2max=N2max;
       
                    
        % Step 1.3: Compute the K associated with max(N2)   
       
        im = InternalModesWKBSpectral(N2=self.N2,zIn=[-Lz 0],zOut=zInitial,latitude=self.latitude,nModes=self.nModes);       
        
        
        %Unit test: how to test if this modes were computed rigth?
        %Unit test: plot FiK,GiK related with this mode

        [FInitial,GInitial,h,k] = im.ModesAtFrequency(0.8*sqrt(N2max));
        Kmax= max(k);

        self.FInitial = FInitial;
        self.GInitial = GInitial;

        % Step 1.4: Define KRadial based on Kmin=0, Kmax and nK        
        minOrder = floor(log10(2*pi/Kmax));
        if minOrder<=0
            minOrder=1;
        end

       % KRadial equally spaced in log scale
        wavelengthLog=logspace(minOrder,5,self.nK);     
        KRadialLog=fliplr((2*pi)./wavelengthLog);
        self.KRadialLog = KRadialLog;

        % KRadial equally spaced in linear scale
        wavelengthLin=linspace(10^minOrder,10^5,10000);
        KRadialLin=fliplr(linspace((2*pi)./wavelengthLin(1),(2*pi)./wavelengthLin(end),10000))';        
        self.KRadialLin = KRadialLin; 


        
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
            
            
            %%%%%%%% IMPORTANT!!!!! %%%%%%%%%%%%%

            % the format for FiK and GiK is [depth (nZ), vertical modes (nModes), horiz wavenumber (nK)]
            % the format for hiK and omegaiK is [vertical modes (nModes), horiz wavenumber (nK)]
            % the format for hiK and omegaiK is [vertical modes (nModes), horiz wavenumber (nK)]

            
        end
        self.zPerModeLog =zPerModeLog;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             INTERPOLATION IN LINEAR SCALE              %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Step 1: Interp2 Quadrature points per K
        % ndgrid
        
 
        zNew = linspace(-Lz,0,501);
        GLin=self.scatteredInterpolation(GiK,zNew,self.KRadialLin,(1:self.nModes));
        FLin=self.scatteredInterpolation(FiK,zNew,self.KRadialLin,(1:self.nModes));

        hLin= self.interp2D(hiK,self.KRadialLin,(1:self.nModes));
        omegaLin= self.interp2D(omegaiK,self.KRadialLin,(1:self.nModes));

        %
        self.zNew = zNew;
        self.zPerMode =zPerModeLog;
        self.F = FLin;
        self.G = GLin;
        self.h = hLin;
        self.omega = omegaLin;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 3: Computation of the energy coeficients based on 
        % the squared equations (Jeffrey's paper)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % 
        % % Step 1: Compute coeficients in 2D [nModes,nKLin]

        HKEcoef = (1/4)*(1+ (f0^2./(self.omega.^2)));
        VKEcoef= (1/4)* ((self.KRadialLin)'.^2 .* self.h.^2);        
        PEcoef =  (1/4)* (((self.KRadialLin)'.^2.*self.h.^2)./self.omega.^2);


        % Step 2: Make the coeficients in 3D [nModes,nK,nZ]
        % (not necessary to create the 3D matriz if using write index on .* )

        for i = 1:size(self.G,1)
             HKEcoef3D(i,:,:)=HKEcoef;
             VKEcoef3D(i,:,:)=VKEcoef;
             PEcoef3D(i,:,:)=PEcoef;
        end

        self.HKEcoef=HKEcoef3D;
        self.VKEcoef=VKEcoef3D;
        self.PEcoef=PEcoef3D;


      %   %%%%%%% TEST!!!!! %%%%%%%%%%%
      % 
      % 
      %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %   % Step 4: Computation of energy distribution according 
      %   % with the alternative Internal Wave Spectrum
      %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % 
      %   % Should I define the function here or do something apart?
      %   % for now I will past here and just adapt:
      % 
      %   j_star=3;
      %   slope=1;
      %   GMAmplitude =1;
      % 
      %   % GM Parameters. We will use the same??           
      %   L_gm = 1.3e3; % thermocline exponential scale, meters
      %   invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
      %   E_gm = 6.3e-5; % non-dimensional energy parameter
      %   E_T = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm*GMAmplitude;
      % 
      %   % Compute the proper M normalization
      %   M = (j_star^2 +(1:1024).^2).^((-5/4));
      %   M_norm = sum(M);
      % 
      %   %Create the energy matrix 2D 
      %   totalEnergy = zeros(options.nModes,length(self.KRadialLin));        
      % 
      % 
      % 
      %   % Step 4.1: Distributing the energy %%%
      %   for j=(1:options.nModes-1)    %I need to think better about the inds here!!!
      % 
      %       for i=(1:length(self.KRadialLin)-1)            
      % 
      %           %Defining LR                
      %           LR= sqrt(self.g*self.h(j,i))/f0;
      % 
      %           %Defining Bfunc and B_norm
      %           fun = @(k) (1./(k.^2*LR^2 + 1).^(1*slope))*LR;
      %           B_norm = integral(fun,self.KRadialLin(1),self.KRadialLin(end));
      % 
      %           % Integrate the energy btw 2 Kh
      %           E = E_T*(integral(fun,self.KRadialLin(i),self.KRadialLin(i+1))/B_norm)*(((j^2 + j_star^2).^((-5/4)))/M_norm);
      %           totalEnergy(j,i) = E;     
      %           clear E
      % 
      %       end               
      % 
      %   end 
      % 
      %   A2D = sqrt((totalEnergy./self.h)/2);
      % 
      %  for i = 1:size(zPerMode,1)
      %      A(i,:,:)=A2D; 
      %  end
      % 
      %  self.A = A;
      % 
      %  N2atQuadPoints=self.N2(self.zPerMode);
      % 
      %  for i = 1:length(self.nModes)
      %       N2atQuadPoints3D(:,i,:)=N2atQuadPoints;
      %  end 
      % 
      %  HKE = self.A.^2.*self.HKEcoef.*self.F.^2;
      %  VKE = self.A.^2.*self.VKEcoef.*self.G.^2;
      %  PE= self.A.^2.*self.PEcoef.*self.G.^2.*N2atQuadPoints3D;
      % 
      %  self.HKE=HKE;
      %  self.VKE=VKE;
      %  self.PE=PE;
      % 
      end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Kinetic Energy 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function HKEVariance =  HKEVariance(self,options)
            %
            arguments
                self
                options.zVector = linspace(-self.Lz,0,4000) 
                options.plot logical = true
            end
            
            %Befome summing over K I need to interpolate to the same z
            %Vector. This loop is definitely not the best way of doing it            

            for i= 1:self.nK
                for j=1:self.nModes
                    HKESamez(:,j,i)=interp1(self.zPerMode(:,i),self.HKE(:,j,i),options.zVector);  
                end
            end
            

            %Also can call HKEAtHorizontalWavenumber 

            %Sum over modes and k
            HKEVariance=sum(sum(HKESamez,3),2);

            %%% plot %%%%
            if options.plot ==1

                figure(10)

                plot(HKEVariance*100,options.zVector)
                title("HKE")
                ylabel("Depth [m]")
                xlabel("Variance [cm^2/s^2]")
                grid on
            else
            end
        end


    
        %%%%%%%%%%%%%%%
        
        function HKEAtHorizontalWavenumber = HKEAtHorizontalWavenumber(self, z, options)
                    % 
            arguments
                self
                z (1,1) double %or vector
                options.KRadial double = self.KRadialLinear
                options.plot logical = true

            end


           %Sum over modes
           HKEatzK=sum(self.HKE,2);

            % interp the matriz [nz, nK] for the same position on the
            % vertical (z)

            for i= 1:length(self.KRadialLinear)
                HKEatk(i)=interp1(self.zPerMode(:,i),HKEatzK(:,i),z);             
            end

            %interp on the KRadial vector specified by the user

            HKEAtHorizontalWavenumber= interp1(self.KRadialLinear,HKEatk,options.KRadial);            

            if options.plot ==1

                figure(20)

                loglog(options.KRadial, HKEAtHorizontalWavenumber*100)    
                title("HKE")
                ylabel("Variance [cm^2/s^2]")
                xlabel("Horizontal Wavenumber [m]")
                grid on
            else
            end

             

        end



        %%%%%%%%%%%%%%%

        function HKEAtVerticalMode = HKEAtVerticalMode(self, z, options)
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
            
            HKEAtVerticalMode = HKEatzMode;

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

        function S = HKEAtFrequencies(self,z,options)
            arguments
                self
                z
                options.omegaVector = linspace(self.f0,0.8*sqrt(self.N2max),self.nK);
                options.spectrumType
                options.plot = true
            end
            
            
            
            for i= 1:length(self.KRadialLinear)
                for j=1:self.nModes
                    HKEGivez(j,i)=interp1(self.zPerMode(:,i),self.HKE(:,j,i),z);  
                end
            end

            
            for j=1:self.nModes
                    HKEOmega(j,:)=interp1(self.omega(j,:),HKEGivez(j,:), options.omegaVector);  
            end
           
            S = sum(HKEOmega,1); 

            
            if options.plot ==1

                figure(60)

                loglog((options.omegaVector)*(24*3600)/(2*pi), S*100)    
                title("HKE")
                ylabel("Variance [cm^2/s^2]")
                xlabel("Frequency [cycle/day]")
                grid on
            else
            end     
              
            
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


        %%%%%%%%%%%


        function checkOrthogonalityPlot(self)
            arguments
                self               
            end 

            N2atQuadPoints=self.N2(self.zPerMode(:,2));
            self.g=9.80665;

            figure()
            subplot(1,3,1)
            suptitle('Orthogonality Check Same Mode')

            plot(self.G(:,2,2),self.zPerMode(:,2),'k',LineWidth=1.5)            
            ylabel('depth')  
            ylim([-1000 0])
            title('G_{j=2}')

            subplot(1,3,2)
            plot(self.G(:,2,2),self.zPerMode(:,2),'k',LineWidth=1.5)         
            title('G_{j=2}')          
            ylim([-1000 0])

            subplot(1,3,3)
            plot((N2atQuadPoints.*self.G(:,2,2).*self.G(:,2,2)),self.zPerMode(:,2),'k',LineWidth=1.5)           
            title('N_2G_{j=2}G_{j=2}')
            ylim([-1000 0])
          
            
            
            delFuncSameMode = trapz(self.zPerMode(:,2),N2atQuadPoints.*self.G(:,2,2).*self.G(:,2,2),1);
            disp(['delFunc= ', num2str(delFuncSameMode/self.g)])
            
            figure()
            subplot(1,3,1)
            suptitle('Orthogonality Check Diff Mode')

            plot(self.G(:,2,2),self.zPerMode(:,2),'k',LineWidth=1.5)            
            ylabel('depth')  
            ylim([-1000 0])
            title('G_{j=2}')

            subplot(1,3,2)
            plot(self.G(:,3,2),self.zPerMode(:,2),'k',LineWidth=1.5)         
            title('G_{j=3}')          
            ylim([-1000 0])

            subplot(1,3,3)
            plot((N2atQuadPoints.*self.G(:,3,2).*self.G(:,2,2)),self.zPerMode(:,2),'k',LineWidth=1.5)           
            title('N_2G_{j=2}G_{j=3}')
            ylim([-1000 0])
          

            delFuncDiffMode = trapz(self.zPerMode(:,2),N2atQuadPoints.*self.G(:,3,2).*self.G(:,2,2),1);
            disp(['delFunc= ', num2str(delFuncDiffMode/self.g)])
        end

        %%%%%%%%%%%%%%
        function checkOrthogonalityAllModes(self)

            arguments
                self               
            end   

            for indK=1:length(self.KRadial)
                disp(indK)    
                N2atQuadPoints=self.N2(self.zPerMode(:,indK));
                dz=gradient(self.zPerMode(:,indK));
    
                B=N2atQuadPoints.*self.G(:,:,indK).*dz;  
                BT= transpose(self.G(:,:,indK));
                
                self.delFuncAll = (BT*B)./self.g;

                isidentity=@(a,tol) all(abs(a-eye(size(a)))<tol);

                tol=1; %I think this tolerance is too big

                isOrthogonal=isidentity(self.delFuncAll,tol);
                self.test=isOrthogonal;
                
                if sum(isOrthogonal)<self.nModes-1
                   disp("Vetical structure G is not orthogonal")
                   return                   
                end    
                 disp("Vetical structure G is orthogonal") 
            end
        end

        %%%%%%%%%%%%%%%%
        function plotQuadraturePoints(self,Mode)
            arguments
                self 
                Mode (1,1) integral
            end
        end

        %%%%%%%%%
        function checkEnergySum(self)
            arguments
                self                 
            end

            N2atQuadPoints=self.N2(self.zPerMode);

            for i = 1:length(self.nModes)
                N2atQuadPoints3D(:,i,:)=N2atQuadPoints;
            end

  
           % Integrating in the vertical and summing over modes
           % for each K, the vertical grid is different, so I am doing
           % this computation in a loop, but probably there is a better way
           
           
           for i = 1:length(self.KRadialLinear)
                HKEatK(:,i)= trapz(self.zPerMode(:,i),self.HKE(:,:,i));
                VKEatK(:,i)= trapz(self.zPerMode(:,i),self.VKE(:,:,i));
                PEatK(:,i)= trapz(self.zPerMode(:,i),self.PE(:,:,i));
           end

            self.HKEatK =sum(HKEatK);
            self.VKEatK =sum(VKEatK);
            self.PEatK =sum(PEatK);
            self.TEatK = sum(squeeze(self.A(1,:,:)).^2.*self.h)/2;

            disp(["Total Energy: ",num2str(sum(self.TEatK)), "and the Total " + ...
                "Energy by summation of Energy pieces is: ", num2str(sum(self.HKEatK+self.VKEatK+self.PEatK))])
            
            figure(2)
            
            loglog(self.KRadialLinear,self.HKEatK,LineWidth=1.5) 
            hold on
            loglog(self.KRadialLinear,self.VKEatK,LineWidth=1.5) 
            loglog(self.KRadialLinear,self.PEatK,LineWidth=1.5) 
            loglog(self.KRadialLinear,self.TEatK,LineWidth=1.5) 
            xlim([min(self.KRadialLinear) max(self.KRadialLinear)])
            %xticks(log10(2*pi./[1e5 1e4 1e3 1e2 1e1]))

            ylabel('Energy')  
            xlabel("log(KRadial)")
            legend("HKE","VKE","PE","TE" )

            %%% KRadial needs to be evenly spaced.
            %%% Kmax related to 80% of the maximum stratification is very large 
            %%% (1.4 which generates a wavelength of 4m!!). 
            %%% Equally spacing the vector from 0 to Kmax with 64 generates the following result:
            
            %%% L1=inf
            % L2=473m
            % All the wavelengths are small. How can this be resolved?
            % - Decrease Kmax?
            % - Increase the number of points?


        end

         %%%%%%%%%
         function SIModes = initScatteredInterpolant(self,data)
             arguments                
                 self 
                 data
             end

                % ZLog
                ZVectorLog=reshape(self.zPerModeLog,[],1);

                % KLog
                KVectorRepLog = reshape(repmat(self.KRadialLog,[self.nZ 1]),[],1);
                lambdaVectorLog= (2*pi)./KVectorRepLog;

             if size(data,3)>1
                 for n=1:self.nModes
                     dataVector = reshape(data(:,n,:),[],1);
    
                     SIModes{n}=scatteredInterpolant(ZVectorLog,lambdaVectorLog,dataVector);                          
                 end
             else
                 dataVector = reshape(data,[],1);
    
                 SIModes=scatteredInterpolant(ZVectorLog,lambdaVectorLog,dataVector); 
             end

         end           

         
   

        %%%%%%%%%
        function DataInterpMat = scatteredInterpolation(self, data, zVectorNew, KVectorNew, verticalMode)
             arguments                
                 self 
                 data
                 zVectorNew
                 KVectorNew
                 verticalMode
             end

             SIModes = initScatteredInterpolant(self,data);

            
             %zLin
             ZVectorLin=reshape(repmat(zVectorNew,[1, length(KVectorNew)]),[],1);
            
             %Klin
             KVectorRepLin = reshape(permute(repmat(KVectorNew,[1 length(zVectorNew)]),[2 1]),[],1);
             lambdaVectorLin= (2*pi)./KVectorRepLin;

             if size(data,3)>1

                 for i=1:length(verticalMode)
    
                     DataInterp = SIModes{verticalMode(i)}(ZVectorLin, lambdaVectorLin);
    
                     DataInterpMat1(:,:,i) =reshape(DataInterp, length(zVectorNew),length(KVectorNew));
                 end

             else

                 for i=1:length(verticalMode)
    
                     DataInterp = SIModes(ZVectorLin, lambdaVectorLin);
    
                     DataInterpMat1(:,:,i) =reshape(DataInterp, length(zVectorNew),length(KVectorNew));
                 end

             end
             DataInterpMat = permute(DataInterpMat1,[1,3,2]);

             % Interpolating
                  
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
        



         end


    end

    

