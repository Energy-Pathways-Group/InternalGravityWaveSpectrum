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
