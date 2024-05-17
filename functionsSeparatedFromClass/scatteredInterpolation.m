function DataInterpMat = scatteredInterpolation(self,data, zVectorNew, KVectorNew, verticalMode)
             arguments                
                 self
                 data
                 zVectorNew
                 KVectorNew
                 verticalMode
             end

             % WHY IS THAT BEEN CALLED SO MANY TIMES?
             SIModes = initScatteredInterpolant(self,data);

            
             %zLin
             ZVectorLin=reshape(repmat(zVectorNew,[1, length(KVectorNew)]),[],1);
            
             %Klin
             KVectorRepLin = reshape(permute(repmat(KVectorNew,[1 length(zVectorNew)]),[2 1]),[],1);
             lambdaVectorLin= (2*pi)./KVectorRepLin;

             if size(data,3)>1

                 for i=1:length(verticalMode)


                     %WHY IS THAT BEEN CALLED SO MANY TIMES?   
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
