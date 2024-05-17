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