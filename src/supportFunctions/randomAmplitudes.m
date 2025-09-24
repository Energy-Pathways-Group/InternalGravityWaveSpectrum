function [Ap,Am] = randomAmplitudes(self,options)
            % returns random amplitude for a valid flow state
            %
            % Returns Ap, Am, A0 matrices initialized with random amplitude
            % for this flow component.  
            %
            % - Topic: Quadratic quantities
            % - Declaration: [Ap,Am,A0] = randomAmplitudes()
            % - Returns Ap: matrix of size [Nj Nkl]
            % - Returns Am: matrix of size [Nj Nkl]
            
            arguments (Input)
                self {mustBeNonempty}
                options.shouldOnlyRandomizeOrientations (1,1) double {mustBeMember(options.shouldOnlyRandomizeOrientations,[0 1])} = 0
            end
            arguments (Output)
                Ap double
                Am double                
            end
            spectralMatrixSize=([self.nModes, self.nK]);

            Ap = zeros(spectralMatrixSize);
            Am = zeros(spectralMatrixSize);
            

            Ap = ((randn(spectralMatrixSize) + sqrt(-1)*randn(spectralMatrixSize))/sqrt(2));          
            if options.shouldOnlyRandomizeOrientations == 1
                Ap = Ap./abs(Ap);
            end


            Am = ((randn(spectralMatrixSize) + sqrt(-1)*randn(spectralMatrixSize))/sqrt(2));          
            if options.shouldOnlyRandomizeOrientations == 1
                Am = Am./abs(Am);
            end

            %Am = conj(Ap); NOPE
               
            

            % First Im gonna try only with the conj

            % Am = ((randn(self.wvt.spectralMatrixSize) + sqrt(-1)*randn(self.wvt.spectralMatrixSize))/sqrt(2));
            % 
            % if options.shouldOnlyRandomizeOrientations == 1
            %     Am(logical(validModes)) = Am(logical(validModes)) ./ abs(Am(logical(validModes)));
            % end

end




            