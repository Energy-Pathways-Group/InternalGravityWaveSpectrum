function [totalEnergyPerComponent] = amplitudesWithSpectrum(self,spectrum)
            arguments (Input)
                self {mustBeNonempty}
                spectrum {mustBeNonempty}
                
            end
            arguments (Output)
               totalEnergyPerComponent 
            end
                                    
            totalEnergyPerComponent=zeros(self.nModes,self.nK);

                                    
            for iJ = 1:length(self.j)
                for iK = 1:length(self.KRadialLog)
                    if iK == 1
                        % Primeiro intervalo
                        lowerBound = 0.1 * self.KRadialLog(1);  % Use um valor pequeno para evitar zero
                        upperBound = sqrt(self.KRadialLog(iK) * self.KRadialLog(iK + 1));
                    elseif iK == length(self.KRadialLog)
                        % Último intervalo
                        lowerBound = sqrt(self.KRadialLog(iK - 1) * self.KRadialLog(iK));
                        upperBound = self.KRadialLog(iK);
                    else
                        % Intervalos intermediários
                        lowerBound = sqrt(self.KRadialLog(iK - 1) * self.KRadialLog(iK));
                        upperBound = sqrt(self.KRadialLog(iK) * self.KRadialLog(iK + 1));
                    end

                    % Calcula o intervalo dk
                    dk = upperBound - lowerBound;

                    % Integral direta
                    integralValue = integral(@(k) spectrum(k, self.j(iJ)), lowerBound, upperBound)/dk;
                    totalEnergyPerComponent(iJ, iK) = integralValue;                    
                    
                end
            end
            % Ap(isnan(Ap)) = 0;
            % Am(isnan(Am)) = 0;
        end
    