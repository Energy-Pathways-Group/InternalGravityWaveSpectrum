function Srnd = randomSpectrum(self)

    n=round(self.KRadialLog/min(self.KRadialLog))*10;
    sigma2 = self.A2;
%no need for loop    
    for indj= self.j
        for indK= 1:self.nK
           thisGauss=normrnd(0,sqrt(sigma2(indj,indK)),[1,n(indK)]); 
           thisChi2=thisGauss.^2;
           Srnd(indj,indK) = mean(thisChi2); 
    
        end
    end
    Srnd=(Srnd.*self.h)/2;
%Just proving for to self that this make sense
% a=normrnd(0,2,[1 20]);
% b=a.^2;
% c=mean(b);

end