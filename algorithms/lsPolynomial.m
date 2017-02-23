% -------------------------------------------------------------------
% METODO: lsPolynomial
% Polynomial interpolation method for Line Search
% Author: Thiago Lima Silva (thiagolims@gmail.com)
% -------------------------------------------------------------------

function [s,xs, fxs] = lsPolynomial(f, x, d)   
   delta    = 0.5;     % fixed increment   
   ratio    = 1.618;    % golden ratio      
   tolerance= 0.5;      % accuracy      

   alphaL = 0; % beginning of the interval
   alphaU = 99; % end of the interval

   f0 = f(x)
   fq = f(x+delta.*d);   
   fq1 = -1;
   fq2 = -1;

   % FASE I: Find the interval containing the optimal solution   
   found = 0;
   q  = 0;      
   if  fq < f0
      q        = q+1;            
      fq1      = fq;
      
      alphaq   = sum(delta*ratio.^(0:q));      
      fq = f(x + alphaq.*d);      
      if fq1 < f0 && fq1 < fq
           found = 1;      
           alphaU = alphaq;
           alphaI = delta;
      else 
         q = q+1;
      end
   else
      found =  1;      
      alphaU = delta;
      alphaI = (alphaU - alphaL)/2;
   end
   
   
   while found ~= 1
      alphaq   = sum(delta*ratio.^(0:q));

      fq2      = fq1;
      fq1      = fq;
      fq = f(x+alphaq.*d);      

      if (fq1 < fq2) && (fq1 < fq)
         alphaL = sum(delta*ratio.^(0:q-2));
         alphaU = alphaq;
         alphaI = fq1;
         found = 1;
      else
         q = q+1;
      end
   end        
   
   fB    =999;
   fBpre =0; 
   
   k = 1;
   alphaI = (alphaU-alphaL)/2;      
   while(abs(fB-fBpre) >= tolerance) && k <= 50
      fL = f(x+ alphaL.*d);
      fI = f(x+ alphaI.*d);
      fU = f(x + alphaU.*d);

      a2 = (1/(alphaU-alphaI))*((fU - fL)/(alphaU-alphaL) - (fI - fL)/(alphaI-alphaL));     
        
      a1 = (fI - fL)/(alphaI-alphaL) - a2*(alphaL + alphaI);
      a0 = fL - a1*alphaL - a2*alphaL^2;
      %alphaB = (-1/2*a2)*a1;
      
      xc = alphaI; xl = alphaL; xr = alphaU;
      fl = fL; fc = fI; fr = fU;
      
      acr = xc - xr; bcr = xc^2 - xr^2;
      arl = xr - xl; brl = xr^2 - xl^2;
      alc = xl - xc; blc = xl^2 - xc^2;
    
      alphaB = 0.5*(bcr*fl + brl*fc + blc*fr)/(acr*fl + arl*fc + alc*fr);
           
      fBpre  = fB;
      
      %qB     = a0 + a1*alphaB + a2*alphaB^2;
      fB = f(x + alphaB.*d);

      if alphaI < alphaB
         if fI < fB
            alphaU = alphaB;
         elseif fB <= fI
            alphaL = alphaI;
            alphaI = alphaB;
         end
      else
         if fI < fB            
            alphaL = alphaB;            
         else
            alphaU = alphaI;
            alphaI = alphaB;                 
         end 
      end
      k = k+1;
  end 
   
  s  = alphaB;
  xs = x + alphaB*d;
  fxs = f(xs);  
end