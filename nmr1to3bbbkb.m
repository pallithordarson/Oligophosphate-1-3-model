function b=nmr1to3bbbkb(htot,Ltot,K11,K12,K13);%this is a subfunction within the
%nmr1to3fitbindXXX functions, it calculates b = concentration of free ligand
%The inputs are Ltot = total conc. of ligand (guest) and htot = total
%conc. of host, K11, K12 and K13, the stepwise binding constants for 1:3 complex
r = size(Ltot,1);%determines the number of datapoints
uu = ones(r,1);%creates a vector full of ones (1) of the same length as data
%this is done to reserve memory and speed up calculation
%Next the coeffiencts for the cubic equation to be solved are determined,
%the cubic equation is a modified from K. A. Connors in Binding Constants, John
%Wiley and Sons, New York, 1987, p. 161 eq. 4.29. See also Tsukube, H.; Furuta, 
%H.; Odani, A.; Takeda, Y.; Kudo, Y.; Liu, Y.; Sakamoto, H.; Kimura, K 
%in Comprehensive Supramolecular Chemistry, ed. Atwood, J. L.; Davis, 
%J. E.D.; MacNicol, D. D.; Vögtle, F.; Lehn, J-M.; Elsevier Science, Oxford, %
%New York, Tokyo, 1996, Vol 8. pp 435
%The vector uu, as explained above, has the value 1 for all points.
%The equation becomes: a1*x^4 + a2*x^3 + a3*x^2 + a4*x + a5 = 0 (eq. 1)
a1 = (uu.*(K11.*K12.*K13));
a2 = (uu.*((K11.*K12)-(Ltot.*K11.*K12.*K13)+(3.*htot.*K11.*K12.*K13)));
a3 = (uu.*(K11-(Ltot.*K11.*K12)+(2.*htot.*K11.*K12)));
a4 = (uu.*(1-(Ltot.*K11)+(htot.*K11)));
a5 = (uu.*(-1.*Ltot));
%a is the matrix were the rows corresponds to the data points and
%the column to the coeffients in the cubic equation
a = [a1 a2 a3 a4 a5];
r = size(a,1);%determs the size of a
b = zeros(r,1);%creates a vector the size of r full of zeros (to reserve memory)
b(1) = 1e-12;%sets the first solution (at Ltot=0) =1e-12 since b = 0 could 
%cause problems due to x/0 errors.
for n = 2:r;%starts a loop which solves the cubic equation row by row
   x = roots(a(n,:));%for each row, x = the three solution of eq. 1
  % display(n);
   idx = imag(x)==0 & real(x)>0;
   if all(idx ~=0);
       b(n) = min(x(idx));
   else
       ff=real(x);
       idxf = imag(ff)==0 & real(ff)>0;
       if all(idxf ==0);
           b(n) = 1e-21;
       else
           b(n) = min(ff(idxf));
       end            
   end
         
end
