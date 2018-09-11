function y=hfilta(x,fs,fmin,fmax);
%Author: J. Halamek, ISI of the CAS, Czechia, 2000-2018
N=length(x);
    y=fft(x);
    if fmin>0
        imin=fix(fmin/(fs/N));
        if imin>0
        y(1:1+imin)=zeros(size(y(1:1+imin)));
        y(end-imin+1:end)=0;
        else
            y(1)=0;
        end
    else 
        imin=1;
    end
    if fmax<fs/2
        imax=fix(fmax/(fs/N));
  
        y(imax+1:end-imax+1)=0;
       
    end
    y=real(ifft(y));
    

