function [FWHM] = findFWHM(x, fx)
%findFWHM: Finds the full width half max (FWHM) of a function
%
%	FWHM = findFWHM(x, fx);

if max(abs(fx)) == 1e-31
    FWHM = 0;

else
    fx_n = fx(~isnan(fx));
    [m, n] = max(fx_n(1:end-1));		%	Find maximum value and index
    ind = find(fx>=m/2);	%	Find indicies where I>=max(I)/2
    nl = min(ind);			%	Leftmost index
    nr = max(ind);			%	Rightmost index

    if length(ind)==1
        FWHM = 0;
    else
    
        if fx(nl) ~= fx(nl-1)
            %	Linear interpolate x positions
            xl = (x(nl)-x(nl - 1))*(m/2-fx(nl - 1))/(fx(nl)-fx(nl - 1)) + x(nl - 1);
        else
            xl = x(nl);
        end
        
        if nr == length(fx) || fx(nr) == fx(nr + 1) || isnan(fx(nr + 1))
            xr = x(nr);
        else 
            %	Linear interpolate x positions
            xr = (x(nr + 1)-x(nr))*(m/2-fx(nr))/(fx(nr + 1)-fx(nr)) + x(nr);
        end
        
        %	Get FWHM
        FWHM = xr-xl;
    end
end
end