function [FWHM, E_max] = findFWHM(x, fx, parttype)
%findFWHM: Finds the FWHM and E_max of the first valid peak
%
%	[FWHM, E_max] = findFWHM(x, fx, parttype);

if max(abs(fx)) == 1e-31
    FWHM = 0;
    E_max = NaN;
    return;
end

% 1. Assign threshold based on particle type
if parttype == 0
    noise_threshold = 0.002;
elseif parttype == 1
    noise_threshold = 10^-2;
else
    noise_threshold = 0;
end

% 2. Early-Exit Search for the FIRST valid peak
n = []; % Initialize empty index
for i = 2:(length(fx)-1)
    % Check if point is above threshold AND strictly greater than neighbors
    if fx(i) > noise_threshold && fx(i) > fx(i-1) && fx(i) > fx(i+1)
        n = i;
        m = fx(i);
        break; % Exit the loop immediately! We found the first peak.
    end
end

% Fallback if no true peak exists above the threshold (e.g., a single rising slope)
if isempty(n)
    [m, n] = max(fx);
end

% ---------------------------------------------------------
% Capture the corresponding energy maximum
if parttype == 1
    % PROTONS: Average the max energy with neighbors if their fx is within 0.03 of m
    E_vals = x(n);
    
    % Check left neighbor
    if n > 1 && (m - fx(n-1)) <= 0.03
        E_vals = [E_vals, x(n-1)];
    end
    
    % Check right neighbor
    if n < length(fx) && (m - fx(n+1)) <= 0.03
        E_vals = [E_vals, x(n+1)];
    end
    
    % Calculate the average energy of the qualifying plateau points
    E_max = mean(E_vals);
else
    % ELECTRONS: Use the exact peak
    E_max = x(n);
end
% ---------------------------------------------------------

% 3. Find the left bound by walking backward from the peak
nl = n;
while nl > 1 && fx(nl-1) >= m/2
    nl = nl - 1;
end

% 4. Find the right bound by walking forward from the peak
nr = n;
while nr < length(fx) && fx(nr+1) >= m/2
    nr = nr + 1;
end

% 5. Perform Linear Interpolation
if nl == nr
    FWHM = 0;
else
    % Left interpolate
    % SAFEGUARD: Check if nl is 1 to prevent fx(0) index errors
    if nl == 1 || fx(nl) == fx(nl-1)
        xl = x(nl);
    else
        xl = (x(nl)-x(nl - 1))*(m/2-fx(nl - 1))/(fx(nl)-fx(nl - 1)) + x(nl - 1);
    end
    
    % Right interpolate
    if nr == length(fx) || fx(nr) == fx(nr + 1) || isnan(fx(nr + 1))
        xr = x(nr);
    else 
        xr = (x(nr + 1)-x(nr))*(m/2-fx(nr))/(fx(nr + 1)-fx(nr)) + x(nr);
    end
    
    FWHM = xr - xl;
end

end