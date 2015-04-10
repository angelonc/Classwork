function y = linearSW(xr, a, b)

% Get slope and intercept
m = 1 / (b - a);
int = 0 - (m*a);

% Set exceptions
if m < 0
    setL = 1;
    setR = 0;
else
    setL = 0;
    setR = 1;
end
    

% Calculate discrete function
y = zeros(1,length(xr));
for i = 1:length(xr)
    if xr(i) <= a
        y(i) = setL;
    elseif xr(i) >= b
        y(i) = setR;
    else
        y(i) = m*xr(i) + int;
    end
end