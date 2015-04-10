function y = linearSW(xr, a, b)

% Get slope and intercept
m = 1 / (b - a);
int = 0 - (m*a);

% Calculate discrete function
y = zeros(1,length(xr));
for i = 1:length(xr)
    if m > 0
        if xr(i) <= a
            y(i) = 0;
        elseif xr(i) >= b
            y(i) = 1;
        else
            y(i) = m*xr(i) + int;
        end
    elseif m < 0
        if xr(i) >= a
            y(i) = 0;
        elseif xr(i) <= b
            y(i) = 1;
        else
            y(i) = m*xr(i) + int;
        end
    end
end