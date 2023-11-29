function y = customFunction(x)
    if x >= 1 && x <= 100
        y = 10 * (x - 1) - 1000;
    elseif x >= 101 && x <= 200
        y = 10 * (x - 100);
    else
        error('Input x is out of the specified range.');
    end
end