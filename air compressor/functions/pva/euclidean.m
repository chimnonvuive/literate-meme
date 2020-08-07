function distances = euclidean(rX, rY)
    distances = sqrt(sum((rX' - rY').^2));
end

