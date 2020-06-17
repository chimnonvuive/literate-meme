function angle_rad = radk2pi(angle_rad)
    % angle_rad must actually be a non-decreasing array.
    % it does not detect error/the case where angle_rad(i)<angle_rad(i+1)
    % mathematically
    left_max = {};
    for i=1:length(angle_rad)-1
        if angle_rad(i)>angle_rad(i+1)
            left_max(end+1) = {i};
        end
    end
    left_max = cell2mat(left_max);
    for k=1:length(left_max)
        if k<length(left_max)
            angle_rad(left_max(k)+1:left_max(k+1)) = angle_rad(left_max(k)+1:left_max(k+1)) + 2*pi*k;
        end
    end
    angle_rad(left_max(end)+1:end) = angle_rad(left_max(end)+1:end) + 2*pi*(k+1);
end