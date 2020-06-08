function [a]=animate(phi,so_bien)
    for i=1:so_bien
        text(cos(phi(i)),sin(phi(i)),[num2str(i)]);
        hold on, axis equal
    end
    a=1;
end