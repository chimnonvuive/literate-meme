function animation(rAs, rBs, rCs, rDs)
    figure
    sys(length(rAs)) = struct('cdata',[],'colormap',[]);
    
    v = VideoWriter('./outputs/animation.avi');
    v.Quality = 100; open(v);
    
    for i = 1:length(rAs)

%     ----- plot position -----

        link1 = [0,rAs(i)]; 
        link2_AB = [rAs(i),rBs(i)];
        link2_AC = [rAs(i),rCs(i)];
        link2_BC = [rBs(i),rCs(i)]; 
        link4 = [rCs(i),rDs(i)];
        
        hold on; axis equal;
        plot(real(link1), imag(link1), 'bo-');
        plot(real(link2_AB), imag(link2_AB), 'ro-');
        plot(real(link2_AC), imag(link2_AC), 'ro-');
        plot(real(link2_BC), imag(link2_BC), 'ro-');
        plot(real(link4), imag(link4), 'ko-');
        xlim([-150,150]); ylim([-150,150]);
        sys(i) = getframe; clf;
    end
    writeVideo(v,sys); close(v);
end