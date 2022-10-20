function [retdroplet, bCritical] = CalcPos(droplet,ndrops,mtens,btens,Cspec,gamma)

    %Set flags.
    bRun = 1;
    bCritical = 1;
    KREJ = 0;
    
%     if c == 0
%         %do nothing
%     elseif c == 1
%         droplet(idrop).x = droplet(idrop).x + perturb;
%     elseif c == 2
%         droplet(idrop).x = droplet(idrop).x - perturb;
%     elseif c == 3
%         droplet(idrop).y = droplet(idrop).y + perturb;
%     elseif c == 4
%         droplet(idrop).y = droplet(idrop).y - perturb;
%     elseif c == 5
%         droplet(idrop).z = droplet(idrop).z + perturb;
%     elseif c == 6
%         droplet(idrop).z = droplet(idrop).z - perturb;
%     end
    while bRun == 1
        %This is an iterative process as described in the SI.  The while
        %loop will continue until a criteria is met.
            
        %First, reset all areas.  Cover the droplet in monolayer, and
        %define the bilayer area as zero.
        %The bilayer area is split into a second array corresponding to how
        %the droplets are connected.
        for i = 1:ndrops
            droplet(i).marea = 4*pi*((droplet(i).volapp*3/(4*pi))^(1/3))^2;
            for j = 1:ndrops
                droplet(i).barea(j) = 0;
            end
        end
        
        %Next, check for overlaps.
        for i = 1:ndrops
            for j = i:ndrops  %Don't repeat connections - start at i.
                if i ~= j     %And don't do this for i == j.
                    %Distance between droplets centers is calculated.
                    d = ((droplet(i).x-droplet(j).x)^2+(droplet(i).y-droplet(j).y)^2+(droplet(i).z-droplet(j).z)^2)^(1/2);
                    
                    %The current radius is calculated from the apparent
                    %volume of each droplet.
                    r1app = (droplet(i).volapp*3/(4*pi))^(1/3);
                    r2app = (droplet(j).volapp*3/(4*pi))^(1/3);
                    
                    %We have an overlap if the summed apparent radii in the
                    %droplet pair is greater than the distance between
                    %their centers.
                    if d < r1app+r2app
                        
                        %Details for these equations are provided in the
                        %supplementary material.
                        a = 1/(2*d)*((d+r1app-r2app)*(d-r1app+r2app)*(r1app-d+r2app)*(d+r1app+r2app))^(1/2);
                        h1 = r1app-(r1app^2-a^2)^(1/2);
                        h2 = r2app-(r2app^2-a^2)^(1/2);
                        
                        %Increase volapp by the volume of the spherical
                        %cap.
                        droplet(i).volapp2 = droplet(i).volapp2 + pi*h1/6*(3*a^2+h1^2);
                        droplet(j).volapp2 = droplet(j).volapp2 + pi*h2/6*(3*a^2+h2^2);
                        
                        %Reduce marea by the area of the spherical cap
                        droplet(i).marea = droplet(i).marea - pi*(a^2+h1^2);
                        droplet(j).marea = droplet(j).marea - pi*(a^2+h2^2);
                        
                        %And replace with the bilayer area.
                        droplet(i).barea(j) = a^2*pi;
                        droplet(j).barea(i) = a^2*pi;

                    end
                end
            end
        end
        volerr = 0;
        %Check how much this changed from the previous step.  Have we
        %converged?
        for i = 1:ndrops
            volerr = volerr + ((droplet(i).volapp2-droplet(i).volapp)/droplet(i).vol)^2;
        end
        volerr=(volerr/ndrops)^(1/2);
        
        %If we have too much change
        if volerr > 1.0E-8
            KREJ = KREJ + 1;
            for i = 1:ndrops
                %Reset volapp2 to the original volume.
                droplet(i).volapp = droplet(i).volapp2;
                droplet(i).volapp2 = droplet(i).vol;
            end
            if KREJ > 200           %What if we can't find a converged solution?
                bCritical = 0;      %Occasionally this will blow up.
                bRun = 0;           %Set critical fail to true, exit.
            end
        else
            %If we're in good shape, use this as a starting point for the
            %next integration step.
            droplet(i).volapp = droplet(i).volapp2;
            bRun = 0;
            %Calculate energy
            for i = 1:ndrops
                %First add the bilayer energy
                benergy = 0;
                for j = 1:ndrops
                    if i ~= j
                        %Go through the row.  Find areas.
                        if droplet(i).barea(j) ~= 0
                            dV = (droplet(i).voltage - droplet(j).voltage);
                            dV = dV + (droplet(i).type - droplet(j).type);
                            Cs = Cspec*(1+gamma*dV^2);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %This is the electrowetting term.  Reduce the
                            %bilayer tension by 1/2*Cspec*V^2.
                            benergy = benergy + 1/2*(btens-1/2*dV^2*Cs)*droplet(i).barea(j);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end
                    end
                end
                %Monolayer tension is assumed invariant, just keep it here.
                droplet(i).energy = mtens*droplet(i).marea+benergy;
            end
        end
    end
    retdroplet = droplet;

end