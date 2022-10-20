function [Y] = WriteVals(droplet)
    ndroplet = size(droplet,2);  %So how many droplets do we have?
    
    %We are tracking for each droplet:
    %{
    xlocation  (D = v)
    ylocation
    zlocation
    voltage   (D = mess of electrical equations)
    
    %We will also have dC/dt terms.  There are ndroplet*(ndroplet-1) potential
    membranes to manage.  Let's add these after fixing the first terms.
    
    %}
    %initialize the vector.
    Y = zeros(1,ndroplet*7);
    for i = 1:ndroplet
        Y(i) = droplet(i).x;
    end
    for i = 1:ndroplet
        Y(i+ndroplet) = droplet(i).y;
    end
    for i = 1:ndroplet
        Y(i+ndroplet*2) = droplet(i).z;
    end
    for i = 1:ndroplet
        Y(i+ndroplet*3) = droplet(i).voltage;
    end   
end