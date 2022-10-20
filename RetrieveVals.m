function [droplet] = RetrieveVals(Y, droplet0)
    ndroplet = size(Y,2)/7;  %So how many droplets do we have?
    
    %We are tracking for each droplet:
    %{
    xlocation  (D = v)
    ylocation
    zlocation
    voltage   (D = mess of electrical equations)
    %}
    
    droplet = droplet0;
    
    for i = 1:ndroplet
        droplet(i).x = Y(i);
    end
    for i = 1:ndroplet
        droplet(i).y = Y(i+ndroplet*1);
    end
    for i = 1:ndroplet
        droplet(i).z = Y(i+ndroplet*2);
    end
    for i = 1:ndroplet
        droplet(i).voltage = Y(i+ndroplet*3);
    end   
end