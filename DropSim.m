clear all
close all
%Droplet Dynamics - Coupled Electrowetting and Energetics
%Script for Morphogenesis-Inspired Two-Dimensional Electrowetting in
%Droplet Networks

%Written by:
%Joyce El-Beyrouthy, University of Georgia
%Michelle Mansour, University of Georgia
%Jesse Gulle, University of Georgia
%Eric Freeman, University of Georgia
%Contact Author:  Eric Freeman, ecfreema@uga.edu

%The purpose of this script is to generate a series of droplets, boundary
%conditions, and so on.  Next these are used to create a system of
%differential equations and solve them for droplet motion.

%Droplets are primarily stored in a custom struct, where each droplet
%contains information on its location, velocity, voltage, volume, and any
%other variables of interest.
droplet = struct;

%Most inputs are read from text files.
%The inputs for each case presented in the manuscript are included in
%the supplementary information.  Place them in the same folder as this
%script and execute.

%Required Functions (place in same folder as this script):
%CalcPos
%CalcForces
%WriteVals
%RetrieveVals

%Read sim information
%Read solution information - start time, stop time, etc
[tstart, tstop, DT, DTmax, perturb, nframe, EPS, writevid] = textread('siminfo.txt', '%f %f %f %f %f %f %f %f');

%Read global information.  These are material properties.
[Cspec, Gspec, perm, T, mtens, btens, density, density2, damp, alpha] = textread('globalinfo.txt', '%f %f %f %f %f %f %f %f %f %f');

%Read individual droplet information
[xlocs ylocs zlocs rapps Ks Cls lipid anchor] = textread('dropinfo.txt', '%f %f %f %f %f %f %f %f');

%Read signal information
[Vamp,Vf,Vi] = textread('vinfo.txt', '%f %f %f');

%initialize the droplet locations and properties
%locations are in microns, convert to m.
ndrops = length(xlocs);

%Initialize KREJ
KREJ = 0;

%Wire conductivity - very useful.  This can be much lower than the
%experimental value for most applications.  Since the current input is
%scaled by Gwire, too high of a value will lead to issues with convergence.
Gwire = 1E-7;

for i = 1:ndrops

    %Read initial locations
    droplet(i).x = xlocs(i)*1.0E-6;
    droplet(i).y = ylocs(i)*1.0E-6;
    droplet(i).z = zlocs(i)*1.0E-6;

    %We have to assume some initial conditions to calculate changes.
    %Use these.
    droplet(i).xhist(1) = xlocs(i)*1.0E-6;
    droplet(i).yhist(1) = ylocs(i)*1.0E-6;
    droplet(i).zhist(1) = zlocs(i)*1.0E-6;

    droplet(i).rapp = rapps(i)*1.0E-6;
    droplet(i).vol = droplet(i).rapp^3*4/3*pi;
    droplet(i).mass = droplet(i).vol*density;

    droplet(i).type = lipid(i)*1E-3;

    droplet(i).marea = 0;
    droplet(i).barea = [];
    
    droplet(i).anchor = anchor(i);

    %Change here, dV in calcpos/calcforce and V assignments
    %Initially, there are no fields across the interiors!  Cancel the lipid
    %asymmetry at first.  Let this develop from Vsource.
    droplet(i).voltage = -1*lipid(i)*1E-3; %Initial conditions for V replicate asymmetry
end

%Initialize external forces
for i=1:ndrops
    droplet(i).extfx = 0;
    droplet(i).extfy = 0;
    droplet(i).extfz = 0;
end

%Initialize velocities
for i = 1:ndrops
    droplet(i).xvel = 0;
    droplet(i).yvel = 0;
    droplet(i).zvel = 0;
end

bRun = 1;
%initialize apparent volumes for iteration
for i = 1:ndrops
    droplet(i).volapp = droplet(i).vol;
    droplet(i).volapp2 = droplet(i).vol;
end

time = tstart;
frame = 1;

%Assume S = 1, G = last unless otherwise specified

%Set indices for voltages
S = 1;
G = ndrops;

%Anchor droplets on electrodes
%droplet(S).anchor = 1;
%droplet(G).anchor = 1;

%Before we begin, convert the initialized droplet struct into Y for easy
%management.
Y = WriteVals(droplet);
Y0 = Y;
time0 = time; 
frameDT = tstop/nframe;        %Determine timestep per frame
frame = 0;                      %Initialize frame counter
bCritical = 1;                  %Did it break?

Cmat = zeros(ndrops);           %Initial size is nxn
Gmat = zeros(ndrops);
Cstore = Cmat;                  %Initialize Cstorage.

pT = 0;
pTstore = pT;

bfirst = 0; %Set a flag for when we are calculating the first set of data.
bretry = 0;
%This is the integration loop here.
%Contained inside of this loop is the RK integration step.  We repeat it
%until we've passed tstop.
while time < tstop
    bCritical = 1;        %Set flag - assume RK integration is innocent until proven guilty.  0 for false, 1 for true.  Trip when false.

    %Before integration, record the current state of the droplet clusters
    %for output and  post-processing

    %Set the voltage input.  Several options here, continuous functions are
    %not necessary.
    Vin = Vfunction(Vamp,time0,tstop,Vf,Vi);
    time0   %This is just to display the current time to track progress in the command window
    if time0 >= frame*frameDT  %Record a select number of frames.
        %Store the history for output to a text file etc.
        frame = frame + 1;
        tplot(frame) = time0;
        %time0
        Vplot(frame) = Vin;
        iplot(frame) = Gwire*(Vplot(frame)-(droplet(S).voltage-droplet(G).voltage));
        Crecord(frame,:,:) = Cstore;

        %update droplet history if needed.
        for i = 1:ndrops
            droplet(i).xhist(frame) = droplet(i).x;
            droplet(i).yhist(frame) = droplet(i).y;
            droplet(i).zhist(frame) = droplet(i).z;
            droplet(i).rhist(frame) = (3/4*1/pi*droplet(i).volapp)^(1/3);
            droplet(i).vhist(frame)=droplet(i).voltage;
        end
    end

    %Now integrate forward in time.
    for RK = 1:7

        %Initialize droplet properties
        [droplet, bCritical] = CalcPos(droplet,ndrops,mtens,btens,Cspec,alpha);
        %See https://doi.org/10.1088/0960-1317/16/8/009 for an example of
        %this approach!

        %Now, perturb each droplet forward and backwards if they are not
        %set as anchors (immobile)
        %Key to parallel implementation is minimizing crosstalk.  Each of
        %these calculations are completely independent from each other.
        %Only use reductive assignments within here.
        retforces = CalcForces(droplet,ndrops,mtens,btens,Cspec,alpha,perturb);
        
        for i = 1:ndrops
            droplet(i).xforce = retforces(i,1);
            droplet(i).yforce = retforces(i,2);
            droplet(i).zforce = retforces(i,3);
        end
        
        %Enable this line to gradually increase the forces for
        %intercalation figures
        %droplet(4).xforce = droplet(4).xforce+1E-8*time;
        
        %Repeat for voltages!
        %Each droplet struct has a collection of areas defined, one for
        %each other droplet in the system.  We can use these to define a
        %capacitance matrix of how they are all connected together.
        %Loop through the collection of droplets, and build a row for each
        %out of their respective areas.
        %Update the source voltage in case it has changed
        Vin = Vfunction(Vamp,time,tstop,Vf,Vi);
        eqnset = [];
        %Initialize the voltage vector.
        for i = 1:ndrops
            V(i,1) = droplet(i).voltage;
        end

        %In theory, every droplet can have a connection with every other
        %droplet.
        %We start with an nxn matrix for C and G.
        Cmat = zeros(ndrops);
        Gmat = zeros(ndrops);
        
        %As we trim away droplets and values that are disconnected, we need
        %to keep track of how the original indices change.  We handle this
        %through a lookup vector.  Initialize this with i, and as we remove
        %rows and columns from Cmat/Gmat we can keep track of what the
        %original indices were.
        lookup = [];

        for i = 1:ndrops
            %Build the library list
            lookup(i) = i;
            %Build the capacitance and conductance matrices from the areas.
            for k = 1:length(droplet(i).barea)
                dV = (droplet(i).voltage - droplet(k).voltage);
                dV = dV + (droplet(i).type - droplet(k).type);
                Cs = Cspec*(1+alpha*dV^2);
                Cmat(i,k) = Cmat(i,k)-droplet(i).barea(k)*Cs; %off diagonal
                Cmat(i,i) = Cmat(i,i)+droplet(i).barea(k)*Cs; %diagonal
                Gmat(i,k) = Gmat(i,k)-droplet(i).barea(k)*Gspec; %off diagonal
                Gmat(i,i) = Gmat(i,i)+droplet(i).barea(k)*Gspec; %diagonal
            end
        end

        %dC/dt cannot be calculated directly - we will compare to the value
        %obtained from the first time step in each RK iteration, then use
        %backwards finite difference.
        %If this is the very first frame, we need to assume the previous
        %frame is the same as the current frame.
        if bfirst == 0
            Cstore = Cmat;
            dC = Cstore;
            bfirst = 1;
        else
            if RK == 1
                dC = (Cmat-Cstore)/(time0-pT);
                pT = time0;
                Cstore = Cmat;      %Let's try to use Cstore
            else
                %Now, use the first to establish dC
                dC = (Cmat-Cstore)/(time-time0);  
            end
        end
        %The original Cmat assumes that all droplets share a link.  If the
        %membrane does not exist, we eliminate the link - we can assume
        %this happens when the listed area on the diagonal is zero.
        for i = 1:ndrops
            %Move from end to beginning - if we trim early, we will cause
            %errors.
            locate = ndrops-(i-1);
            if Cmat(locate,locate) == 0 %Disconnected link!
                %Trim the row and column from Cmat
                Cmat(locate,:) = [];
                Cmat(:,locate) = [];
                
                %Do the same for Gmat
                Gmat(locate,:) = [];
                Gmat(:,locate) = [];
                
                %Do the same for dC
                dC(locate,:) = [];
                dC(:,locate) = [];            
                %Now, remove the values for V, and lookup.
                lookup(locate) = [];
                V(locate) = [];
            end
        end

        test1 = size(Gmat); %Blew up?  Sometimes a droplet flies away.  If this is part of the electrical path, this causes problems.

        if test1(1) == 0
            %Failure?
            bCritical = 0;
        elseif isnan(det(Cmat)) == 1
            bCritical = 0;  %We have a singular matrix - the connection between source and ground does not exist
        else
            RHS = -(Gmat+dC)*V;
            %Get the values for the source and ground to calculate current
            %across the wire
            SVolt = V(find(lookup==S,1));
            GVolt = V(find(lookup==G,1));
            
            Gcheck = size(GVolt);       %On occasion the source or ground will be detached, and cannot be found
            Scheck = size(SVolt);       %These lines check whether this is a 0x1, and sets bCritical below accordingly
            
            if Gcheck(1)==0 || Scheck(1)==0
                bCritical = 0;      %Failed - try a smaller timestep
            else
                %Modify the RHS to include current at the source (input, natural boundary condition)
                RHS(find(lookup==S,1)) = RHS(find(lookup==S,1)) + Gwire*(Vin-(SVolt-GVolt));
                
                
                %Now, remove the ground (dVG/dt = 0)
                %so we can invert C.
                %Be careful here - G may have shifted.  Find the new indice.  Use
                %the lookup vector with the original indice for the source and
                %grounds.
                Cmat(find(lookup==G,1),:) = [];
                Cmat(:,find(lookup==G,1)) = []; 

                V(find(lookup==G,1)) = [];
                RHS(find(lookup==G,1)) = [];
                lookup(find(lookup==G,1)) = [];


                %The linear system of differential equations can now be separated -
                %obtain an equation for dV/dt for each droplet

                eqnfinal = zeros(ndrops,1);
            end
            if cond(Cmat) > 1E6
                bCritical = 0;  %If the matrix is still singular, we can check this here.
            else
                eqnset = inv(Cmat)*RHS; 
                %Now, write these back to a vector that contains zeros if there are
                %no changes.  Expand eqnset to include zeros as well.
                for i = 1:length(eqnset)
                    eqnfinal(lookup(i)) = eqnset(i); 
                end
            end
        end
        %Writevals takes the information from the droplet structure and
        %passes it to an array that is suitable for DOPRi5.  Retrievevals does
        %the opposite, writing values from Y1 to droplet.
        %If so inclined, feel free to rewrite the code to remove all of the
        %structures - they may slow it down.
        Y = WriteVals(droplet);

        %Now we need to convert these into differential equations for the
        %solver.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %The total number of differential equations is ndrops*4.  7
        %variables tracked for each droplet, x, y, z, and V.  Add more
        %terms here as needed for osmotic swelling, change in
        %concentrations, etc.

        %This next section compares the motion of each droplet to their
        %neighbors.  Experimentally, the rate of membrane formation and
        %separation varies - we propose varying the damping factor
        %accordingly.
        dropvec = zeros(ndrops,3);
        forcevec = zeros(ndrops,3);
        adjdamp = zeros(ndrops,1);
        for i = 1:ndrops
            denom = 0;
            for j = 1:ndrops
                if (Cstore(i,j) ~= 0) && (i ~=j)  %We have a link.
                    %Create the vector between the droplets
                    dropdrop = [droplet(j).x-droplet(i).x, droplet(j).y-droplet(i).y, droplet(j).z-droplet(i).z];
                    dropvec(i,:) = dropvec(i,:)-dropdrop*Cstore(i,j);  %Use the membrane dimensions to adjust weight factor
                end
            end
            dropvec(i,:) = dropvec(i,:)/norm(dropvec(i,:));  %This describes the typical location of the 
            
            forcevec(i,:) = [droplet(i).xforce, droplet(i).yforce, droplet(i).zforce];
            forcevec(i,:) = forcevec(i,:)/norm(forcevec(i,:));
            adjdamp(i) = (1-0.0*(dot(dropvec(i,:),forcevec(i,:)))); %This returns a factor that either increases or decreases depending on the motion of the droplet
            adjdamp(i) = 1;
        end
        
        
        
        for i = 1:ndrops  %change in x location
            if droplet(i).anchor == 0
                D(i) = (droplet(i).extfx+droplet(i).xforce)/(adjdamp(i)*damp);  %We can calculate xvel directly.
            else
                D(i) = 0;
            end
        end
        for i = 1:ndrops  %change in y location
            if droplet(i).anchor == 0
                D(i+ndrops) = (droplet(i).extfy+droplet(i).yforce)/(adjdamp(i)*damp);
            else
                D(i+ndrops) = 0;
            end
        end
        for i = 1:ndrops  %change in z location (disable for now)
            if droplet(i).anchor == 0
                D(i+ndrops*2) = (droplet(i).extfz+droplet(i).zforce)/(adjdamp(i)*damp);
            else
                D(i+ndrops*2) = 0;
            end
        end
        for i = 1:ndrops  %change in voltage
            D(i+ndrops*3) = eqnfinal(i);
        end   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %We now have a collection of variables for each one.
        %This is Dormand Prince, the default solver for ode45 in MATLAB in all of
        %its glory. We write it manually here so we can provide equations
        %that are not continuous.  Droplet motion is not provided by a
        %known differential equation, so we estimate it numerically and
        %integrate forward, adjusting the equations as necessary.

        %The K values are used in the integration - we remember each one.
        nvar = length(D);
        if RK == 1
            for j = 1:nvar
                K1(j) = D(j);
            end
        elseif RK == 2
            for j = 1:nvar
                K2(j) = D(j);
            end
        elseif RK == 3
            for j = 1:nvar
                K3(j) = D(j);
            end
        elseif RK == 4
            for j = 1:nvar
                K4(j) = D(j);
            end
        elseif RK == 5
            for j = 1:nvar
                K5(j) = D(j);
            end
        elseif RK == 6
            for j = 1:nvar
                K6(j) = D(j);
            end
        elseif RK == 7
            for j = 1:nvar
                K7(j) = D(j);
            end
        end

        %Now integrate forward
        %The values for the Butcher tableau seen here may be found at:
        %https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
        if RK == 1
            for j = 1:nvar
                Y1(j) = Y0(j) + DT*0.2*K1(j);
            end
            time = time0 + 0.2*DT;
        elseif RK == 2
            for j = 1:nvar
                Y1(j) = Y0(j) + DT*(3/40*K1(j)+9/40*K2(j));
            end
            time = time0 + 0.3*DT;
        elseif RK == 3
            for j = 1:nvar
                Y1(j) = Y0(j) + DT * ((44 / 45) * K1(j) - (56 / 15) * K2(j) + (32 / 9) * K3(j));
            end
            time = time0 + 0.8*DT;
        elseif RK == 4
            for j = 1:nvar
                Y1(j) = Y0(j) + DT * ((19372 / 6561) * K1(j) - (25360 / 2187) * K2(j) + (64448 / 6561) * K3(j) - (212 / 729) * K4(j));
            end
            time = time0 + 8/9*DT;
        elseif RK == 5
            for j = 1:nvar
                Y1(j) = Y0(j) + DT * ((9017 / 3168) * K1(j) - (355 / 33) * K2(j) + (46732 / 5247) * K3(j) + (49 / 176) * K4(j) - (5103 / 18656) * K5(j));
            end
            time = time0 + DT;
        elseif RK == 6     
            for j = 1:nvar
                Y1(j) = Y0(j) + DT * ((35 / 384) * K1(j) + (500 / 1113) * K3(j) + (125 / 192) * K4(j) - (2187 / 6784) * K5(j) + (11 / 84) * K6(j));
            end
            time = time0 + DT;
        elseif RK == 7  
            %Regular Solution
            for j = 1:nvar
                Y1(j) = Y0(j) + DT* ((35 / 384) * K1(j) + (500 / 1113) * K3(j) + (125 / 192) * K4(j) - (2187 / 6784) * K5(j) + (11 / 84) * K6(j) + 0*K7(j));
            end
            %Alternative Solution
            for j = 1:nvar
                AltSoln(j) = Y0(j)+DT*((5179 / 57600) * K1(j) + (7571 / 16695) * K3(j) + (393/640) * K4(j) - (92097 / 339200) * K5(j) + (187 / 2100) * K6(j) + (1 / 40) * K7(j));
            end
            %Compare the two for the initial error estimate
            for j = 1:nvar
                ErrEstimate(j) = (Y1(j) - AltSoln(j));
            end
        end

        %Check for failures:
        if bCritical == 0
            RK = 8;  %Jump to the end of the RK process.
        else

        for j = 1:nvar
            Y(j) = Y1(j);
        end
        droplet = RetrieveVals(Y,droplet);      %send back to droplet for recalculating D
            %Call function retrievevals to write back into droplet
            %This allows us to calculate new derivatives.
        end
    end
    %We now exit the RK routine.  Check the error - do we advance in time?

    %Calculate the error.  Normalize by a set value - take precautions
    %to not divide by zero when normalizing.
    Err = 0;    %Initialize to zero.  Sum for each variable in Y.
    for i = 1:nvar
        denom = abs(max(Y1(i), Y0(i))); %Scale the error by the value of the integrated variable.
        denom = max(denom, 1E-5);       %Don't divide by zero!
        Err = Err + (ErrEstimate(i)/denom)^2;   %Normalize
    end
    Err = (Err/nvar)^(1/2);
    if abs(Err) > EPS   %EPS is the defined limit.
        bCritical = 0;
    end

    %Allow the new timestep to vary between 1/2 and 2x the previous
    %timestep.
    if bCritical == 0
        DT = 1/2*DT;
    else
        DT = DT*min(2,max(0.5, (EPS/Err)^(1/5)*0.85));
    end
    if DT > DTmax
        DT = DTmax;
    end
    %Now we examine whether or not to update the values
    if bCritical == 1;  %Did we pass?
        KREJ = 0;
        for i = 1:nvar
            Y(i) = Y1(i);       %Write new values into current values
        end
        droplet = RetrieveVals(Y, droplet); %update droplet conditions
        Y0 = Y;                 %update the new starting points
        time0 = time; 
        pTstore = pT;
    else
        KREJ = KREJ + 1;  %This is a rejected step - many possible causes here.
        if KREJ > 1000  %Try 1000 times before giving up entirely.
            disp('Unable to continue - try a smaller time step or adjust the EPS')
            time = tstop+1;
        else  %Reset everything to the last working point.
            for i = 1:nvar
                Y(i) = Y0(i);       %Restore original values.  Reset time?
            end
            droplet = RetrieveVals(Y, droplet); %reset droplet conditions
            time = time0;
            %Need to fix pT here somehow.
            pT = pTstore;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%
end

%Write the results to a space delimited file for analysis later.
header = 'time, ';

for i = 1:ndrops
    strdrop = sprintf('x%d, y%d, z%d, V%d, r%d, ',i,i,i,i,i);
    header = [header strdrop];
end
header = [header 'i '];
header = [header newline];
writefile = fopen('data.txt','w');
fprintf(writefile,header);

%It's possible that we didn't achieve all of the frames - use frame
%instead of nframe.

for j = 1:frame
    line = sprintf('%0.5g, ', tplot(j));
    for i = 1:ndrops
        strvals = sprintf('%0.5g, %0.5g, %0.5g, %0.5g, %0.5g, ',droplet(i).xhist(j), droplet(i).yhist(j), droplet(i).zhist(j), droplet(i).vhist(j),droplet(i).rhist(j));
        line = [line strvals];
    end
    line = [line sprintf('%0.5g', iplot(j))];
    line = [line newline];
    %write
    fprintf(writefile,line);
end
fclose(writefile);

%The capacitance is a little trickier since it is variable in size - save
%this to a .mat file.  Also save the time and voltage inputs for plotting.
savefile = 'output.mat';
save(savefile,'tplot','Vplot','Crecord','Cspec','iplot')


%If we want a video, we have to create figures for each frame then stitch
%them together into a timelapse.
%This is very time consuming, and you cannot click away from the figures as
%they are produced.
%Leave this disabled unless a video is desired.
if writevid == 1
    nframe = frame;
    for i = 1:ndrops
        for j = 1:nframe
            xplot(i,j) = droplet(i).xhist(j);
            yplot(i,j) = droplet(i).yhist(j);
            zplot(i,j) = droplet(i).zhist(j);
            radii(i,j) = droplet(i).rhist(j);
        end
    end

    writerObj = VideoWriter('video.avi');
    writerObj.FrameRate = 10;
    open(writerObj);
    rmax = 0;
    for j = 1:nframe
        for i = 1:ndrops
            if radii(i,j)>rmax
                rmax = radii(i,j);
            end
        end
    end

    %find limits for plot
    xmax = 0;
    xmin = 0;
    ymax = 0;
    ymin = 0;
    zmax = 0;
    zmin = 0;
    for j = 1:nframe
        for i = 1:ndrops
            if xplot(i,j)+rmax >= xmax
                xmax = xplot(i,j)+rmax;
            end
            if xplot(i,j)-rmax <= xmin
                xmin = xplot(i,j)-rmax;
            end
            if yplot(i,j)+rmax >= ymax
                ymax = yplot(i,j)+rmax;
            end
            if yplot(i,j)-rmax <= ymin
                ymin = yplot(i,j)-rmax;
            end
            if zplot(i,j)+rmax >= zmax
                zmax = zplot(i,j)+rmax;
            end
            if zplot(i,j)-rmax <= zmin
                zmin = zplot(i,j)-rmax;
            end
        end
    end

map = [0 1 1];
    for j = 1:nframe
        fig = figure;
        hold on
        colormap(map)
        material shiny
        for i = 1:ndrops
            r = radii(i,j)*1E6;
            a = xplot(i,j)*1E6;
            b = yplot(i,j)*1E6;
            c = zplot(i,j)*1E6;
            [x,y,z] = sphere(100);
            surfl(x*r+a, y*r+b, z*r+c)
            view(2)
        end
        alpha 0.5
        shading interp
        light('Position',[0 0 1])
        axis equal
        xlabel('X (microns)', 'fontsize', 16)
        ylabel('Y (microns)', 'fontsize', 16)
        zlabel('Z (microns)', 'fontsize', 16)
        xlim([xmin xmax]*1E6*1.1)
        ylim([ymin ymax]*1E6*1.1)
        zlim([zmin zmax]*1E6*1.1)
        titlestr = sprintf('Droplet Dynamics t = %4.2f seconds',tplot(j));
        title(titlestr,'FontSize',14, 'fontsize', 16)
        daspect = [1 1 1];
    %     pbaspect = [1 1 1];    
    %     img(j) = figure;
        hold off
        filename = [sprintf('%03d',j) '.png'];
        workingDir = pwd;
        fullname = fullfile(workingDir,'images',filename);
        print(fullname,'-dpng')
    %     imwrite(img(j),fullname)
        close(fig)
    end

    %Now turn into a movie
    outputVideo = VideoWriter(fullfile(workingDir,'video.avi'));
    outputVideo.FrameRate = 10;
    open(outputVideo)
    for j = 1:nframe
        filename = [sprintf('%03d',j) '.png'];
        workingDir = pwd;
        fullname = fullfile(workingDir,'images',filename);
        img = imread(fullname);
        writeVideo(outputVideo,img)
    end
    close(outputVideo)
end







