function [V] = Vfunction(A,t,tstop,f,i)
    switch i
        case{1}
            V = A*sin(2*pi*f*(t))+135E-3;
        case{2}
            V = A*sawtooth(2*pi*f*(t+1/(4*f)), 1/2);  %Start at zero?
        case{3}
            if t < tstop/5
                V = 0;
            elseif (t < 2/5*tstop) && (t >= tstop/5)
                V = -A;
            elseif (t < 3/5*tstop) && (t >= 2/5*tstop)
                V = 0;
            elseif (t < 4/5*tstop) && (t >= 3/5*tstop)
                V = A;
            elseif (t >= 4/5*tstop)
                V = 0;
            end
        case{4}
            if t < tstop/10
                V = 0;
            elseif (t < 2/10*tstop) && (t >= tstop/10)
                V = A;
            elseif (t < 3/10*tstop) && (t >= 2/10*tstop)
                V = 0;
            elseif (t < 4/10*tstop) && (t >= 3/10*tstop)
                V = -A;
            elseif (t < 5/10*tstop) && (t >= 4/10*tstop)
                V = 0;
            elseif (t < 6/10*tstop) && (t >= 5/10*tstop)
                V = A;
            elseif (t < 7/10*tstop) && (t >= 6/10*tstop)
                V = 0;
            elseif (t < 8/10*tstop) && (t >= 7/10*tstop)
                V = -A;
            elseif (t < 9/10*tstop) && (t >= 8/10*tstop)
                V = 0;
            elseif (t >= 9/10*tstop)
                V = A;
            end
        case{5}
            if t < tstop/20
                V = 0;
            elseif (t < 2/20*tstop) && (t >= tstop/20)
                V = A;
            elseif (t < 3/20*tstop) && (t >= 2/20*tstop)
                V = 0;
            elseif (t < 4/20*tstop) && (t >= 3/20*tstop)
                V = -A;
            elseif (t < 5/20*tstop) && (t >= 4/20*tstop)
                V = 0;
            elseif (t < 6/20*tstop) && (t >= 5/20*tstop)
                V = A;
            elseif (t < 7/20*tstop) && (t >= 6/20*tstop)
                V = 0;
            elseif (t < 8/20*tstop) && (t >= 7/20*tstop)
                V = -A;
            elseif (t < 9/20*tstop) && (t >= 8/20*tstop)
                V = 0;
            elseif (t < 10/20*tstop) && (t >= 9/20*tstop)
                V = A;
            else
                V = 0;
            end
        case{6}  %Single Pulse
            if t < tstop/5
                V = 0;
            elseif (t < 1/2*tstop) && (t >= tstop/5)
                V = A;
            else
                V = 0;
            end
         case{7}  %Repeated pulses
            if t < tstop/20
                V = 0;
            elseif (t < 2/20*tstop) && (t >= tstop/20)
                V = A;
            elseif (t < 3/20*tstop) && (t >= 2/20*tstop)
                V = 0;
            elseif (t < 4/20*tstop) && (t >= 3/20*tstop)
                V = A;
            elseif (t < 5/20*tstop) && (t >= 4/20*tstop)
                V = 0;
            elseif (t < 6/20*tstop) && (t >= 5/20*tstop)
                V = A;
            elseif (t < 7/20*tstop) && (t >= 6/20*tstop)
                V = 0;
            elseif (t < 8/20*tstop) && (t >= 7/20*tstop)
                V = A;
            elseif (t < 9/20*tstop) && (t >= 8/20*tstop)
                V = 0;
            elseif (t < 10/20*tstop) && (t >= 9/20*tstop)
                V = A;
            else
                V = 0;
            end
        case{8}  %Fixed voltage
            V=A;
        case{9}  %Short Spike
            if t < tstop/5
                V = 0;
            elseif (t < (tstop*0.20+0.1)) && (t >= tstop*0.20)
                V = A;
            elseif (t < (tstop*0.60+0.1)) && (t >= tstop*0.60)
                V = -A;
            else
                V = 0;
            end
        case{10}  %Read voltage from abf file
            load('VoltageSignal.mat');  %May want to load this in the main file for a single read
            if t < t_seg(end)
                V = interp1(t_seg,V_seg,t)*1E-3;  %Convert from mV to V!
            else
                V = V_seg(end)*1E-3;
            end
        otherwise
            V = 0;
    end
end