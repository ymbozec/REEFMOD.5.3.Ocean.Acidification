%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: 09/2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ALGALREMOVAL2 = f_algal_removal(algal_cm2, ALGALREMOVAL, total_area_cm2)

% This determines the balance amount of grazing needed at every time step
% -> algal_cm2(:,a) gives the cover (cm2) of each algal type (a) over the grid
% -> ALGALREMOVAL gives the global balance amount of grazing for each algae
% -> total_area_cm2 is the total area of the grid (cm2)

% Need to define this to fit with the original code:
total_algal_cm2 = sum(algal_cm2,1) ; % total cover of each algae over the whole grid
turfinit = total_algal_cm2(1)/total_area_cm2 ;
dictinit = total_algal_cm2(2)/total_area_cm2 ;
lobinit = total_algal_cm2(3)/total_area_cm2 ;

ALGALREMOVAL=ALGALREMOVAL';

% @@@@@@@@@@ Below is orginal code @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

shortfallturf=0;
shortfalldict=0;
shortfalllob=0;
excessturf=0;
excessdict=0;
excesslob=0;

ALGALREMOVAL2=ALGALREMOVAL;% use to put in final values

if turfinit>ALGALREMOVAL(1,1); % enough turf available
    excessturf=turfinit-ALGALREMOVAL(1,1);
else
    shortfallturf=ALGALREMOVAL(1,1)-turfinit;
    ALGALREMOVAL2(1,1)=turfinit;
end

if dictinit>ALGALREMOVAL(1,2); % enough dict available
    excessdict=dictinit-ALGALREMOVAL(1,2);
else
    shortfalldict=ALGALREMOVAL(1,2)-dictinit;
    ALGALREMOVAL2(1,2)=dictinit;
end

if lobinit>ALGALREMOVAL(1,3); % enough lob available
    excesslob=lobinit-ALGALREMOVAL(1,3);
else
    shortfalllob=ALGALREMOVAL(1,3)-lobinit;
    ALGALREMOVAL2(1,3)=lobinit;
end

if shortfallturf>0 %need to look at next most preferred food
    if excessdict>0 % if excess dict available
        if excessdict>shortfallturf %if enough dict to take it
            ALGALREMOVAL2(1,2)=ALGALREMOVAL2(1,2)+shortfallturf;
            excessdict=excessdict-shortfallturf;
            shortfallturf=0;
        else %take what can and keep looking to lob
            ALGALREMOVAL2(1,2)=ALGALREMOVAL2(1,2)+excessdict;
            shortfallturf=shortfallturf-excessdict;
            excessdict=0;
        end
    end
    if excesslob>0 & shortfallturf>0
        if excesslob>shortfallturf %if enough lob to take
            ALGALREMOVAL2(1,3)=ALGALREMOVAL2(1,3)+shortfallturf;
            excesslob=excesslob-shortfallturf;
        else %shouldn't happen but if not enough lob
            ALGALREMOVAL2(1,3)=ALGALREMOVAL2(1,3)+(shortfallturf-excesslob);
            shortfallturf=shortfallturf-excesslob;
            excesslob=0;
        end
    end
end

if shortfalldict>0 %need to look at next most preferred food
    if excessturf>0 % if excess turf available
        if excessturf>shortfalldict %if enough turf to take it
            ALGALREMOVAL2(1,1)=ALGALREMOVAL2(1,1)+shortfalldict;
            excessturf=excessturf-shortfalldict;
            shortfalldict=0;
        else %take what can and keep looking to lob
            ALGALREMOVAL2(1,1)=ALGALREMOVAL2(1,1)+excessturf;
            shortfalldict=shortfalldict-excessturf;
            excessturf=0;
        end
    end
    if excesslob>0 & shortfalldict>0
        if excesslob>shortfalldict %if enough lob to take
            ALGALREMOVAL2(1,3)=ALGALREMOVAL2(1,3)+shortfalldict;
            excesslob=excesslob-shortfalldict;
        else %shouldn't happen but if not enough lob
            ALGALREMOVAL2(1,3)=ALGALREMOVAL2(1,3)+(shortfalldict-excesslob);
            shortfalldict=shortfalldict-excesslob;
            excesslob=0;
        end
    end
end

if shortfalllob>0 %need to look at next most preferred food
    if excessturf>0 % if excess turf available
        if excessturf>shortfalllob %if enough turf to take it
            ALGALREMOVAL2(1,1)=ALGALREMOVAL2(1,1)+shortfalllob;
            excessturf=excessturf-shortfalllob;
            shortfalllob=0;
        else %take what can and keep looking to dict
            ALGALREMOVAL2(1,1)=ALGALREMOVAL2(1,1)+excessturf;
            shortfalllob=shortfalllob-excessturf;
            excessturf=0;
        end
    end
    if excessdict>0 & shortfalllob>0
        if excessdict>shortfalllob %if enough dict to take
            ALGALREMOVAL2(1,2)=ALGALREMOVAL2(1,2)+shortfalllob;
            excessdict=excessdict-shortfalllob;
        else %shouldn't happen but if not enough lob
            ALGALREMOVAL2(1,2)=ALGALREMOVAL2(1,2)+(shortfalllob-excessdict);
            shortfalllob=shortfalllob-excessdict;
            excessdict=0;
        end
    end
end
