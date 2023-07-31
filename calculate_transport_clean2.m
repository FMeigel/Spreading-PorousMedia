function [S_norm,S_avg,S_all,edge]=calculate_transport_clean(C,edges,diameter,U_pore,NT,dt,dx,initialedgescell,finaledgescell,Ufac,kvalue,radiusfac,iscut)
%% Calculate spread in a network using Taylor Dispersion
% Code is optimised for fast calculation of large networks
% Needed: Network: C-Connectivity, edges, vertices, U-velocities
% To calculate the network from the image use porous media code
%
% INPUT         C       : connectivity matrix of network; NxN matrix
%               edge    : edges of the network; Ex2 matrix
%               diameter: diameter of the edges at every point; 1xE cell
%               U_pore  : velocities of the edges at every point
%               NT      : Number of Timesteps
%               dt      : time discrimential
%               dx      : space discrimential
%               initialedgescell    : cell with the initial boundary lists
%               finaledgescell      : cell with final boundary lists
%               Ufac    : scaling of the velocity
%               kvalue  : diffusion coefficient
%               iscut   : spatial fill up or first passage time?
%
% Output        S_norm  : Last profile of all space points
%               S_avg   : For each time point the spatial porfiles averaged
%                         over all edges
%               S_all   : Total amount of solute for each time point. This
%                         can be understood as survival probability.
% 



%% 

%inflow and outflow
initialedges=initialedgescell{1};
finaledges=finaledgescell{1};

%boundary edges
initialedgescoef=initialedgescell{2};
finaledgescoef=finaledgescell{2};

%%
%[E V C nn ee es] = calc_connectivity(edge,C);
E = size(edges,1);

k =kvalue;

%phi = [];
%dphi = phi;

A= cell(1,E);
A0 = cell(1,E);               %Base radii
DADt = cell(1,E);             %time derivative of edge radius
DADx = cell(1,E);             %spatial derivative of radius
DUDx = cell(1,E);             %spatial derivative of velocity
DADDADx = cell(1,E);          %second derivative of edge radius
a0=1;
rdx=dx;
sumTotalpoints=0;
for countEdges=1:1:E
    A{countEdges}=diameter{countEdges}/2*radiusfac; %/2
    A0{countEdges}=a0*ones(size(A{countEdges}));
    DADt{countEdges}=zeros(size(A{countEdges}));
    DADx{countEdges} = [(A{countEdges}(2)-A{countEdges}(1))/rdx ((A{countEdges}(3:end)-A{countEdges}(1:end-2))/(2*rdx))' (A{countEdges}(end)-A{countEdges}(end-1))/rdx]';
    DUDx{countEdges} = [(U_pore{countEdges}(2)-U_pore{countEdges}(1))/rdx ((U_pore{countEdges}(3:end)-U_pore{countEdges}(1:end-2))/(2*rdx))' (U_pore{countEdges}(end)-U_pore{countEdges}(end-1))/rdx]';
    DADDADx{countEdges} = [2*(DADx{countEdges}(2)-DADx{countEdges}(1))/rdx ((DADx{countEdges}(3:end)-DADx{countEdges}(1:end-2))/(2*rdx))' 2*(DADx{countEdges}(end)-DADx{countEdges}(end-1))/rdx]';
    sumTotalpoints=sumTotalpoints+size(A{countEdges},1);
end


%%

%Go from cells to matrices


Anc = zeros(1,sumTotalpoints);
DADxnc = zeros(size(Anc));
DUDxnc = zeros(size(Anc));
DADDADxnc = zeros(size(Anc));
Unc = zeros(size(Anc));
edge = zeros(E,4);
cPosition=0;
for countEdges = 1:1:size(edges,1)
    sizeA=size(A{countEdges},1);
    Anc((cPosition+1):(cPosition+sizeA)) = A{countEdges};
    DADxnc((cPosition+1):(cPosition+sizeA)) = DADx{countEdges};
    DUDxnc((cPosition+1):(cPosition+sizeA)) = DUDx{countEdges};
    DADDADxnc((cPosition+1):(cPosition+sizeA)) = DADDADx{countEdges};
    Unc((cPosition+1):(cPosition+sizeA)) = U_pore{countEdges}*Ufac;
    edge(countEdges,:) = [edges(countEdges,:) cPosition (cPosition+sizeA)];
    cPosition=cPosition+sizeA;
end

%clear('A','DADx','U');
%%


%% Calculate dispersion

a0 = Anc(1,1);

 
% THIS IS THE MOST IMPORTANT LINE.
% Here you calculate the matrices for Crank Nicolson
[td, md,Ancmat,Pnc,Dnc] = calculate_dispersion_slim_W2(C,edge,Anc,DADxnc,DUDxnc,Unc,dt,dx,k,initialedgescoef,finaledgescoef);




%% Index version for open outflow boundary condition


S_old=zeros(1,sumTotalpoints);
inflowlist=[];

%Make boundary condition list


for countedges=1:1:size(initialedges,2)
    if iscut == true
        inflowlist=[inflowlist,edge(initialedges(countedges),3)+10];
    else
        inflowlist=[inflowlist,edge(initialedges(countedges),3)+1];
    end
end
outflowlist=[];
for countedges=1:1:size(finaledges,2)
    outflowlist=[outflowlist,edge(finaledges(countedges),4)];
end

specinflow=[];
for countsize=1:1:size(inflowlist,2)
    specinflow=[specinflow,(1.0)/(1.0/(A{initialedges(countsize)}(1).^2))];
end
specoutflow=[];
for countsize=1:1:size(outflowlist,2)
    specoutflow=[specoutflow,0];
end


% Specify the inflow in the inflow nodes

if iscut == false
    inflowlist
    S_old(inflowlist)=specinflow;
    specinflow
else
     S_old(inflowlist)=1;
end



S_avg=zeros(NT,E);
S_all=zeros(NT,1);
for countEdges=1:1:E
    S_avg(1,countEdges)=mean(S_old(1,(edge(countEdges,3)+1):edge(countEdges,4))./((Anc(1,(edge(countEdges,3)+1):edge(countEdges,4))).^2));
end
S_all(1)=sum(S_old(:));

%Use matrix prepocessing to increase speed


[Ltd,Utd,Ptd,Qtd,Dtd]=lu(td);



% Iteratively apply the code to compute time progression 

for j = 2:1:NT
    % Invert matrix to compute dispersion in next time step
    % Apply steady inflow
     if iscut==false
         S_old(1,inflowlist)=specinflow;
     end
     
     if rem(j,1000)==0
        [j/1000] 
     end
    
    rhs = md * S_old(1,:)';
    S_oldold=S_old;
    %ss = td\rhs;
    %ss= S_td*(R_td\(R_td'\(S_td'*rhs)));
    ss=Qtd*(Utd\(Ltd\(Ptd*(Dtd\rhs))));
    %ss=dtd\rhs;
    S_old = ss';
    
    % Apply steady inflow
    if iscut==false
        S_old(inflowlist)=specinflow;
    else
        % S_old(inflowlist)=1;
    end
  
    
    % Open outflow boundary condition 
    specoutflowII=[];
    for countsize=1:1:size(outflowlist,2)
        S_oldend=1/2*(1*S_old(1,(outflowlist(countsize)))+1*S_oldold(1,(outflowlist(countsize)))); 
        hSb=log(S_old(1,(outflowlist(countsize)-3))-S_old(1,(outflowlist(countsize)-2)));
         hSm=log((S_old(1,(outflowlist(countsize)-2))-S_old(1,(outflowlist(countsize)-1)))/(S_old(1,(outflowlist(countsize)-3))-S_old(1,(outflowlist(countsize)-2))));
         expfacDeri=exp(hSm*(3+0.0014-1*Pnc(outflowlist(countsize))*dt)+hSb);
             
         if hSm >=0
            expfacDeri=0;
            %j
            %'hsm>=0'
         end
         
         if isnan(hSm)
            expfacDeri=0;
            %j
            %'isnan'
         end
         
         if isnan(expfacDeri)
             expfacDeri=0;
         end
         if isreal(expfacDeri)==0
             expfacDeri=abs(0);
             %j
         end
         
         
         if Anc(outflowlist(countsize))==Anc(outflowlist(countsize)-1)
             outflowval5=dt/dx*(S_oldend.*Pnc(outflowlist(countsize))+1.0*(Dnc(outflowlist(countsize)))/dx.*(expfacDeri));
         else
             outflowval5=dt/dx*(S_oldend.*Pnc(outflowlist(countsize))+1.0*(Dnc(outflowlist(countsize)))/dx.*(expfacDeri))*(abs(Anc(outflowlist(countsize))-Anc(outflowlist(countsize)-1))./Anc(outflowlist(countsize)));
         end
         
        
                  

        specoutflowII=[specoutflowII,(S_old(1,(outflowlist(countsize)))-outflowval5)];
        specdata = [k, Dnc(outflowlist(countsize)),Unc(outflowlist(countsize)),Pnc(outflowlist(countsize)), (S_old(1,(outflowlist(countsize)))), (S_old(1,(outflowlist(countsize)-1)))-(S_old(1,(outflowlist(countsize)))), -1.5*S_old(1,(outflowlist(countsize)-3))+4*S_old(1,(outflowlist(countsize)-2))-2.5*S_old(1,(outflowlist(countsize)-1)),exp(hSm*(3+0.005+Unc(outflowlist(countsize))*dt)+hSb),exp(hSm*(3+0.005-Unc(outflowlist(countsize))*dt)+hSb),exp(hSm*(3+0.005+100*Unc(outflowlist(countsize))*dt)+hSb)];
    end
    
    if iscut==false
        S_old(1,outflowlist)=specoutflowII;
        %'yes'
    else
        S_old(1,outflowlist)=zeros(size(specoutflowII));
    end
    for countEdges=1:1:E
        S_avg(j,countEdges)=mean(S_old(1,(edge(countEdges,3)+1):edge(countEdges,4))./((Anc(1,(edge(countEdges,3)+1):edge(countEdges,4))).^2));
    end
    S_all(j)=sum(S_old(:));
end
S_norm=S_old./(Anc.^2);
%specdata
end
