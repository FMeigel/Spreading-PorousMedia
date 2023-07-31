%% This driver files computes the spatial profiles as a solute disperses by
% advection and diffusion in porous medium

%NOTE: For this code to run succesfully you need a porous medium to
%simulate. Either load the porous medium in the workspace before running
%the code or load a .mat file

numberfiles=1;

pictmatnames=["../../Vs/Rot/CFP_00_V1.mat"]; %% Fill in correct path




Tmat=[1030];

savematnames=["../../Vs/CFP_test.mat"]; % %Fill in correct path





for countmat=1:1:numberfiles

load(pictmatnames(countmat))
%% Set the values to scale the flow velocity and the diffusion constant


ufac=0.1*1.75/2/8*3;
radiusfac=1.0;
kavlue=0.5*1.5*2;
kavlue=5.0;

%% Initialize your code parameters

isexperimental=true;
% For false you need a porous medium from createNetwork.m
% For true you need a porous medium from a 2d image

iscut=false;
% For this false you calculate spatial profiles
% For this true you calculate the first passage time profiles!


% Define parameters for dispersion code
T=Tmat(countmat);

dt=1;
dx=1;

%% Rename variables for later use
C=full(pore_C_sparse);


edge=network_edges;
vertices=network_vertices;
diameter=diameter_pore;


%% Sorting into bins

% A sorting into bins is performed to have comprehensive output and to
% define inflow nodes, outflow nodes, and boundary nodes


if isexperimental == true
   
        
        
    % correct negative velocities. 

    orderededge=edge;
    U_avg=zeros(1,size(edge,1));
    Dia_avg=zeros(1,size(edge,1));
    for countedge=1:1:size(edge,1)
        U_avg(countedge)=mean(U_pore{countedge});
        Dia_avg(countedge)=mean(diameter{countedge});
    end
    for eleminedge=1:1:size(edge,1)
        if U_avg(eleminedge)<0

            firstent=orderededge(eleminedge,1);
            orderededge(eleminedge,1)=orderededge(eleminedge,2);
            orderededge(eleminedge,2)=firstent;
            mean(U_pore{eleminedge});
            U_pore{eleminedge}=-U_pore{eleminedge};
            mean(U_pore{eleminedge});

        end
    end

    U_avg=zeros(1,size(edge,1));
    Dia_avg=zeros(1,size(edge,1));
    for countedge=1:1:size(edge,1)
        U_avg(countedge)=mean(U_pore{countedge});
        Dia_avg(countedge)=mean(diameter{countedge});
    end

    for countegde=1:1:size(edge,1)
        if U_avg(countegde)<0.5e-13
           U_avg(countegde)=0;
           for elemcount=1:1:length( U_pore{countegde})
               U_pore{countegde}(elemcount)=0;
           end
        end
    end

    U_avg=zeros(1,size(edge,1));
    Dia_avg=zeros(1,size(edge,1));
    for countedge=1:1:size(edge,1)
        U_avg(countedge)=mean(U_pore{countedge});
        Dia_avg(countedge)=mean(diameter{countedge});
    end
    for repeater=1:1:4
        for eleminedge=1:1:size(edge,1)
            if U_avg(eleminedge)==0
                firstent=orderededge(eleminedge,1);
                if ismember(firstent,orderededge(:,2))
                   [~,index]=find(orderededge(:,2)==firstent);
                   %[size(index)]
                   if size(index)==1
                       [U_avg(index)];
                       if U_avg(index)==0

                           orderededge(eleminedge,1)=orderededge(eleminedge,2);
                           orderededge(eleminedge,2)=firstent;
                       end
                   elseif size(index)==2
                       if U_avg(index(1))==0 & (U_avg(index(2))==0)
                           [eleminedge]
                           orderededge(eleminedge,1)=orderededge(eleminedge,2);
                           orderededge(eleminedge,2)=firstent;
                       end
                   end

                end
                firstent=orderededge(eleminedge,1);
                if ~ismember(firstent,orderededge(:,2))
                    orderededge(eleminedge,1)=orderededge(eleminedge,2);
                    orderededge(eleminedge,2)=firstent;
                end
            end
        end
    end


    % Define outflow nodes

    outflowlist=[];
    outflowInd=setdiff(orderededge(:,2),orderededge(:,1));
    for eleminoutflow=outflowInd'
        outflowedge=find(orderededge(:,2)==eleminoutflow);
        for eleminoutflowedge=outflowedge'
            outflowlist=[outflowlist,eleminoutflowedge];
        end
    end

    %Define inflow nodes

    inflowlist=[];
    inflowInd=setdiff(orderededge(:,1),orderededge(:,2));
    for elemininflow=inflowInd'
        inflowedge=find(orderededge(:,1)==elemininflow);
        for elemininflowedge=inflowedge'
            inflowlist=[inflowlist,elemininflowedge];
        end
    end

    %Sort out wrongly identified inflow nodes
    inflowlistcell=cell(2,1);
    outflowlistcell=cell(2,1);
    inflowlistupdate=inflowlist;
    outflowlistupdate=outflowlist;
    indexlistinflow=[];
    indexlistoutflow=[];

    %Discard wrongly identified inflow nodes
    for elemininflow=inflowlist
        if U_avg(elemininflow)<0.0001
            indexlistinflow=[indexlistinflow,find(inflowlist==elemininflow)];
        end        
    end
    for eleminoutflow=outflowlist
        if U_avg(eleminoutflow)==0
            indexlistoutflow=[indexlistoutflow,find(outflowlist==eleminoutflow)];
        end        
    end
    inflowlistupdate(indexlistinflow)=[];
    outflowlistupdate(indexlistoutflow)=[];

    % Define upper outflow nodes
    maxyvertex=max(network_vertices(:,2));
    findupvertices=find(network_vertices(:,2)==maxyvertex);
    edgesUplist=[];
    for countvertices=1:1:length(findupvertices)
        for countedges=1:1:length(orderededge)
            if orderededge(countedges,2)==findupvertices(countvertices)
                edgesUplist=[edgesUplist,countedges];
            end
            if orderededge(countedges,1)==findupvertices(countvertices)
                edgesUplist=[edgesUplist,countedges];
            end
        end
    end

    % In {1} is the inflow/outflow
    % In {2} is the boundary


    edgesinflowlist=[];
    for countedges=1:1:length(orderededge)
        for countvertex=1:1:length(inflowNodelist)
            if orderededge(countedges,1)==inflowNodelist(countvertex)
                edgesinflowlist=[edgesinflowlist,countedges];
            end
        end
    end
    edgesoutflowlist=[];
    for countedges=1:1:length(orderededge)
        for countvertex=1:1:length(outflowNodelist)
            if orderededge(countedges,2)==outflowNodelist(countvertex)
                edgesoutflowlist=[edgesoutflowlist,countedges];
            end
        end
    end

    for countegde=edgesinflowlist
       if length(U_pore{countegde})>20
           for elemcount=1:1:20
               U_pore{countegde}(elemcount)=U_pore{countegde}(21);
               diameter{countegde}(elemcount)=diameter{countegde}(21);
           end
       end
    end

    for countegde=edgesoutflowlist
       if length(U_pore{countegde})>20
           for elemcount=1:1:20
               U_pore{countegde}(end-elemcount+1)=U_pore{countegde}(end-20);
               diameter{countegde}(end-elemcount+1)=diameter{countegde}(end-20);
           end
       end
    end

    inflowlistdouble=[vertices(orderededge(edgesinflowlist,1),1)';edgesinflowlist];
    inflowlistdoublet=inflowlistdouble';
    inflowlistsorted=sortrows(inflowlistdoublet);
    edgesinflowlistc=inflowlistsorted(:,2)';

    for elemcounter=1:1:length(edgesinflowlistc)
        tubeval=edgesinflowlistc(elemcounter);
        moddef=1;
        multifac=1;
        modueval=mod(elemcounter,moddef);
        if  modueval~=0
            fullUcellvec=U_pore{tubeval};
            reducedUcellvec=fullUcellvec(1:40);
            enlargement=min(abs(moddef-modueval),abs(0-modueval));
            enlargedUvec=fullUcellvec;
            for enlarger=1:1:(multifac*enlargement)
                enlargedUvec=[enlargedUvec;reducedUcellvec];
            end
            U_pore{tubeval}=enlargedUvec;
            fullDcellvec=diameter{tubeval};
            reducedDcellvec=fullDcellvec(1:40);
            enlargedDvec=fullDcellvec;
            for enlarger=1:1:(multifac*enlargement)
                enlargedDvec=[enlargedDvec;reducedDcellvec];
            end
            diameter{tubeval}=enlargedDvec;
        end
    end

    inflowlistcell{1}=edgesinflowlist;
    outflowlistcell{1}=edgesoutflowlist;

    inflowlistcell{2}=inflowlist;
    outflowlistcell{2}=outflowlist;
        
    if iscut == true
        
        %This is giving the first passage time through a 2d image porous
        %medium
        
        % In {1} is the inflow/outflow
        % In {2} is the boundary
        
%         orderededge=edge;
%         outflowlistcell{1}=edgesindicesXsideUpNew;        
%         
%         inflowlistcell{2}=edgesindicesXsideDownNew;
%         outflowlistcell{2}=deadendedgesClean;
        
        %find center for initial edge to start
        
        pore_vertices=vertices;
        pore_edges=orderededge;
        
        maxvertx=max(pore_vertices(:,1));
        maxverty=max(pore_vertices(:,2));
        minvertx=min(pore_vertices(:,1));
        minverty=min(pore_vertices(:,2));
        xwindow=[minvertx+(maxvertx-minvertx)*1/3-100,minvertx+(maxvertx-minvertx)*1/3+100];
        ywindow=[minverty+(maxvertx-minvertx)/2-100,minverty+(maxvertx-minvertx)/2+100];
        initaledge=0;
        
        for searchedge=1:1:length(pore_edges)
           if pore_vertices(pore_edges(searchedge,1),1)>  xwindow(1) && pore_vertices(pore_edges(searchedge,1),1)<  xwindow(2)
               if pore_vertices(pore_edges(searchedge,1),2)>  ywindow(1) && pore_vertices(pore_edges(searchedge,1),2)<  ywindow(2)
                   if abs(pore_vertices(pore_edges(searchedge,1),2)-pore_vertices(pore_edges(searchedge,2),2)) >20
                       initaledge=searchedge;
                   end
               end
           end
        end
        
        if initaledge==0
            for searchedge=1:1:length(pore_edges)
               if pore_vertices(pore_edges(searchedge,2),1)>  xwindow(1) && pore_vertices(pore_edges(searchedge,2),1)<  xwindow(2)
                   if pore_vertices(pore_edges(searchedge,2),2)>  ywindow(1) && pore_vertices(pore_edges(searchedge,2),2)<  ywindow(2)
                       if abs(pore_vertices(pore_edges(searchedge,1),2)-pore_vertices(pore_edges(searchedge,2),2)) >20 
                           initaledge=searchedge;
                       end
                   end
               end
            end
        end
        
        inflowlistcell{1}=initaledge;
   
        
     end
    
else 
    % For artificial networks from createNetwork.m, only the option of flow
    % profiles is given in this code
    
    % Note that this part needs adjustment according to chosen artificial
    % network
    
    nbins=30;
    edgeslists=cell(1,nbins);
    for countbins=1:1:nbins
        edgeslists{countbins}=countbins;
    end
    
    % In {1} is the inflow/outflow
    % In {2} is the boundary
    
    orderededge=edge;
    inflowlistcell{1}=1;
    outflowlistcell{1}=nbins;
    inflowlistcell{2}=1;
    outflowlistcell{2}=nbins;
    %outflowlistcell{2}=[nbins:59];
    network_edges=pore_edges;
    network_vertices=pore_vertices;
        
end
%%



%% Calculate Spread
% !! Initial Conditions inside function
tic
[S,S_avg,S_all,saveedge]=calculate_transport_clean2(C,orderededge,diameter,U_pore,T,dt,dx,inflowlistcell,outflowlistcell,ufac,kavlue,radiusfac,iscut); %5000
toc


%% Check correct normalisation of U
U_avg=zeros(1,size(edge,1));
Dia_avg=zeros(1,size(edge,1));
for countedge=1:1:size(edge,1)
    U_avg(countedge)=mean(U_pore{countedge});
    Dia_avg(countedge)=mean(diameter{countedge});
end
testprop=zeros(1,size(edge,1));
testprop(U_avg>mean(U_avg))=1;

%% Sparse out for plotting
dilutionfac=5;
dilutionsteps=1:(dilutionfac-1):T;
dilutionS=S_avg(dilutionsteps,:);

save(savematnames(countmat),'dilutionS','orderededge','vertices','U_pore','ufac','kavlue')

end
