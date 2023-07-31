function [td, md, P,Pout,Dout] = calculate_dispersion_slim_clean(C,edge,A,DADx,DA2Dx,U,rdt,rdx,k,inflowlist,outflowlist)
%
% CALCULATE_DISPERSION - calculates the dispersion in a network
%   The dispersion is calculated from appropriate Taylor dispersion formula
%
%   Input:  C       - connectivity matrix of network
%           edge    - edges building the network
%           A       - contraction pattern of all edges, radius check if for
%                      scale
%           DADx    - space derivative of radius
%           DA2Dx   - space derivative of area
%           initedge- initialization edge
%           rdt      - incremental time step ADJUST IF 
%           rdx      - incremental space step
%           U       - longitudinal flow velocity in all edges
%           NL      - discretization of edge
%           k       - molecular diffusivity
%           inflowlist - contains the edges which are a dead end (left side)
%           outflowlist - contains the edges which are a dead end (right side)
%   Output: td       - output matrix for CN
%           md       - output matrix for CN
%
%           P, Pout, Dout - variables to check during debugging


%% Parameters
                     
ak = 1/k;

knc=k*ones(size(DADx));

aknc=1./knc;


%% Initialize variables used throughout the function
E = size(edge,1);                 % number of edges
T = size(A,1);                    % number of time points
TL = size(A,2);                   % total number of space points
BC = sum(sum(C(edge(:,2),:))) + sum(sum(C(edge(:,1),:))) ...
   + 2 * E - 2*sum(sum(C)==3);      % total number of boundary conditions
dim = size(A);


%% Calculate transport and diffusion term

% Here all terms are included

% First derivative coefficient
  P = (0.25 * rdt/rdx) * ( -U .* ones(dim) ...
      -aknc .* (U .* U .* A .* DADx)/24.0 - 2 * knc .* (DADx ./A));
  


% Constant term (for absorption)
Anc =0* 0.5* rdt* (2*U.*(DADx ./A));

% Second derivative coefficient
D = (0.5 * rdt/(rdx)^2) * (knc .* ones(dim) + aknc .* (U .* U .* A .* A)/48.0);


%for debugging

Pout=-P/(0.25 * rdt/rdx);
Dout=D/(0.5 * rdt/(rdx)^2);
Ancmat=[Pout(end),Dout(end)];


% Calculate matrix entries

% Use these to sum over all contribution from different edges/ create
% matrixes
or = zeros(dim);
ol = zeros(dim);
dir = zeros(dim);
dil = zeros(dim);
cn = 0;
tde = zeros(T,BC+3*TL-6*E);
mde = zeros(T,BC+3*TL-6*E);
row = zeros(1,BC+3*TL-6*E);
column = zeros(1,BC+3*TL-6*E);
pmatr=[];
pmatl=[];
for i = 1:1:E
    % Compute tridiagonal matrix elements
    

    Ah = 0.5 * (Anc(:,edge(i,3)+1:edge(i,4)-1)+Anc(:,edge(i,3)+2:edge(i,4)));
    Ph = 0.5 * (P(:,edge(i,3)+1:edge(i,4)-1) + P(:,edge(i,3)+2:edge(i,4)));
    Dh = 0.5 * (D(:,edge(i,3)+1:edge(i,4)-1) + D(:,edge(i,3)+2:edge(i,4)));
    or(:,edge(i,3)+1:edge(i,4)-1) = Dh - Ph; 
    ol(:,edge(i,3)+1:edge(i,4)-1) = Dh + Ph; 
    dir(:,edge(i,3)+1:edge(i,4)) = [1 + Ph(:,1) - Dh(:,1)-Ah(:,1) ...
        ones(size(A,1),size(Ph,2)-1) + (Ph(:,2:end) - Ph(:,1:end-1)) ...
        - (Dh(:,1:end-1) + Dh(:,2:end))-Ah(:,(1:end-1)) 1 - Ph(:,end) - Dh(:,end)-Ah(:,end)];
    dil(:,edge(i,3)+1:edge(i,4)) = [1 - Ph(:,1) + Dh(:,1)+Ah(:,1)...
        ones(size(A,1),size(Ph,2)-1) - (Ph(:,2:end) - Ph(:,1:end-1)) ...
        + (Dh(:,1:end-1) + Dh(:,2:end))+Ah(:,2:end) 1 + Ph(:,end) + Dh(:,end)+Ah(:,end)];

    len = edge(i,4) - edge(i,3) - 2;
    
    % diag 0
    row(cn+1:cn+len) = edge(i,3)+2:1:edge(i,4)-1;
    column(cn+1:cn+len) = edge(i,3)+2:1:edge(i,4)-1;
    tde(:,cn+1:cn+len) = dil(:,edge(i,3)+2:1:edge(i,4)-1);
    mde(:,cn+1:cn+len) = dir(:,edge(i,3)+2:1:edge(i,4)-1);
    cn = cn + len;
    
    % diag -1
    row(cn+1:cn+len) = edge(i,3)+2:1:edge(i,4)-1;
    column(cn+1:cn+len) = edge(i,3)+1:1:edge(i,4)-2;
    tde(:,cn+1:cn+len) = -or(:,edge(i,3)+1:1:edge(i,4)-2);
    mde(:,cn+1:cn+len) = or(:,edge(i,3)+1:1:edge(i,4)-2);
    cn = cn + len;
    
    % diag +1
    row(cn+1:cn+len) = edge(i,3)+2:1:edge(i,4)-1;
    column(cn+1:cn+len) = edge(i,3)+3:1:edge(i,4);
    tde(:,cn+1:cn+len) = -ol(:,edge(i,3)+2:1:edge(i,4)-1);
    mde(:,cn+1:cn+len) = ol(:,edge(i,3)+2:1:edge(i,4)-1);
    cn = cn + len;
    
    
    % Calculate entries for bifurcation points
    % left boundary bifurcation point
    if sum(C(edge(i,1),:)) ~= 1
        brrin = find(edge(i,1) == edge(:,2))';
        brl1 = find(edge(i,1) == edge(1:i-1,1));         
        brl2 = find(edge(i,1) == edge(i+1:end,1));
        if ~isempty(brl2)
            brl2 = (i) * ones(size(brl2,1),1) + brl2;

        end
        brlin = [brl1' brl2'];
        brl = edge(brlin,3)' + ones(size(brlin));
        brr = edge(brrin,4)';

        lbrin = [brrin];
        lbr2 = [brl brr];
        lbr = [brr];
        mdi = ones(T,1);%ones
        tdi = ones(T,1);%ones
        nn = size(lbr,2);
        nn2 = size(lbr2,2);
        %setting the values for neighbors
       


        pfori=0;
        aBoundi=0;
        sumA=0;

        
        for k = 1:1:nn
            
            porgi = 0.5 * (P(:,edge(i,3)+2) + P(:,edge(i,3)+1));
            dadxBoundi = -1*(A(:,lbr(k))^2-A(:,edge(i,3)+2)^2)/(2*rdx);
            dadxBoundk = -1*(A(:,lbr(k)-1)^2-A(:,edge(i,3)+1)^2)/(2*rdx);
            kncin=knc(:,lbr(k))+U(:,lbr(k)).^2.*A(:,lbr(k)).^2./(48*knc(:,lbr(k)));
            
            if nn==1
                pBoundi = (0.25 * rdt/rdx) * ( -U(:,edge(i,3)+1) + 2*  kncin * (A(:,lbr(k))^2-A(:,edge(i,3)+1)^2)/(A(:,lbr(k))^2+A(:,edge(i,3)+1)^2));
                pBoundk = (0.25 * rdt/rdx) * ( -U(:,edge(i,3)+1)*A(:,edge(i,3)+1)^2/A(:,lbr(k))^2 + 2*  kncin * (A(:,lbr(k))^2-A(:,edge(i,3)+1)^2)/(A(:,lbr(k))^2+A(:,edge(i,3)+1)^2));
            else
                pBoundk = (0.25 * rdt/rdx) * ( -U(:,lbr(k)) + 2*  kncin * (A(:,lbr(k))^2-A(:,edge(i,3)+1)^2)/(A(:,lbr(k))^2+A(:,edge(i,3)+1)^2));
                pBoundi = (0.25 * rdt/rdx) * ( -U(:,lbr(k))*A(:,lbr(k))^2/A(:,edge(i,3)+1)^2 + 2*  kncin * (A(:,lbr(k))^2-A(:,edge(i,3)+1)^2)/(A(:,lbr(k))^2+A(:,edge(i,3)+1)^2));
            end
            
            p= 1*(0.5 *(pBoundk));


            
            d = (0.5 * rdt/(rdx)^2)*kncin;
            
            aBoundi = aBoundi + 0.5* 1/4/rdx* ( (1+ (A(:,lbr(k))^2)/(A(:,edge(i,3)+1)^2)) * (2*U(:,lbr(k))) );
            aBoundk=0;
            
                tdi = tdi + pBoundi + d +aBoundk;
                tde(:,cn + k) = pBoundk - d; 
                mdi = mdi - pBoundi - d -aBoundk;
                mde(:,cn + k) = - pBoundk + d; 
            column(cn+k) = lbr(k);
        end
        dadxBoundifull =  2*0.5*(sumA-A(:,edge(i,3)+2)^2)/(2*rdx);
        
        %setting the values for end point itself
        p = 0.5 * (P(:,edge(i,3)+2) + P(:,edge(i,3)+1));
        
        
        aBoundi=0;
        
        
        
        d = 0.5 * (D(:,edge(i,3)+2) + D(:,edge(i,3)+1));
        tde(:,cn + nn + 1) = - p - d; 
        tde(:,cn + nn + 2) = tdi - p + d +aBoundi;
        mde(:,cn + nn + 1) = p + d; 
        mde(:,cn + nn + 2) = mdi + p - d -aBoundi;
        column(cn+nn+1) = edge(i,3)+2;
        column(cn+nn+2) = edge(i,3)+1;
        row(cn+1:cn+nn+2) = (edge(i,3)+1) * ones(1,nn+2);
        cn = cn + nn + 2;


    else
        % reflecting boundary
        if ismember(i,inflowlist)==0
            
            
            if edge(i,3)~=0
                column(cn+1:cn+3) = [edge(i,3) edge(i,3)+1 edge(i,3)+2];
                row(cn+1:cn+3) = (edge(i,3)+1) * ones(1,3);
                tde(:,cn+1:cn+3) = [-or(:,edge(i,3)) dil(:,edge(i,3)+1) -ol(:,edge(i,3)+1)];
                mde(:,cn+1:cn+3) = [or(:,edge(i,3)) dir(:,edge(i,3)+1) ol(:,edge(i,3)+1)];
                cn = cn + 3;
            else
                column(cn+1:cn+3) = [edge(i,3)+1 edge(i,3)+1 edge(i,3)+2];
                row(cn+1:cn+3) = (edge(i,3)+1) * ones(1,3);
                tde(:,cn+1:cn+3) = [-or(:,edge(i,3)+1) dil(:,edge(i,3)+1) -ol(:,edge(i,3)+1)];
                mde(:,cn+1:cn+3) = [or(:,edge(i,3)+1) dir(:,edge(i,3)+1) ol(:,edge(i,3)+1)];
                cn = cn + 3;
                
            end
        else
            
            column(cn+1:cn+2) = [edge(i,3)+1 edge(i,3)+2];
            row(cn+1:cn+2) = (edge(i,3)+1) * ones(1,2);
            tde(:,cn+1:cn+2) = [dil(:,edge(i,3)+1)+2*Ph(:,1) -ol(:,edge(i,3)+1)];
            mde(:,cn+1:cn+2) = [dir(:,edge(i,3)+1)-2*Ph(:,1) ol(:,edge(i,3)+1)];
            cn = cn + 2;
        end
    end
    
    %set right boundary
    %bifurcation point
    if sum(C(edge(i,2),:)) ~= 1
        brlin = find(edge(i,2) == edge(:,1))';
        brr1 = find(edge(i,2) == edge(1:i-1,2));         
        brr2 = find(edge(i,2) == edge(i+1:end,2));
        if ~isempty(brr2)
            brr2 = (i) * ones(size(brr2,1),1) + brr2;

        end
        brrin = [brr1' brr2'];
        brr = edge(brrin,4)';
        brl = edge(brlin,3)' + ones(size(brlin));

        brin = [brlin];
        rbr2 = [brl brr];
        rbr = [brl];
        mdi = ones(T,1);%ones
        tdi = ones(T,1);%ones
        nn = size(rbr,2);
        nn2 = size(rbr2,2);
        %setting the values for neighbors

        pfori=0;
        aBoundi=0;
        sumA =0;

        for k = 1:1:nn

            porgi = 0.5 * (P(:,edge(i,4)-1) + P(:,edge(i,4)));
            dadxBoundi = -1*(A(:,edge(i,4)-1)^2-A(:,rbr(k))^2)/(2*rdx);
            dadxBoundk = -1*(A(:,edge(i,4))^2-A(:,rbr(k)+1)^2)/(2*rdx);
            kncin=knc(:,edge(i,4))+U(:,edge(i,4)).^2.*A(:,edge(i,4)).^2./(48*knc(:,edge(i,4)));
            if nn==1
                pBoundi = (0.25 * rdt/rdx) * ( -U(:,edge(i,4)) - 2*  kncin * (A(:,rbr(k))^2-A(:,edge(i,4))^2)/(A(:,rbr(k))^2+A(:,edge(i,4))^2));
                pBoundk = (0.25 * rdt/rdx) * ( -U(:,edge(i,4))*A(:,edge(i,4))^2/A(:,rbr(k))^2 - 2*  kncin * (A(:,rbr(k))^2-A(:,edge(i,4))^2)/(A(:,rbr(k))^2+A(:,edge(i,4))^2));
            else
                pBoundi = (0.25 * rdt/rdx) * ( -U(:,rbr(k))*A(:,rbr(k))^2/A(:,edge(i,4))^2 - 2*  kncin * (A(:,rbr(k))^2-A(:,edge(i,4))^2)/(A(:,rbr(k))^2+A(:,edge(i,4))^2));
                pBoundk = (0.25 * rdt/rdx) * ( -U(:,rbr(k)) - 2*  kncin * (A(:,rbr(k))^2-A(:,edge(i,4))^2)/(A(:,rbr(k))^2+A(:,edge(i,4))^2));
            end
            
            p= 1*(0.5 *(pBoundi+pBoundk));
           
            pfori=pfori+p;
            pmatr=[pmatr;edge(i,4) rbr(k) k p pBoundi pBoundk -U(:,edge(i,4))];
            
            d = (0.5 * rdt/(rdx)^2)*kncin;
            aBoundi = aBoundi - 0.5*1/4/rdx*( (1+(A(:,rbr(k))^2/(A(:,edge(i,4))^2))) * (2*U(:,edge(i,4))) );

            aBoundk=0;

                tdi = tdi - pBoundi + d+aBoundk;
                tde(:,cn + k) = - (pBoundk + d); 
                mdi = mdi + pBoundi - d - aBoundk;
                mde(:,cn + k) = pBoundk + d; 
           
            column(cn+k) = rbr(k);
        end
        %setting the values for end point itself
        dadxBoundifull = 2*0.5*(A(:,edge(i,4)-1)^2-sumA)/(2*rdx);

       
        p = 0.5 * (P(:,edge(i,4)-1) + P(:,edge(i,4)));
      

        
         aBoundi=0;
        
        d = 0.5 * (D(:,edge(i,4)-1) + D(:,edge(i,4)));
        tde(:,cn + nn + 1) = p - d; 
        tde(:,cn + nn + 2) = tdi + p + d +aBoundi;
        mde(:,cn + nn + 1) = - p + d; 
        mde(:,cn + nn + 2) = mdi - p - d -aBoundi;
        column(cn+nn+1) = edge(i,4)-1;
        column(cn+nn+2) = edge(i,4);
        row(cn+1:cn+nn+2) = edge(i,4) * ones(1,nn+2);
        cn = cn + nn + 2;


    else
        % reflecting boundary
        if ismember(i,outflowlist)==0
            [outflowlist]
            column(cn+1:cn+3) = [edge(i,4)-1 edge(i,4) edge(i,4)+1];
            row(cn+1:cn+3) = edge(i,4) * ones(1,3);
            tde(:,cn+1:cn+3) = [-or(:,edge(i,4)-1) dil(:,edge(i,4)) -ol(:,edge(i,4))];
            mde(:,cn+1:cn+3) = [or(:,edge(i,4)-1) dir(:,edge(i,4)) ol(:,edge(i,4))];
            cn = cn + 3;
        else
            %i
            column(cn+1:cn+2) = [edge(i,4)-1 edge(i,4)];
            row(cn+1:cn+2) = edge(i,4) * ones(1,2);
            tde(:,cn+1:cn+2) = [-or(:,edge(i,4)-1) dil(:,edge(i,4))];
            mde(:,cn+1:cn+2) = [or(:,edge(i,4)-1) dir(:,edge(i,4))];
            cn = cn + 2;
        end
    end
end



%% Sparse matrices


size(row);
size(column);
max(column);
find(column==max(column));
find(tde==0);
row(find(tde==0))=[];
column(find(tde==0))=[];
mde(find(tde==0))=[];
tde(find(tde==0))=[];
row(find(column==max(column)));
max(row);
size(TL);
TL;
td = sparse(row, column,tde(:),TL,TL);
md = sparse(row, column,mde(:),TL,TL);
Ancmat= md;

    

end
