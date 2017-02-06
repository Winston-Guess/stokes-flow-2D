clear all ; close all ; clc ;
npoints = 100;
% Set the number of elements
n_el = 16; n_np = 45;
dof = 2 * n_np;
f_ = 10;
height = 2;
Length = 10;              % length of bar
% Set the properties of the material and problem 
D_visc = 10;
w = 6;          % weight of element

ICA = [1 3 13 2 8 7 % Interconnectivity array
       1 13 11 7 12 6
       3 5 15 4 10 9
       3 15 13 9 14 8
       11 13 23 12 18 17
       11 23 21 17 22 16
       13 15 25 14 20 19
       13 25 23 19 24 18
       21 23 33 22 28 27
       21 33 31 27 32 26
       23 25 35 24 30 29
       23 35 33 29 34 28
       31 33 43 32 38 37
       31 43 41 37 42 36
       33 35 45 34 40 39
       33 45 43 39 44 38];

dofICA = [ICA*2-1,ICA*2];
   

GA1=[0  0   0 0  0 1.25 1.25 1.25 1.25 1.25 2.5 2.5 2.5 2.5 2.5 3.75 3.75 3.75;
     2  1.5 1 0.5 0 2    1.5  1.   0.5  0    1.8   1.5 1   0.5 0    2   1.5 1.];
 
GA2 =[3.7500  3.7500  5.0000  5.0000  5.0000  5.0000  5.0000  6.2500  6.2500;
      0.5000    0     2.0000  1.5000  1.0000  0.5000     0    2.0000  1.5000];

GA3 =[6.2500  6.2500  6.2500  7.5000  7.5000  7.5000  7.5000  7.5000  8.7500;
      1.0000  0.5000   0      2.0000  1.5000  1.0000  0.5000     0.2    2.0000];

GA4 =[8.7500  8.7500  8.7500  8.7500  10.0000  10.0000  10.0000  10.0000  10;
      1.5000  1.0000  0.5000    0     2.0000   1.5000   1.0000   0.5000   0];

GA = [GA1, GA2, GA3, GA4];

A_e=inline('.5*((xe2*ye3-xe3*ye2)-(xe1*ye3-xe3*ye1)+(xe1*ye2-xe2*ye1))'...
    ,'xe1','xe2','xe3','ye1','ye2','ye3');

psi1GP = [.1666666666 .6666666666 .1666666666];
psi2GP = [.1666666666 .1666666666 .6666666666];
WGP = [.1666666666 .1666666666 .1666666666];

elx = zeros(n_el,7);                         % element x-coords var
ely = zeros(n_el,7);                         % element y-coords var
elxd = zeros(n_el,7);                         % element x-coords var
elyd = zeros(n_el,7);                         % element y-coords var
pres = zeros(n_el,7);                         % element y-coords var

draw = [1 4 2 5 3 6 7 1];

K_DD = zeros(dof , dof);
G = zeros(dof , dof/2);
Nil = zeros(dof/2 , dof/2);
fv_Gamma = zeros(dof,1);

for e = 1:n_el
    K_DDe = zeros(12 , 12);
    Ge = zeros(12,6);
    for gpn = 1:3
        psi1 = psi1GP(gpn);
        psi2 = psi2GP(gpn);
        psi3 = 1-psi1-psi2;

        for i=1:6
            elx(e,i) = GA(1,ICA(e,i));
            ely(e,i) = GA(2,ICA(e,i));

            if i==1
                elx(e,7) = GA(1,ICA(e,i));
                ely(e,7) = GA(2,ICA(e,i));
            end
        end
        xy_e = [elx(e,1:6)',ely(e,1:6)'];

        N1 = psi1*(2*psi1-1);
        N2 = psi2*(2*psi2-1);
        N3 = psi3*(2*psi3-1);
        N4 = 4*psi1*psi2;
        N5 = 4*psi2*psi3;
        N6 = 4*psi1*psi3;               % Table 7.5 Pg 175 Belytschko
        
        N = [N1 N2 N3 N4 N5 N6];
        
%%%%%%%%%%
        GN_e = [4*psi1-1, 0, -3+4*psi1+4*psi2, 4*psi2, -4*psi2, 4-4*psi2-8*psi1;
                0, 4*psi2-1, -3+4*psi2+4*psi1, 4*psi1, 4-4*psi1-8*psi2, -4*psi1];
%%%%%%%%%%%check
        Je = GN_e*xy_e;

        BB = Je\GN_e;
        B_DDe = [BB(1,1),0,BB(1,2),0,BB(1,3),0,BB(1,4),0,BB(1,5),0,BB(1,6),0;
              0,BB(1,1),0,BB(1,2),0,BB(1,3),0,BB(1,4),0,BB(1,5),0,BB(1,6);
              BB(2,1),0,BB(2,2),0,BB(2,3),0,BB(2,4),0,BB(2,5),0,BB(2,6),0;
              0,BB(2,1),0,BB(2,2),0,BB(2,3),0,BB(2,4),0,BB(2,5),0,BB(2,6)];
          
        B_De = [BB(1,1) BB(2,1), BB(1,2) BB(2,2), BB(1,3) BB(2,3),...
            BB(1,4) BB(2,4), BB(1,5) BB(2,5), BB(1,6) BB(2,6)]; 

        K_DDe = K_DDe + D_visc*WGP(gpn)*(B_DDe)'*B_DDe*det(Je);
        Ge = Ge + WGP(gpn)*(B_De)'*N*det(Je);
    end
        
    K_DD(dofICA(e,1:12),dofICA(e,1:12)) = K_DD(dofICA(e,1:12),dofICA(e,1:12)) + ...
        K_DDe([1 3 5 7 9 11 2 4 6 8 10 12],[1 3 5 7 9 11 2 4 6 8 10 12]);
    
    G(dofICA(e,1:12),ICA(e,1:6)) = G(dofICA(e,1:12),ICA(e,1:6)) + ...
        Ge([1 3 5 7 9 11 2 4 6 8 10 12],[1 4 2 5 3 6]);

    hold on

end

f = fv_Gamma;
H = zeros(dof/2,1);

F = [f;H];
K = [K_DD G; G' Nil];

GT = G';

% Inlet Velocity (Dirichlet) BCs
% x-Direction
vx_nodes = [2 3 4];
vx_dof = 2*vx_nodes-1;
vx = [.4 .1 0];

for i = 1:max(size(vx_dof))
    F = F - vx(i)*K(:,vx_dof(i));
end
% disp(K)
% disp(F)

% y-Direction
vy_nodes = [2 3 4 42 43 44];
vy_dof = 2*vy_nodes;
vy = [0 0 0 0 0 0];

for i = 1:max(size(vy_dof))
    F = F - vy(i)*K(:,vy_dof(i));
end


% Walls (Dirichlet) BCs
noslip_nodes = [1 6 11 16 21 26 31 36 41,5 10 15 20 25 30 35 40 45,2 3 4];
noslip = [noslip_nodes*2-1,noslip_nodes*2,[42 43 44]*2];

for i = 1:max(size(noslip))
    F = F - 0*K(:,noslip(i));
end

noslip = sort(noslip);

for i = 1:max(size(noslip))
    K(noslip(i),:) = 0;
    K(:,noslip(i)) = 0;
    K(noslip(i),noslip(i)) = 1;
end

K_DD_BC = K(1:dof,1:dof);
G_BC = K(1:dof,dof+1:dof+dof/2);
GT_BC = K(dof+1:dof+dof/2,1:dof);
f_BC = F(1:dof);
H_BC = F(dof+1:dof+dof/2);

% % Outlet Pressure (Nuemann) BCs
% px_nodes = [41 42 43 44 45];
% px = [0 0 0 0 0];

% for i = 1:size(vx_dof)
%     H = H - vx(i)*GT(:,px_nodes(i));
% end

ph = (GT_BC*(K_DD_BC\G_BC))\(GT_BC*(K_DD_BC\f_BC)-H_BC);
vh = K_DD_BC\(f_BC-G_BC*ph);
d = vh;
d(noslip(1:max(size(noslip)))) = 0;
d(vx_dof(1:3)) = vx(1:3);

% Calculate unknown nodal tempretures using partitioning method
for el = 1:n_el
    plot(elx(el,draw), ely(el,draw),'Color','k')
    set(gcf, 'Visible', 'off')
    hold on
end
set(gcf, 'Visible', 'on')


for el = 1:n_el
    for i=1:6
        elxd(el,i) = d(ICA(el,i)*2-1);
        elyd(el,i) = d(ICA(el,i)*2);
        pres(el,i) = ph(ICA(el,i));
        if i==1
            elxd(el,7) = d(ICA(el,i)*2-1);
            elyd(el,7) = d(ICA(el,i)*2);
            pres(el,7) = ph(ICA(el,i));
        end
    end
end
% contourf(elx,ely,pres), hold on
quiver(elx,ely,elxd,elyd,'Color','b')

set(gcf, 'Visible', 'on')
% Create xlabel
xlabel({'x (m)'},'FontWeight','normal','FontSize',14);
% Create ylabel
ylabel({'y (m)'},'FontWeight','normal'...
    ,'FontSize',14);
