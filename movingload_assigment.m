clc
clear 
close all
addpath("FEM Functions\")
addpath("Strutture\")
f_max=20;
n_modes=4;
%% load file
[file_i, xy, nnod, sizee, idb, ndof, incidenze, l, gamma, m, EA, EJ, posiz, nbeam, pr] = loadstructure;

disp('node degree of freedom')
[[1:size(idb,1)]' idb]
ndof

%% mass and stiffness
[M,K] = assem(incidenze,l,m,EA,EJ,gamma,idb);

MFF = M(1:ndof,1:ndof);
KFF = K(1:ndof,1:ndof);
MCF = M(ndof+1:end,1:ndof);
KCF = K(ndof+1:end,1:ndof);
MFC = M(1:ndof,ndof+1:end);
MCC = M(ndof+1:end,ndof+1:end);
KFC = K(1:ndof,ndof+1:end);
KCC = K(ndof+1:end,ndof+1:end); 

%% dumping
ab=[0.1 2.0e-4];
C = ab(1)*M + ab(2)*K; 
CFF = C(1:ndof,1:ndof);
CCF= C(ndof+1:end,1:ndof);
CCC = C(ndof+1:end,ndof+1:end);
CFC = C(1:ndof,ndof+1:end);

%% mode shapes and natural frequency

[modi, d] = eig(MFF \ KFF);
fre = sqrt(diag(d))/2/pi; 
scale_factor=10;
[Y, ind] = sort(fre);     % Sort in ascending order
ndis = sum(Y<f_max);
modi_ind = ind(1:ndis);
tipoEF= ones(1,length(EA));
for ii = 1:ndis
    figure
    diseg2(modi(:,modi_ind(ii)), scale_factor, incidenze, l, gamma, posiz, idb, xy);
end
%% modal approach

Phi = modi(:,modi_ind(1:n_modes)); 
Mmod = Phi'*MFF*Phi; 
Kmod = Phi'*KFF*Phi; 
Cmod = Phi'*CFF*Phi;  

%% moving load
Mass=100; %[Kg]
v_M=2; %[m/s]
dt=0.001;
dx=v_M*dt;
Fn_global=[];
l=l(1,5:11);
incid=incidenze(5:11,:);
for idx = 1:length(l)
    a = dx:dx:l(idx);
    b = l(idx) - a;
    Fn_local = Mass*9.81*[zeros(size(a));
            a.*b.*(b-a)/l(idx)^3 + b./l(idx);
            a.*b.^2/l(idx)^2;
            zeros(size(a));
            a.*b.*(a-b)/l(idx)^3 + a./l(idx);
            -a.^2.*b/l(idx)^2];
Temp=zeros(size(M,1),length(a));
Temp(incid(idx, :),:) = Fn_local; % poi cambiare incidenzecon incid se voglio farlo su meno nodi
clear Fn_local
Fn_global = [Fn_global Temp];
end

L_part1=size(Fn_global,2);
Fn_global(:,end+2/dt) = 0;
Fn_global(ndof+1:end, :) = [];

t_vect=(1:size(Fn_global, 2))*dt;
figure
plot(t_vect, Fn_global')
ylabel ( 'Nodal forces')
xlabel('Time [s]')
%%
% passaggio coordinate modali

Qn_global = Phi'*Fn_global;

figure
plot(t_vect,Qn_global);
grid on
xlabel('Q');
ylabel('Time [s]')
%%
x0=zeros(1,n_modes*2);
[t, X]=ode45(@(t,X) odefn(t,X,t_vect,Qn_global,Mmod,Kmod,Cmod), t_vect, x0);
TH = Phi*X(:,end/2+1:end)';
X_displ= TH(idb(12,2),:);
figure('Color','w') % sfondo bianco
plot(t, X_displ, 'b', 'LineWidth', 1.5), grid on
xlabel('Time [s]')
ylabel('Vertical displacement of node A [m]')
title('Displacement of node A along y over time due to a 1 kg moving load at 2 m/s')
%%
figure;

% Plot delle tre componenti
plot(t, Fn_global(idb(12,1), :), 'r', 'DisplayName', 'Force X','LineWidth',2); hold on;
plot(t, Fn_global(idb(12,2), :), 'g', 'DisplayName', 'Force Y','LineWidth',2);
plot(t, Fn_global(idb(12,3), :), 'b', 'DisplayName', 'Moment','LineWidth',2);

% Aggiunta del titolo e delle etichette
title('Force applied on node 12');
xlabel('Tempo [s]');
ylabel('Force / Moment [N / Nm]');

% Aggiunta della legenda
legend show;
grid on;
set(gcf, 'Color', 'w');




function dXdt = odefn(t, X, t_vect, Qn_global, Mmod, Kmod, Cmod)
    n = length(Mmod);          % numero di modi
    q = X(1:n);                % coordinate modali
    dq = X(n+1:end);           % velocità modali

    % Interpolazione della forza a tempo t
    Qn = interp1(t_vect, Qn_global', t, 'linear', 0)';  % forza interpolata

    % Equazione del moto nel dominio modale:
    ddq = Mmod \ (Qn - Cmod*dq - Kmod*q);

    % Sistema in forma di stato
    dXdt = [dq; ddq];
end
