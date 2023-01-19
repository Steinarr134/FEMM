% Solution to Example 1.7 in Rao: The finite element method in Eng., 5.ed. 
% Indriði Sævar Ríkharðsson / Jonas Thor Snaebjornsson 2018
% For four types of analyses 

clc; 
close all; 
clear
Texti = char('Pipe pressure','Flow in pipes','Flow in ends',...
             'Temperature','Heat flow','Total heat flow in nodes',....
             'Voltage','Current in resistors','Total currents in nodes',...
             'Displacements','Forces in springs','Total force in nodes');
    
%Type of analyses 1=pipe flow; 2=Thermal; 3=Electric; 4=springs; 
Type = 1
switch Type
    case 1
        d = [75, 50, 38, 50, 65, 25, 50, 75, 65, 75]/1000;        % Diameter in m
        L = [125, 175, 100, 150, 200, 125, 75, 50, 150, 100];   % Length in m
        rho= 0.8545 ;              % Density in g/cm^3
        rho = rho*(1/1000)/(1/(100^3));
        mu = 206.89/1000;                % Viscosity in Pa*s[kg/(ms)]
        disp('Pressure vector bar')
        P = [1.5; 1.4; 0; 0; 0; 1.3; 0; 0; 1; 1.2];   %  Pressure in bar
        P = P*10^5;   % Pressure in Pa
        disp('Flow vector')
        Q = [0;0;0;0;0;0;0;0;0;0];
    case 2
        kt = [0,0,0,0,0,0,0,0,0,0]  % Heat conductivity of elements;
        t  =  [1,1,1,1,1,1,1,1,1,1] % Length(thicness) of elements;
        % Temperature vector
        disp('Temperature vector')
        T = [1,1,1,1,1,1,1,1,1,1];
        % Heat flow
        disp('Heat flow [W]')
        F = [0,0,0,0,0,0,0,0,0,0];
        P = T';
        Q =F';
    case 3
        % Voltage vector
        disp('Voltage vector')
        V = [0,0,0,0,0,0,0,0,0,0];
        % Current vector
        disp('Current vector')
        I = [0,0,0,0,0,0,0,0,0,0];
        P = V';
        Q = I';
        R = [1,1,1,1,1,1,1,1,1,1]; %Resistance
    case 4
        % Spring constants
        kg = [1,1,1,1,1,1,1,1,1];
        % Displacemets for known values
        disp('Displacement vector')
        U = [0,0,0,0,0,0,0,0,0,0];
        % Force vector known values
        disp('Force vector')
        F = [0,0,0,0,0,0,0,0,0,0];
        P = U';
        Q = F';
end

%%
disp('Number of nodes')
np = size(P,1)              % Calculate number of nodes.
disp('Element connection matrix')
              %    1    2    3    4    5    6    7    8    9    10
connect_matrix = [1,4; 2,3; 3,4; 4,8; 8,7; 3,5; 5,7; 5,6; 7,9; 8,10];

% fra maxim
% connect_matrix = [1,4; 2,3; 3,5; 4,3; 4,8; 5,7; 6,5; 7,8; 8,10; 9,7];
disp('Number of elements')
ne =size(connect_matrix,1)  % Calculates number of elements
disp('matrix for known deegres of freedom 0 = known , 1 = unknown')
ir=[0,0,1,1,1,0,1,1,0,0]
% Zero Global matrix
K = zeros(np);

% Setting up 3D matrices and vectors for element values
K0=[1,-1;-1,1];
Ke =cat(3,K0);      % Concatenate arrays
fe =cat(3,[0;0]);   %
disp('Calculation of element constants')

%%
% Element system matrices defined

% what is k?? why is it defined on the fly in this loop and then only
% used inside this loop? This doesn't make any sense. -SHH
for i=1:ne
    if Type == 1 
        k(i) = (d(i)^4/L(i))*(pi/(128*mu));
    elseif Type == 2
        k(i) = kt(i)/t(i);
    elseif Type == 3
        k(i) = 1/R(i);
    elseif Type == 4
        k(i) = kg(i);
    end
    Ke(:,:,i)=K0*k(i);
end
k

% Global system matrix assembled
for nel=1:ne
    i=connect_matrix(nel,1);
    j=connect_matrix(nel,2);
    K(i,i)=K(i,i)+Ke(1,1,nel);
    K(j,j)=K(j,j)+Ke(2,2,nel);
    K(i,j)=K(i,j)+Ke(1,2,nel);
    K(j,i)=K(j,i)+Ke(2,1,nel);
end

% Let x be vektor that includes unknowns
% alter Q vektor for known values
y = Q - K*P;
%%
% implement boundary values zero out corresponding lines and columns
KR=K
for i=1:np
    for j=1:np
        if ir(i)==0
            KR(i,j)=0;
            KR(j,i)=0;
            KR(i,i)=1;
            y(i)=0;
        end
    end
end
KR

disp('Find unknown values')
x = KR\y;

% Update solution vector with known values
P = P+x;

%Find element solution values
for i=1:ne
    p(i,1)=P(connect_matrix(i,1));
    p(i,2)=P(connect_matrix(i,2));
end

%%
disp('Flow and velocity in pipes')
%Calculate element values qe = [ke]*pe
for i =1:ne
    pe=[p(i,1);p(i,2)];
    qe(:,:,i)=Ke(:,:,i)*pe;
    v(i)=qe(1,i)/(pi*d(i)^2/4);
    Re(i)= abs(v(i))*d(i)*rho/mu; % Reynolds number
end
disp(Texti(Type*3-2,:))
p
disp(Texti(Type*3-1,:))
qe
disp(Texti(Type*3,:))
Q = K*P*1000000
disp('sum of Flow in ends, should be 0')
sum(Q)
if Type == 1
    disp('Velocity in pipes[cm/s]')
     v'.*100
    disp('Reynoldsnumber in pipes')
    Re'
end

