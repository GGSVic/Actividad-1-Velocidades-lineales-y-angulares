%Limpieza de pantalla
clear all
close all
clc

%Declaración de variables simbólicas (No tienen un valor específico)
syms th1(t) th2(t) th3(t) l1 l2 l3 t

%Configuración del robot, 0 para junta rotacional, 1 para junta prismática
RP = [0 0 0]; 

%Creamos el vector de coordeanadas articulares
Q = [th1 th2 th3]; 
disp('Coordenadas articualres'); 
pretty(Q); 

%Creamos el vector de velocidades articulares
Qp = diff(Q, t); 
disp('Velocidades articulares')
pretty(Qp); 

%Número de grado de libertad del robot
GDL = size(RP, 2); %Dimensión o número de columnas
GDL_str = num2str(GDL); 

%ARTICULACIÓN 1

%Posición de la junta 1 respecto a 0
P(:,:,1) = [l1*cos(th1); 
            l1*sin(th1); 
                     0];

%Posición de la junta 2 respecto a 1
P(:,:,2) = [l2*cos(th2) ; 
           l2*sin(th2); 
              0]; 

%Posición de la junta 3 respecto a 1
P(:,:,3) = [l3*cos(th3) ; 
           l3*sin(th3); 
              0]; 

%Matriz de rotación de la articulación 1 respecto a 0
R(:,:,1) = [cos(th1) -sin(th1) 0; 
            sin(th1) cos(th1)  0;
            0           0      1]; 

%Matriz de rotación de la articulación 2 respecto a 1
R(:,:,2) = [cos(th2) -sin(th2) 0; 
            sin(th2) cos(th2)  0; 
            0           0      1]; 

%Matriz de rotación de la articulación 2 respecto a 1
R(:,:,3) = [cos(th3) -sin(th3) 0; 
            sin(th3) cos(th3)  0; 
            0           0      1]; 

%Creamos un vector de ceros
Vector_Zeros = zeros(1, 3); 

%Inicializamos las matrices de transformación Homogénea locales
A(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);

%Inicializamos las matrices de transformación Homogenea globales
T(:,:,GDL) = simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]); 

%Inicializamos los vectores d eposiicón vistos desde el marco de referencia
%inercial
PO(:,:,GDL) = P(:,:,GDL); %...

%Inicializamos las matrices de rotación vistas desde el marco de referencia
%inercial
RO(:,:,GDL) = R(:,:,GDL);  


for i = 1:GDL
   i_str = num2str(i);

   %Declaración de las matrices de transformación local para cada 
   %articulación, concatenando cada uno de sus elementos:  
   disp(strcat('Matriz de Transformación local A', i_str));
   A(:,:,i) = simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);
   pretty(A(:,:,i));

   %Calculamos la matriz de transformación homogénea global multiplicando 
   %la matriz actual con las matrices de las articulaciones previas: 
   try
       T(:,:,i) = T(:,:,i-1)*A(:,:,i);
   catch
       T(:,:,i) = A(:,:,i); %Para la primer articulación, la matriz global
                            %es equivalente a la local 
   end

   disp(strcat('Matriz de Transformación global T', i_str));
   T(:,:, i) = simplify(T(:,:,i));
   pretty(T(:,:,i));

   %Obtenemos la matriz de rotación "RO" y el vector de translación de PO
   %de la matriz de transformación homogénea global T(:,:,GDL);

   RO(:,:,i) = T(1:3,1:3,i);
   PO(:,:,i) = T(1:3,4, i);
   pretty(RO(:,:,i));
   pretty(PO(:,:,i));

end


%Calculamos el jacobiano lineal y angular de forma analitica
%Calculamos jacobianos analiticos ( lineal y angular)
Jv_a(:,GDL)=PO(:,:,GDL);
Jw_a(:,GDL)=PO(:,:,GDL);

for k = 1:GDL
   if (((RP(k)==0))|(RP(k)==1))
   %Para las articulaciones rotacionales: 
       try
           Jv_a(:,k) = cross(RO(:,3,k-1), PO(:,:,GDL)-PO(:,:,k-1));
           Jw_a(:,k) = RO(:,3,k-1);
       catch
           Jv_a(:,k) = cross([0,0,1], PO(:,:,GDL)); 

           Jw_a(:,k) = [0,0,1]; %No hay matriz de rotación previa se  
                                 %obtiene la Matriz identidad
       end
   else
       %Para las articulaciones prismáticas
       try
           Jv_a(:,k)=RO(:,3,k-1);
       catch
           Jv_a(:,k) = [0,0,1];
       end
           Jw_a(:,k) = [0,0,0];
   end
end

Jv_a = simplify (Jv_a);
Jw_a = simplify (Jw_a);
disp('Jacobiano lineal obtenido de forma analítica');
pretty(Jv_a);
disp('Jacobiano angular obtenido de forma analítica');
pretty(Jw_a)

disp('Velocidad lineal obtenida mediante el Jacobiano lineal')
V = simplify(Jv_a*Qp'); 
pretty(V); 
disp('Velocidad angular obtenida mediante el Jacobiano angular')
W=simplify(Jw_a*Qp'); 
pretty(W); 