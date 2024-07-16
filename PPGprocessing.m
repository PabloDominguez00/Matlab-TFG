
clear 
close all force
clc

%addpath(genpath('C:\Users\Pablo\Desktop\Calculos\biomedical-signal-processing'));
%Encintramos los saltos: signals.TimeTrial(diff(signals.TimeTrial)>0.5)
%(Estan redondeados a precisión de segundos)
%Tarea: interpolar la frecuencia de muestreo a 256, desde la frecuencia
%original (calculada con los huecos dejados con la formula anterior)
%% carga y representacion de los datos
filename = 'Paciente1';
crop = 10; %Cuantos minutos ignoraremos del inicio y el final
signals = parquetread(['C:\Users\Pablo\Desktop\Calculos\' filename 'LH.gzip']);
signals = signals(3:end,:);
%referencia = struct2table(load('C:\Users\Pablo\Desktop\Calculos\Referencia.mat'));
referencia = parquetread(['C:\Users\Pablo\Desktop\Calculos\' filename '_Pulso_SpO2.parquet']);

%Limpiamos puntos incoherentes (80<SPO2<100)

for i=1:length(referencia.Saturacion)
    if (referencia.Saturacion(i)> 100)
        referencia.Saturacion(i)=nan;
    elseif (referencia.Saturacion(i) < 80)
        referencia.Saturacion(i)=nan;
    end
end

%Eliminamos los crop primeros minutos y los crop últimos minutos
referencia=referencia.Saturacion(crop*60:end-crop*60);

clear i

%% remuestreo a 256Hz de los diferentes apartados
fo=256;
%Encontramos los saltos: 
resampled = {};
huecos = find(diff(signals.TimeTrial)>0.5);

if ~isempty(huecos)
    fs=length(huecos);
    for i=1:length(huecos) 
        if i==1  
            if huecos(i)>256*120
                fs(i)=huecos(i)/(signals.TimeTrial(huecos(i))-signals.TimeTrial(1));
            else
                fs(i)=nan;
            end
        elseif (huecos(i)-huecos(i-1)+1)>256*120 %Si el ensayo ha durado mas de 2 minutos
            fs(i)=(huecos(i)-huecos(i-1)+1)/(signals.TimeTrial(huecos(i))-signals.TimeTrial(huecos(i-1)+1));
        else
            fs(i)=NaN;
        end
    end
    for i=1:length(fs)
        if ~isnan(fs(i))
            if i==1
                g = signals.Green_Count(huecos(1):huecos(i));
                r = signals.Red_Count(huecos(1):huecos(i));
                ir = signals.IR_Count(huecos(1):huecos(i));
            else
                g = signals.Green_Count(huecos(i-1)+1:huecos(i));
                r = signals.Red_Count(huecos(i-1)+1:huecos(i));
                ir = signals.IR_Count(huecos(i-1)+1:huecos(i));
                
            end
            fi = fs(i);
            ffin = fo; %fo=frecuencia objetivo
            multiplicador = 100;
            GC = resample(double(g),ffin*multiplicador,fix(round(fi*multiplicador,2)));
            RC = resample(double(r),ffin*multiplicador,fix(round(fi*multiplicador,2)));
            IC = resample(double(ir),ffin*multiplicador,fix(round(fi*multiplicador,2)));
            t = signals.TimeTrial(1):1/fo:length(RC)/fo-1/fo;
    
            %Eliminamos los crop primeros minutos y los crop últimos minutos
            GC=GC(crop*60*fo:end-crop*60*fo);
            RC=RC(crop*60*fo:end-crop*60*fo);
            IC=IC(crop*60*fo:end-crop*60*fo);
            t=t(crop*60*fo:end-crop*60*fo);

            resampled{end+1}=table(t',GC,RC,IC, 'VariableNames', {'Time', 'GC', 'RC', 'IC'});
        end
    end
    remuestreado =  resampled{1}; %Los intervalos a analizar se encuentran en 1
    %segundos= 1:1:length(referencia.ans); %Referencia sin decimales
else 
    fs=length(signals.TimeTrial)/signals.TimeTrial(end);
    g = signals.Green_Count;
    r = signals.Red_Count;
    ir = signals.IR_Count;

    fi = fs;
    ffin = fo; %fo=frecuencia objetivo
    multiplicador = 100;
    GC = resample(double(g),ffin*multiplicador,fix(round(fi*multiplicador,2)));
    RC = resample(double(r),ffin*multiplicador,fix(round(fi*multiplicador,2)));
    IC = resample(double(ir),ffin*multiplicador,fix(round(fi*multiplicador,2)));
    t = signals.TimeTrial(1):1/fo:length(RC)/fo-1/fo;
    
    %Eliminamos los crop primeros minutos y los crop últimos minutos
    GC=GC(crop*60*fo:end-crop*60*fo);
    RC=RC(crop*60*fo:end-crop*60*fo);
    IC=IC(crop*60*fo:end-crop*60*fo);
    t=t(crop*60*fo:end-crop*60*fo);

    remuestreado=table(t',GC,RC,IC, 'VariableNames', {'Time', 'GC', 'RC', 'IC'});
    %resampled = [resampled ; table(t',GC,RC,IC, 'VariableNames', {'Time', 'GC', 'RC', 'IC'})]
        
    %figure
    %hold on
    %plot(t,GC,'m')
    %plot(t,RC,'c')
    %plot(t,IC, 'y')
    %plot(signals.TimeTrial(huecos(i-1)+1:huecos(i))+17.5,signals.Green_Count(huecos(i-1)+1:huecos(i)),'g')
    %plot(signals.TimeTrial(huecos(i-1)+1:huecos(i))+17.5,signals.Red_Count(huecos(i-1)+1:huecos(i)),'r')
    %plot(signals.TimeTrial(huecos(i-1)+1:huecos(i))+17.5,signals.IR_Count(huecos(i-1)+1:huecos(i)),'b')
end

%Borramos las variable que no vamos a usar mas adelante
clear g r ir fi ffin i GC RC IC t fs huecos multiplicador %crop
%% preprocesado: filtrado y deteccion de artefactos
% signal filtering -- Declaramos la pared para filtrar el ruido (Frecuencia Corte) 
ord = 4;
fc = [0.1 15]./(fo/2);

%Limpia los sonidos que no pertenecen, manteniendo los que quieres altos y
%claros
[z,p,k]         =   cheby2( ord , 20 , fc , 'bandpass' ); 

%Simplifica la señal para que sea mas sencilla de entender
[sos,filterGain] = zp2sos(z,p,k);

%Aplica el filtro a la señal
Setup.plotflag  =	false; 
PPGGreen      =   nanfiltfilt( sos , filterGain , -remuestreado.GC(:) , Setup );
PPGIR         =   nanfiltfilt( sos , filterGain , -remuestreado.IC(:) , Setup );  
PPGRed        =   nanfiltfilt( sos , filterGain , -remuestreado.RC(:) , Setup );  

Setup.plotflag  =	false;
isArtifactG      =   energyArtifacts ( PPGGreen , fo , Setup );
isArtifactIR      =   energyArtifacts ( PPGIR , fo , Setup );
isArtifactRed      =   energyArtifacts ( PPGRed , fo , Setup );

PPGIR(isArtifactIR)=nan; 
PPGGreen(isArtifactG)=nan;
PPGRed(isArtifactRed)=nan;

%Obtenemos el instante de tiempo de 
% los maximos y los minimos en la señal filtrada
Setup.plotflag = false;

%Instantes de tiempo de: D->Latido, A->Max, B->Min 
%(partiendo desde el instante de tiempo de medicion inicial)
%Es decir, que son isntantes desde que inicia la tabla válida
[ gD , gA , gB ]    =   pulseDelineation ( PPGGreen , fo , Setup );
[ iD , iA , iB ]    =   pulseDelineation ( PPGIR , fo , Setup );
[ nD , nA , nB ]    =   pulseDelineation ( PPGRed , fo , Setup ); 


%creamos un eje temporal que se ajuste a nuestra señal
t = 0:1/fo:(length(PPGRed)-1)/fo;

% a_iA=iA; a_iB=iB; %para parte de apnea
% a_nA=nA; a_nB=nB;
% 
% n_iA=iA; n_iB=iB; %para parte noche normal
% n_nA=nA; n_nB=nB;  

%Con tiempos arreglados:
irA=iA+remuestreado.Time(1); irB=iB+remuestreado.Time(1);
rA=nA+remuestreado.Time(1); rB=nB+remuestreado.Time(1);
grA=gA+remuestreado.Time(1); grB=gB+remuestreado.Time(1);

clear k p z sos filterGain

%% Prueba referencia verde para obtencion maximos rojo
MaxRojo=zeros(size(rA)); MaxInfra=zeros(size(irA));
MinRojo=zeros(size(rB)); MinInfra=zeros(size(irB));
i=1;
umbral=0.1;
%referencia irA
while  i <= length(rA) && i <= length(irA)
    inf=irA(i)-umbral;
    sup=irA(i)+umbral;
    [idx,]=find(rA(1:end)>inf & rA(1:end)<sup,1);
    if ~isempty(idx)
        MaxRojo(i)=rA(idx);
        MaxInfra(i)=irA(i);
    else
        MaxRojo(i)=nan;
        MaxInfra(i)=nan;
    end
    i=i+1;
end

i=1;
%Referencia irB
while i <= length(rB) && i <= length(irB)
    inf=irB(i)-umbral;
    sup=irB(i)+umbral;
    [idx,]=find(rB(ant:end)>inf & rB(ant:end)<sup,1);
    if ~isempty(idx)
        MinRojo(i)=rB(idx);
        MinInfra(i)=irB(i);
    else
        MinRojo(i)=nan;
        MinInfra(i)=nan;
    end
    i=i+1;
end

% Quitar 0 al final para dejar ambos indices del mismo tamaño
last_nonzero_idx = find(MaxRojo ~= 0, 1, 'last');
if ~isempty(last_nonzero_idx)
    MaxRojo = MaxRojo(1:last_nonzero_idx);
end
last_nonzero_idx = find(MaxInfra ~= 0, 1, 'last');
if ~isempty(last_nonzero_idx)
    MaxInfra = MaxInfra(1:last_nonzero_idx);
end
last_nonzero_idx = find(MinRojo ~= 0, 1, 'last');
if ~isempty(last_nonzero_idx)
    MinRojo = MinRojo(1:last_nonzero_idx);
end
last_nonzero_idx = find(MinInfra ~= 0, 1, 'last');
if ~isempty(last_nonzero_idx)
    MinInfra = MinInfra(1:last_nonzero_idx);
end

irA=MaxInfra;rA=MaxRojo;irB=MinInfra;rB=MinRojo;
%% Referencia de la verde para crear un vector coherente:
% irMax=nan(length(irA),1); rMax = nan(length(rA),1);
% irMin=nan(length(irB),1); rMin = nan(length(rB),1);
% umbral = 0.15;
% for i=1:length(grA)
%     k=0;
%     while rA(k+i)<=grA(i)+umbral && k+i<length(rA)
%         k=k+1;
%     end
%     if abs(grA(i)-rA(i+k-1))<=umbral
%             rMax(i+k-1)=rA(i+k-1);
%         k=0;
%         while irA(k+i)<=grA(i)+umbral && k+i<length(irA)
%             k=k+1;
%         end
%         if abs(grA(i)-irA(i+k-1))<=umbral
%                 irMax(i+k-1)=irA(i+k-1);
%         end    
%     end
% 
% end
% 
% for i=1:length(grB)
%     k=0;
%     while rB(k+i)<=grB(i)+umbral && k+i<length(rB)
%         k=k+1;
%     end
%     if abs(grB(i)-rB(i+k-1))<=umbral
%             rMin(i+k-1)=rB(i+k-1);
%     end
%     k=0;
%     while irB(k+i)<=grB(i)+umbral && k+i<length(irB)
%         k=k+1;
%     end
%     if abs(grB(i)-irB(i+k-1))<=umbral
%             irMin(i+k-1)=irB(i+k-1);
%     end
% end
% 
% sum(~isnan(irMin))
% sum(~isnan(irMax))
% sum(~isnan(rMax))
% sum(~isnan(rMin))
%% Crear vectores definitivos lo mas parecidos posible
% i=1; j=1;
% while i <= length(irMax) && j <= length(rMax)
%     if isnan(irMax(i)) && isnan(rMax(j))
%         irMax(i)=[]; rMax(j)=[];
%     end
%     i=i+1; j=j+1;
% end
% 
% i=1; j=0;
% while i <= length(irMax)
%     if ~isnan(irMax(i))
%         irMax(i)
% 
%     end
%     % j=0;
%     % while irMax(i)<=irMax(i+j)+0.1 && i <= length(irMax)-1
%     % end
%     i=i+1;
% end
%% Prueba aunar latidos
%Manejar casos: si el anterior es nan entonces feo,  
% pi=irA(1:100); pr=rA(1:100);
% auxI=[]; auxR=[];
% i = 1; j=1;

%En algun momento se descompensan los indices o borramos algo que no
%deberiamos

% while i <= length(pi) && j <= length(pr)
%     %Son los primeros y además hay una distancia coherente entre ellos
%     if i<=2 && j<=2 && abs(pi(i)-pr(j))<=0.075 
%         i=i+1; j=j+1;
%     %Son los primeros pero no hay distancia coherente entre ellos
%     elseif i<=2 && j<=2 && abs(pi(i)-pr(j))>0.075
%         %La distancia con el siguiente es coherente
%         if pi(i+1)-pi(i)>=0.4
%             %Resultado no valido, anulamos ambos
%             pi(i)=nan; pr(j)=nan;
%             i=i+1;
%             j=j+1;
%         %La distancia con el siguiente no es coherente
%         elseif pi(i+1)-pi(i)<0.4
%             %Borramos hasta lograr coherencia
%             while abs(pi(i)-pr(j))>0.075
%                 pi(i)=[];
%             end
%             i=i+1; j=j+1;
%         %Fallo en rojo        
%         elseif pr(i+1)-pr(i)<0.4
%             while abs(pr(j)-pi(i))<0.075
%                 pr(j)=[];
%             end
%             i=i+1; j=j+1;
%         end
%     %Hay una distancia coherente entre ellos
%     elseif abs(pi(i)-pr(j))<=0.075
%         i=i+1; j=j+1;
%     %No hay distancia coherente entre ellos
%     elseif abs(pi(i)-pr(j))>0.075
%         %Hay una distancia correcta entre latidos, pero estan demasiado
%         %separados entre ellos, por tanto invalidamos ambos
%         if abs((pi(i)-pi(i-1))-(pr(i)-pr(i-1)))>=0.09 
%             pi(i)=nan;pr(j)=nan;
%             i=i+1; j=j+1;
%         %Evitamos problemas con nan
%         elseif isnan(pi(i)-pi(i-1))
%             pi(i)=nan;
%             i=i+1;
%             %Revisamos tambien en la otra tabla
%             if isnan(pr(j)-pr(j-1))
%                 pr(j)=nan;
%             end
%             j=j+1;
%         elseif isnan(pr(j)-pr(j-1))
%             pr(j)=nan;
%             j=j+1;
%             i=i+1;
%         %La distancia con el siguiente no es coherente
%         elseif pi(i+1)-pi(i)<0.4
%             %Borramos hasta lograr coherencia
%             while abs(pi(i)-pr(j))<0.075
%                 pi(i)=[];
%             end
%             i=i+1; j=j+1;
%         %Fallo en rojo        
%         elseif pr(i+1)-pr(i)<0.4
%             while abs(pr(j)-pi(i))<0.075
%                 pr(j)=[];
%             end
%             i=i+1; j=j+1;
%         %Si una distancia es coherente y la otra se lleva demasiado,
%         %entonces hemos perdido un dato
%         %Elegimos lado
%         elseif pi(i)>pr(j)
%             aux=1;
%             %Hemos encontrado un hueco, lo llenamos con nan
%             while (pi(i)-pr(j))>-0.075
%                 aux=aux+1;
%             end
%             if abs(pi(i)-pr(j))>0.075 %dato bueno
%                     %añadir aux Nans
%                     for p=1:aux
%                         pj=[pj(1:j-1);nan;pj(j:end)];
%                     end
%             else
%                 %En caso de fallo de lo anterior
%             end
%         elseif pr(j)<pi(i)
%             aux=1;
%             %Hemos encontrado un hueco, lo llenamos con nan
%             while (pj(j)-pi(i))>-0.075
%                     aux=aux+1;
%             end
%             if abs(pi(i)-pr(j))>0.075 %dato bueno
%                 %añadir aux Nans
%                 for p=1:aux
%                     pi=[pi(1:i-1);nan;pi(i:end)];
%                 end
%             else
%                 %en caso de fallo de lo anterior
%             end
%         end
%     end
% end

%Propuesta Rodrigo
% while k < max(length(pi),length(pr))-5
%     rojos=(abs(pi(k:k+4)-pr(k))<0.075)'
%     infra=(abs(pr(k:k+4)-pi(k))<0.075)'
% 
%     if find(rojos,1)==1 && find(infra,1)==1 %correcto
%         k=k+1
%     elseif k==1
% 
%     elseif k<3
% 
% 
%     elseif sum(infra)>0 && sum(rojos)==0
%         if pr(k+1)-pr(k)<0.6*(pr(k-1)-pr(k-2))%El espacio es pequeño y por tanto los numeros estan mal
%             for i=find(infra,1)-1
%                 pi(i)=[]
%             end
%         elseif pi(k)-pi(k-1)>1.2*(pi(k-1)-pi(k-2)) %Hay un hueco en ir
%             pi=[pi(1:find(infra,1)-1);nan;pi(find(infra,1):end)]
%             % pi=[pi(1:k);nan;pi(k:end)]
%         end
%     elseif  sum(rojos)>0 && sum(infra)==0
%         if pi(k+1)-pi(k)<0.6*(pi(k-1)-pi(k-2))%El espacio es pequeño y por tanto los numeros estan mal
%             for i=find(rojos,1)-1
%                 pr(i)=[]
%             end
%         elseif pr(k)- pr(k-1)>1.2*(pr(k-1)-pr(k-2)) %Hay un hueco en rojos
%             pr=[pr(1:find(rojos,1)-1);nan;pr(find(rojos,1):end)]
%             % pr=[pr(1:k);nan;pr(k:end)]
% 
%         end
% 
%     end
% end

%aux = 1;

%Intento 2 (iba decente pero falla en caso de que los dos datos sean malos)
% while  k < max(length(pi),length(pr))-1
%     if abs(pi(k)-pr(k)) > 0.075 
%          if pi(k+1)-pi(k) < 0.3 %error en IR ¿Umbral correcto?
%             timeDiff=pr(k)-pr(k-1);
%             while pi(k) < pi(k-1)+timeDiff -0.075 %Sabemos que el anterior es correcto
%                 pi(k)=[];
%             end
%          elseif pr(k+1)-pr(k) < 0.3 %error en R ¿Umbral correcto?
%             timeDiff=pi(k)-pi(k-1);
%             while pr(k) < pr(k-1)+timeDiff-0.075
%                 pr(k)=[];
%             end
%          elseif abs(pi(k)-pr(k+1))<0.075 %Falta un dato en IR
%             pi=[pi(1:k-1),nan,pi(k:end)];
%             k=k+1;
%          elseif abs(pr(k)-pi(k+1))<0.075 %Falta un dato en R
%             pr=[pr(1:k-1),nan,pr(k:end)];
%             k=k+1;
%          elseif 
%             pi(k)
%             pr(k)
%             pi(k)=nan;
%             pi(k)=nan;
%             k=k+1;
%          else
%             disp("achucha")
%             k=k+1;
%          end
%     else
%         k=k+1;
%     end
% end
% 
% sum(abs(pr-pi)<0.075,"omitnan")

%Intento 1
% while k < max(length(pi),length(pr))-1
%     if abs(pi(k)-pr(k)) > 0.05 
%         if abs(pi(k)-pi(k+1)) < 0.3 %error en IR
%             tiempOk=pr(k+1)-pr(k);
%             %Mientras la diferencia entre las distancias de tiempos entre
%             %la columna correcta y la incorrecta no esté entre +-0.05
%             %seguimos contando
%             while abs((pi(k)+tiempOk)-(pi(k+aux)-pi(k)))>0.05 %& (abs(pr(k+1)-pi(k+aux))>) 
%                 aux=aux+1;
%             end
%             for i = aux
%                 pi(k+aux)=[];
%             end
%             aux=1;
%         elseif abs(pr(k)- pr(k+1)) < 0.3 %error en R
% 
%         end
%     else 
%         k=k+1;
%     end
% end

%% Creamos nuevo filtro para mantener la respiración
%IR
fc = 15/(fo/2);
[z,p,k]         =   cheby2( ord , 20 , fc , 'low' ); 
[sos,filterGain] = zp2sos(z,p,k);
FiltroPBajoIR = nanfiltfilt( sos , filterGain , -remuestreado.IC(:) , Setup );

figure;
hold on;
%Dibujamos la señal ahora con la respiracion
plot(t,  FiltroPBajoIR , 'k','LineWidth',1);
%Marcamos los puntos maximos y minimos en funcion de los tiempos obtenidos
%por la totalmente filtrada en la parcialmente filtrada (solo con la
%respiracion), lo que nos da una precision de m +-1 por culpa de la
%aproximación usada en round
plot(iD(~isnan(iD)), FiltroPBajoIR(1+round(iD(~isnan(iD))*fo)), 'ro','LineWidth',1);
plot(irA(~isnan(irA))-crop*60, FiltroPBajoIR(1+round(irA(~isnan(irA))*fo)-crop*60*fo), 'b*','LineWidth',1);
plot(irB(~isnan(irB))-crop*60, FiltroPBajoIR(1+round(irB(~isnan(irB))*fo)-crop*60*fo), 'b*','LineWidth',1);
title('PPGIR');

%Roja
fc = 15/(fo/2);
[z,p,k]         =   cheby2( ord , 20 , fc , 'low' ); 
[sos,filterGain] = zp2sos(z,p,k);
FiltroPBajoRoja = nanfiltfilt( sos , filterGain , -remuestreado.RC(:) , Setup ); %Este nuevo filtro solo quta el ruido, deja l respiración

figure;
hold on;
plot(t,  FiltroPBajoRoja , 'k','LineWidth',1);
plot(nD(~isnan(nD)), FiltroPBajoRoja(1+round(nD(~isnan(nD))*fo)), 'ro','LineWidth',1);
plot(rA(~isnan(rA))-crop*60, FiltroPBajoRoja(1+round(rA(~isnan(rA))*fo)-crop*60*fo), 'b*','LineWidth',1);
plot(rB(~isnan(rB))-crop*60, FiltroPBajoRoja(1+round(rB(~isnan(rB))*fo)-crop*60*fo), 'b*','LineWidth',1);
title('PPGRED');

clear k p z sos filterGain


%% Toda la noche (creacion de las tablas)
noche=table(remuestreado.Time, FiltroPBajoRoja, FiltroPBajoIR,'VariableNames', {'Time', 'Red', 'IR'});

MaxsIr=table(irA, nan(length(irA),1), 'VariableNames', {'Time', 'Value'});
MinsIr=table(irB, nan(length(irB),1), 'VariableNames', {'Time', 'Value'});
MaxsR=table(rA, nan(length(rA),1), 'VariableNames', {'Time', 'Value'});
MinsR=table(rB, nan(length(rB),1), 'VariableNames', {'Time', 'Value'});

%find(MaxsIr.Time<=noche.Time(1), 1); ??
%clear nA nB iA iB

%% Toda la noche (asignacion de valores a las tablas de tiempos)
% Falta aunar en una sola tabla las muestras resampleadas
%Este apartado cuesta un cojon de ejecutar xdddd
interv = 20; %Buscamos coincidencias en los 3 siguientes segundos
tic
prev = 1;
%Maximos de IR
for i=1:length(MaxsIr.Value)
    if ~isnan(MaxsIr.Time(i))
        ini=prev;
         % índice final de la búsqueda, asegurando no exceder el tamaño del array
        fin= min(ini + interv - 1, length(noche.Time));
        idx = find(noche.Time(ini:fin) == MaxsIr.Time(i), 1);
        prev=idx;
        if isempty(idx) %aproximar al valor correcto más cercano
            %Optimizacion para la busqueda, sumamos los elementos true del
            %vector para hayar el indice que buscamos
            ant=sum(noche.Time<MaxsIr.Time(i));
            sig = ant + 1;
            % asegurar que sig no exceda el tamaño del array
            if sig > length(noche.Time)
                sig = length(noche.Time); 
            end
            if noche.IR(i)>noche.IR(sig)
                MaxsIr.Time(i)=noche.Time(ant);
                MaxsIr.Value(i)=noche.IR(ant);
            else
                MaxsIr.Time(i)=noche.Time(sig);
                MaxsIr.Value(i)=noche.IR(sig);
            end
        else %tenemos valor para el tiempo dado
            MaxsIr.Value(i)=noche.IR(idx);
        end
    end
end
toc

tic
prev = 1;
%Mínimos de IR
for i=1:length(MinsIr.Value)
    if ~isnan(MinsIr.Time(i))
        ini=prev;
        fin= min(ini + interv - 1, length(noche.Time));
        idx = find(noche.Time(ini:fin) == MinsIr.Time(i), 1);
        prev=idx;
        if isempty(idx) 
            ant=sum(noche.Time<MinsIr.Time(i));
            sig=ant+1;
            if sig > length(noche.Time)
                sig = length(noche.Time); 
            end
            if noche.IR(ant)<noche.IR(sig)
                MinsIr.Time(i)=noche.Time(ant);
                MinsIr.Value(i)=noche.IR(ant);
            else
                MinsIr.Time(i)=noche.Time(sig);
                MinsIr.Value(i)=noche.IR(sig);
            end
        else
            MinsIr.Value(i)=noche.IR(idx);
        end
    end
end
toc

tic
prev = 1;
%Maximos de R
for i=1:length(MaxsR.Value)
    if ~isnan(MaxsR.Time(i))
        ini=prev;
        fin= min(ini + interv - 1, length(noche.Time));
        idx = find(noche.Time(ini:fin) == MaxsR.Time(i), 1);
        prev=idx;
        if isempty(idx) 
            ant=sum(noche.Time<MaxsR.Time(i));
            sig = ant + 1;
            if sig > length(noche.Time)
                sig = length(noche.Time); 
            end
            if noche.Red(ant)>noche.Red(sig)%
                MaxsR.Time(i)=noche.Time(ant);
                MaxsR.Value(i)=noche.Red(ant);
            else
                MaxsR.Time(i)=noche.Time(sig);
                MaxsR.Value(i)=noche.Red(sig);
            end
        else 
            MaxsR.Value(i)=noche.Red(idx);
        end
    end
end
toc

tic
prev = 1;
%Mínimos de R
for i=1:length(MinsR.Value)
    if ~isnan(MinsR.Time(i))
        ini=prev;
        fin= min(ini + interv - 1, length(noche.Time));
        idx = find(noche.Time(ini:fin) == MinsR.Time(i), 1);
        prev=idx;
        if isempty(idx) 
            ant=sum(noche.Time<MinsR.Time(i));
            sig=ant+1;
            if sig > length(noche.Time)
                    sig = length(noche.Time); 
            end
            if noche.Red(ant)<noche.Red(sig)
                MinsR.Time(i)=noche.Time(ant);
                MinsR.Value(i)=noche.Red(ant);
            else
                MinsR.Time(i)=noche.Time(sig);
                MinsR.Value(i)=noche.Red(sig);
            end
        else
            MinsR.Value(i)=noche.Red(idx);
        end
    end
end
toc


%% Toda la noche (calculo de R y SPO2)

%Revisar como podriasmos solucionar el tema de los latido totales:
latidos=min(height(MaxsR),height(MaxsIr));
%Declaramos por optimización
red_dc=zeros(latidos-1,1);
ir_dc=zeros(latidos-1,1);
red_ac=zeros(latidos-1,1);
ir_ac=zeros(latidos-1,1);
aR=zeros(latidos-1,1);
tic
for i=1:latidos-1
    %Cogemos el minimo de un latido y el siguiente y sacamos la media del
    %valor de la señal entre esos puntos

    %Prueba 
    % red_dc(i) = mean(MinsR.Value(i):MinsR.Value(i+1));
    % ir_dc(i) = mean(MinsIr.Value(i):MinsIr.Value(i+1));

    red_dc(i) = mean(noche.Red(find(noche.Time==MinsR.Time(i)):find(noche.Time==MinsR.Time(i+1))));
    ir_dc(i) = mean(noche.IR(find(noche.Time==MinsIr.Time(i)):find(noche.Time==MinsIr.Time(i+1))));

    red_ac(i) = ((MaxsR.Value(i)-MinsR.Value(i)+(MaxsR.Value(i+1)-MinsR.Value(i+1))))/2;
    ir_ac(i) = ((MaxsIr.Value(i)-MinsIr.Value(i)+(MaxsIr.Value(i+1)-MinsIr.Value(i+1))))/2;

    a=(red_ac(i)/red_dc(i));
    b=(ir_ac(i)/ir_dc(i));
    c=a/b;
    aR(i)=c;
    %j=j+1;
end
toc
%% Toda la noche (calculo FFT y limpieza de la señal)
% fourierApn=fft(aR);
% aL=length(aR);
% 
% A2=abs(fourierApn/aL);
% A1=A2(1:aL/2+1);
% A1(2:end-1)= 2*A1(2:end-1);
% fApn=fo/aL*(0:(aL/2));

%En vistas de las graficas nos interesa establecer el corte en 15Hz, ya que
%es el punto suficientemente alejado como para no perder precisión, pero
%suficientemente cercano para limpiar bien la señal.
fc = 20/(fo/2); 
[z,p,k]         =   cheby2( ord , 20 , fc , 'low' ); 
[sos,filterGain] = zp2sos(z,p,k);
aRFilt = nanfiltfilt( sos , filterGain , -aR, Setup );

SPO2=zeros(1,length(aR));
SPO2Filt= zeros(1,length(aR));

for i=1:length(SPO2)
    SPO2Filt(i)=(110-aRFilt(i));
    SPO2(i)=(110-aR(i));
end

%% Toda la noche (dibujo saturación con apnea)
%Juntamos tiempo y saturación de una apnea
% inicio_estudio=floor(noche.Time(1));
% fin_estudio=ceil(noche.Time(end));
% referencia_noche=referencia(1:end-1,1);
% segundos = 1:1:length(referencia)-1;
% segundos=segundos';
% ref_noche=table(referencia_noche, segundos, 'VariableNames', {'referencia_noche', 'segundos'});

inicio_estudio=floor(noche.Time(1));
fin_estudio=ceil(noche.Time(end));
referencia_noche=referencia(inicio_estudio:fin_estudio-1,1);
segundos = inicio_estudio:1:fin_estudio-1;
segundos=segundos';
ref_noche=table(referencia_noche, segundos, 'VariableNames', {'referencia_noche', 'segundos'});
% figure
% subplot(2,1,1)
% hold on
% plot(MaxsR.Time(1:end-5), red_dc,'DisplayName', 'Red DC'); %* revisar las longitudes
% plot(MaxsR.Time(1:end-5), red_ac,'DisplayName', 'Red AC'); %*
% plot(MaxsIr.Time(1:end-1), ir_dc,'DisplayName', 'Infra DC');
% plot(MaxsIr.Time(1:end-1), ir_ac,'DisplayName', 'Infra AC');
% title('Grafico apnea con respiración media 2')
% legend
 
% prueba1=MaxsR.Time+31;

figure
% subplot(2,1,2)
hold on
plot(ref_noche.segundos, ref_noche.referencia_noche, 'DisplayName', 'Sat Referencia');
yyaxis right
%plot(MaxsIr.Time(1:end-1),SPO2,'DisplayName', 'SPO2');
%plot(prueba1(1:end-1),SPO2,'DisplayName', 'SPO2');
%plot(MaxsR.Time(1:end-1),aR+91,'DisplayName', 'R+91');
% ylim([111 112])
plot(MaxsIr.Time(1:end-1)-600+40,-aRFilt,'DisplayName', 'R apnea low-15');%El -560 esta puesto a ojo
% plot(MaxsR.Time(1:end-1),-aRFilt,'DisplayName', 'R apnea low-15');
legend


%% preprocesado: filtrado y deteccion de artefactos

% signal filtering -- Declaramos la pared para filtrar el ruido (Frecuencia Corte) 
ord = 4;
fc = [0.1 15]./(fo/2);

% [bb,aa]           =   butter ( ord , fc , 'bandpass' );
% Setup.plotflag    =	true;
% PPGsignal         =   nanfiltfilt ( bb , aa , signals.Green_Count(:) );

%Limpia los sonidos que no pertenecen, manteniendo los que quieres altos y
%claros
[z,p,k]         =   cheby2( ord , 20 , fc , 'bandpass' ); 

%Simplifica la señal para que sea mas sencilla de entender
[sos,filterGain] = zp2sos(z,p,k);

Setup.plotflag  =	false;
% PPGsignal       =   nanfiltfilt( sos , g , -signals.Green_Count(:) , Setup );  

%Filtra la señal, y le da la vuelta y la vuelve a filtrar
PPGIR       =   nanfiltfilt( sos , filterGain , -IC(:) , Setup );  
PPGRed      =   nanfiltfilt( sos , filterGain , -RC(:) , Setup );  

disp('CHECK IF PPG SIGNAL IS INVERTED!!!');


% Eliminar Artefactos --> detecta cambios bruscos en la señal, para
% borrarlos
Setup.plotflag  =	false;
isArtifact      =   energyArtifacts ( PPGIR , fo , Setup );
isArtifact      =   energyArtifacts ( PPGRed , fs , Setup );


%% deteccion y delineacion de la señal de PPG

Setup.isArtifact    =   false;
Setup.plotflag      =   true;
plotflag = true;
[ snD , snA , ssnB ]    =   pulseDelineation ( PPGIR , fo , Setup );  % PREGUNTAR POR ESE REDONDE A FS SI AFECTA
% [ snD , snA , snB ]    =   pulseDelineation ( PPGsignal , fs , Setup );  % PREGUNTAR POR ESE REDONDE A FS SI AFECTA


%% EDF de apnealink treatment
% data = edfread(['C:\BSICoS\EnsayosApnealink\' filename '.edf']);
% info = edfinfo(['C:\BSICoS\EnsayosApnealink\' filename '.edf']);
% 
% fs_AL = info.NumSamples/seconds(info.DataRecordDuration);
% 
% 
% for fd = data.Properties.VariableNames
%     signals_apnealink.(genvarname(fd{:}))=[];
% end
% 
% 
% for fd = data.Properties.VariableNames
%     aux=  [data.(fd{:}){:}];
%     signals_apnealink.(fd{:}) =aux(:);
% 
% end
% 
% t60 =  seconds(0:1:(info.NumDataRecords*60-1));
% t600 = seconds(0:0.1:(info.NumDataRecords*60-0.1));
% t6000 = seconds(0:0.01:(info.NumDataRecords*60-0.01));


%% Plot Saturation
% figure
% hold on 
% plot(t60,signals_apnealink.SignalLabel5_Saturacion);
% legend(info.SignalLabels(5))
% xlabel('Tiempo')
% xtickformat('hh:mm:ss')
% legend
% hold off

%% Plot Apnea 

% figure
% hold on 
% plot(t600,signals_apnealink.SignalLabel8_ApneaObstructiv);
% plot(t600,signals_apnealink.SignalLabel9_ApneaCentral);
% plot(t600,signals_apnealink.SignalLabel10_ApneaMixta);
% legend(info.SignalLabels(8),info.SignalLabels(9),info.SignalLabels(10))
% xlabel('Tiempo')
% xtickformat('hh:mm:ss')
% legend
% hold off


%% calculo del HRV (en realidad, como trabajamos con PPPG, se llama PRV, pulse rate variability)

% % , utilizando como punto fiducial el nD.
% se podria repetir este analisis usando, por ejemplo, el punto nB
fr = 4; % frecuencia de resampling para remuestreo a fs uniforme.
Setup.plotflag = true;
[ m , NN , tNN, mt , mtt ,a,b,c,s,tHR] = computeHRVsignals ( nD , fr , Setup );


% Find the nans
nanLocations = isnan(mtt);
% Pick some other value to set the nans to.
alternateValue =mean(mtt,"omitnan");
% Do the replacemenet.
mtt(nanLocations) = alternateValue;
[s,lags] = xcorr(signals_apnealink.SignalLabel4_Pulso,downsample(mtt*60,4));
offset = lags(find(s==max(s)));

%% Opcion 2 ploteo HR

figure
hold on

%plot(t60,signals_apnealink.SignalLabel4_Pulso, 'DisplayName','Reference HR', color = '#0072BD');
% plot(tHR+offset, mt*60, 'DisplayName','instantaneous HRFS [bpm]','LineWidth',1,'Color', "#FF5599");
%plot(tHR+offset, mtt*60, 'DisplayName','mean HR [bpm]','LineWidth',1,'Color', "#000000");

% plot(tHR-10+floor(tHR./60)*60*(1-(fs+0.1)/256), mt*60, 'DisplayName','instantaneous HR [bpm]','LineWidth',1,'Color', "#D95319");
% plot(tHR-10+floor(tHR./60)*60*(1-(fs+0.1)/256), mtt*60, 'DisplayName','mean HR [bpm]','LineWidth',1,'Color', "#000000");

xlabel('Tiempo')
xtickformat('hh:mm:ss.SS'); %set(gca,'xtickformat','hh:mm:ss')

legend
hold off

%% Ploteo esfuerzo y flujo
figure
hold on
for recnum =1:1:info.NumDataRecords
    for signum = [1,3]
        t = (0:info.NumSamples(signum)-1)/fs_AL(signum);
        y = data.(signum){recnum};
        if signum == 1
            plot(seconds(t+recnum*60),y, color = "#0072BD")
        else
            plot(seconds(t+recnum*60),y, color = "#D95319")
        end

    end
end
legend(info.SignalLabels(1),info.SignalLabels(3))
xlabel('Tiempo')
xtickformat('hh:mm:ss')
legend
hold off

figure
%plot(cumsum(signals_apnealink.SignalLabel1_Flujo))

%% el analisis frecuencial conviene hacerlo a cachos de 10 mins del HRV
% 
% ini = 4*60*60*fr; % cogemos un cacho de 5mins, que empieza en la hora 4 del registro
% fin = ini + 10*60*fr;
% 
% Setup.plotflag = true;
% [ pVLF, pLF, pHF, pLFn, LF_HF ] = frequencyIndices ( m(ini:fin) , fr , Setup );

% %% analisis temporal del HRV
% [ mHR , SDNN , SDSD , RMSSD , pNN50  ] = temporalIndices ( NN )
% 
% 
% %% Sobreponer HRV durante una hora
% 
% for h=1:1:7
%     a = figure;
%     hold on
%     lVLF                                =   0;
%     uVLF                                =   0.04;
%     lLF                                 =   0.04;
%     uLF                                 =   0.15;
%     lHF                             =	0.15;
%     uHF                             =	0.4;
% 
% 
%     for ini= 60*60*fr*h : 10*60*fr : 60*60*fr*(h+1)
% 
%         fin = ini + 10*60*fr;
% 
%         Setup.plotflag = false;
%         [ pVLF, pLF, pHF, pLFn, LF_HF, TP,P, f ] = frequencyIndices ( m(ini:fin) , fr , Setup );
%         % pVLF, pLF, pHF, pLFn, LF_HF, TP, P, f, P_VLF
%         plot( f , P , 'DisplayName', 'Signal 1' );
% 
% 
%     end
%     legend('0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60')
    % legend('0-10','10-20','20-30','30-40','40-50','50-60')
%     xlim([0 uHF+0.1]);   
%     xlabel('Freq [Hz]');ylabel('PSD $[ms^{-2}]$','Interpreter','latex');
%     title (['HRV Freq Indices Hour ' num2str(h,'%d')]); 
%     % title (['HRV Freq Indices Hour ']); 
%     xline( lLF , 'k:');     xline( uLF , 'k:');
%     xline( lHF , 'k--');	xline( uHF , 'k--');    
% 
% 
%     hold off
%     saveas(a, ['H' num2str(h,'%d') '.png'],'png')
% end
% 
% toc

