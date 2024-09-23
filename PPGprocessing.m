
clear 
close all force
clc

addpath(genpath('C:\Users\Pablo\Desktop\Calculos\biomedical-signal-processing'));
%cargaBio;

%% carga y representacion de los datos
filename = 'Paciente0';
crop = 10; %Cuantos minutos ignoraremos del inicio y el final
signals = parquetread(['C:\Users\Pablo\Desktop\Calculos\' filename 'LH.gzip']);signals = signals(3:end,:);
referencia = parquetread(['C:\Users\Pablo\Desktop\Calculos\' filename '_Pulso_SpO2.parquet']);

%Limpiamos puntos incoherentes (80<SPO2<100)

for i=1:length(referencia.Saturacion)
    if (referencia.Saturacion(i)> 100)
        referencia.Saturacion(i)=nan;
    elseif (referencia.Saturacion(i) < 80)
        referencia.Saturacion(i)=nan;
    end
end

%Eliminamos los CROP primeros minutos y los CROP últimos minutos
referencia=referencia.Saturacion(crop*60:end-crop*60);
%sumamos los lugares donde sea menor que el tiempo buscado para hallar el
%índice rápido
signals=signals(sum(signals.TimeTrial<crop*60)+1:sum(signals.TimeTrial<signals.TimeTrial(end)-crop*60)-1,:);

clear i

%% Resampled
fo=256;
huecos = find(diff(signals.TimeTrial)>0.5);
i=1;

if ~isempty(huecos) 
    if huecos(i)-signals.TimeTrial(1)>120*256
    
        g = signals.Green_Count(1:huecos(i));
        r = signals.Red_Count(1:huecos(i));
        ir = signals.IR_Count(1:huecos(i));            
        
        t = signals.TimeTrial(1):1/fo:signals.TimeTrial(huecos(i));
        GC = interp1(signals.TimeTrial(1:huecos(i)), double(g), t,'pchip' );
        RC = interp1(signals.TimeTrial(1:huecos(i)), double(r), t ,'pchip' );
        IC = interp1(signals.TimeTrial(1:huecos(i)), double(ir), t ,'pchip' );

        clear g r ir fs
    
        resampled = table(t',GC',RC',IC', 'VariableNames', {'Time', 'GC', 'RC', 'IC'});
    else
        %Time=nan(huecos(1),1);
        Time=signals.TimeTrial(1):1/fo:signals.TimeTrial(huecos(i));
        GC=nan(huecos(1),1);
        RC=nan(huecos(1),1);
        IC=nan(huecos(1),1);
        resampled = table(Time,GC,RC,IC);
    end
    %frecuencia de muestreo tramos intermedios
    while i <= length(huecos)-1
        %Llenamos de nan todos los instantes de tiempo donde no hay datos
        % n=(signals.TimeTrial(huecos(i)+1)-signals.TimeTrial(huecos(i)));
        % n=floor(n*256);
        % Time=nan(n,1);
        Time=signals.TimeTrial(huecos(i)):1/fo:signals.TimeTrial(huecos(i)+1);
        Time=Time';
        [n,~]=size(Time);
        GC=nan(n,1);
        RC=nan(n,1);
        IC=nan(n,1);
        resampled = [resampled;table(Time,GC,RC,IC)];

        %Cumple condición de durar mas de 2 minutos
        if huecos(i+1)-huecos(i)+1>120*256
            g = signals.Green_Count(huecos(i)+1:huecos(i+1));
            r = signals.Red_Count(huecos(i)+1:huecos(i+1));
            ir = signals.IR_Count(huecos(i)+1:huecos(i+1));            

            t = signals.TimeTrial(huecos(i)+1):1/fo:signals.TimeTrial(huecos(i+1));
            GC = interp1(signals.TimeTrial(huecos(i)+1:huecos(i+1)), double(g), t,'pchip' );
            RC = interp1(signals.TimeTrial(huecos(i)+1:huecos(i+1)), double(r), t ,'pchip' );
            IC = interp1(signals.TimeTrial(huecos(i)+1:huecos(i+1)), double(ir), t ,'pchip' );

            clear g r ir
        
            resampled = [resampled; table(t',GC',RC',IC', 'VariableNames', {'Time', 'GC', 'RC', 'IC'})];

        %No la cumple, añadir nan a la tabla
        else
            % Time=nan(length(signals.TimeTrial(huecos(i)+1:huecos(i+1))),1);
            Time=signals.TimeTrial(huecos(i)+1):1/fo:signals.TimeTrial(huecos(i+1));
            Time=Time';
            [n,~]=size(Time);
            GC=nan(n,1);
            RC=nan(n,1);
            IC=nan(n,1);
            resampled = [resampled;table(Time,GC,RC,IC)];
        end

        i=i+1; 
    end
    %Llenamos de nan todos los instantes de tiempo donde no hay datos
    % n=(signals.TimeTrial(huecos(i)+1)-signals.TimeTrial(huecos(i)));
    % n=floor(n*256);
    % Time=nan(n,1);
    Time=signals.TimeTrial(huecos(i)):1/fo:signals.TimeTrial(huecos(i)+1);
    Time=Time';
    [n,~]=size(Time);
    GC=nan(n,1);
    RC=nan(n,1);
    IC=nan(n,1);
    resampled = [resampled;table(Time,GC,RC,IC)];

    if (height(signals)-huecos(i))-huecos(i)>120*256
        g = signals.Green_Count(huecos(i)+1:end);
        r = signals.Red_Count(huecos(i)+1:end);
        ir = signals.IR_Count(huecos(i)+1:end);            
        
        t = signals.TimeTrial(huecos(i)+1):1/fo:signals.TimeTrial(end); %Incluye el crop
        GC = interp1(signals.TimeTrial(huecos(i)+1:end), double(g), t,'pchip' );
        RC = interp1(signals.TimeTrial(huecos(i)+1:end), double(r), t ,'pchip' );
        IC = interp1(signals.TimeTrial(huecos(i)+1:end), double(ir), t ,'pchip' );

        clear g r ir fs
    
        resampled = [resampled;table(t',GC',RC',IC', 'VariableNames', {'Time', 'GC', 'RC', 'IC'})];
    else
        %Time=nan(length(signals.TimeTrial(huecos(i)+1):end),1);
        Time=signals.TimeTrial(huecos(i)+1):1/fo:signals.TimeTrial(end);
        Time=Time';
        [n,~]=size(Time);
        GC=nan(n,1);
        RC=nan(n,1);
        IC=nan(n,1);
        resampled = [resampled;table(Time,GC,RC,IC)];
    end

else 
    %Un único tramo, todo es continuo, no hay que añadir nada

    %Version antigua, borrar si funciona el codigo sin comentar
    % fs=size(signals)/(signals.TimeTrial(end)-crop*60);
    % fs=fs(1);
    % 
    % g = signals.Green_Count;
    % r = signals.Red_Count;
    % ir = signals.IR_Count;            
    % 
    % multiplicador = 100;
    % GC = resample(double(g),fo*multiplicador,fix(round(fs*multiplicador,2)));
    % RC = resample(double(r),fo*multiplicador,fix(round(fs*multiplicador,2)));
    % IC = resample(double(ir),fo*multiplicador,fix(round(fs*multiplicador,2)));
    % t = (0:length(GC)-1) / fo; %signals.TimeTrial(1):1/fo:length(RC)/fo-1/fo; 
    % 
    % clear g r ir fs
    
    fo=256;

    Time = crop*60:1/fo:signals.TimeTrial(end);
    GC = interp1(signals.TimeTrial, double(signals.Green_Count), Time,'pchip' );
    RC = interp1(signals.TimeTrial, double(signals.Red_Count), Time ,'pchip' );
    IC = interp1(signals.TimeTrial, double(signals.IR_Count), Time ,'pchip' );

    resampled = table(Time',GC',RC',IC', 'VariableNames', {'Time', 'GC', 'RC', 'IC'});
end

remuestreado=resampled;

figure
hold on
plot(signals.TimeTrial,signals.Red_Count,'DisplayName', 'Original');
plot(remuestreado.Time,remuestreado.RC,'DisplayName', 'Remuestreado');
legend

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

%% Prueba referencia Infra para obtencion maximos rojo (PROBAR al reves a ver cual da mejores resultados)
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
    [idx,]=find(rB(1:end)>inf & rB(1:end)<sup,1);
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

%% Creamos nuevo filtro para mantener la respiración
%IR
fc = 15/(fo/2);
[z,p,k]         =   cheby2( ord , 20 , fc , 'low' ); 
[sos,filterGain] = zp2sos(z,p,k);
%Este nuevo filtro solo quta el ruido, deja la respiración
FiltroPBajoIR = nanfiltfilt( sos , filterGain , -remuestreado.IC(:) , Setup );


%Roja
fc = 15/(fo/2);
[z,p,k]         =   cheby2( ord , 20 , fc , 'low' ); 
[sos,filterGain] = zp2sos(z,p,k);
FiltroPBajoRoja = nanfiltfilt( sos , filterGain , -remuestreado.RC(:) , Setup ); 

clear k p z sos filterGain


%% Toda la noche (creacion de las tablas)
noche=table(remuestreado.Time, FiltroPBajoRoja, FiltroPBajoIR,'VariableNames', {'Time', 'Red', 'IR'});

MaxsIr=table(irA, nan(length(irA),1), 'VariableNames', {'Time', 'Value'});
MinsIr=table(irB, nan(length(irB),1), 'VariableNames', {'Time', 'Value'});
MaxsR=table(rA, nan(length(rA),1), 'VariableNames', {'Time', 'Value'});
MinsR=table(rB, nan(length(rB),1), 'VariableNames', {'Time', 'Value'});

%find(MaxsIr.Time<=noche.Time(1), 1); ??
%clear nA nB iA iB

%% Prueba reduccion de tiempos con busqueda knn:

windowSize=100;
i=1;
tic
%Maximos de IR
nearestIdx = knnsearch(noche.Time, MaxsIr.Time(1));
MaxsIr.Time(i) = noche.Time(nearestIdx);
MaxsIr.Value(i) = noche.IR(nearestIdx);
for i = 2:length(MaxsIr.Value)
    if ~isnan(MaxsIr.Time(i)) 
        aux = knnsearch(noche.Time(nearestIdx:nearestIdx+256*2), MaxsIr.Time(i));
        nearestIdx = aux + nearestIdx;
        % Asegurarse de que no salimos de los límites de la señal
        lft = max(1, nearestIdx-windowSize);  % Índice inicial(>=1)
        rgt = min(length(noche.Time), nearestIdx+windowSize);  % Índice final (<= que el largo de la señal)
        % Encontrar el valor minimo en esa ventana porque estan en negativo
        [~, max_index] = max(noche.IR(lft:rgt));
        %Max index devuelve un indice relativo al vector recibido, por
        %tanto, ya que el vector va [i-5,i+5], restamos tamañoVentana y le
        %sumamos el indice que ha encontrado en el subvector 
        nearestIdx = lft + max_index;
        MaxsIr.Time(i) = noche.Time(nearestIdx);
        MaxsIr.Value(i) = noche.IR(nearestIdx);     
    end
end
toc

tic
% Mínimos de IR
nearestIdx = knnsearch(noche.Time, MinsIr.Time(1));
MinsIr.Time(1) = noche.Time(nearestIdx);
MinsIr.Value(1) = noche.IR(nearestIdx);
for i = 2:length(MaxsIr.Value)
    if ~isnan(MinsIr.Time(i)) 
        aux = knnsearch(noche.Time(nearestIdx:nearestIdx+256*2), MinsIr.Time(i));
        nearestIdx = aux + nearestIdx;
        lft = max(1, nearestIdx-windowSize);
        rgt = min(length(noche.Time), nearestIdx+windowSize);
        % Encontrar el valor maximo en esa ventana porque los valores son
        % negativos
        [~, max_index] = min(noche.IR(lft:rgt));
        nearestIdx=lft + max_index;
        MinsIr.Time(i) = noche.Time(nearestIdx);
        MinsIr.Value(i) = noche.IR(nearestIdx);     
    end
end
toc

tic
%Maximos de R
nearestIdx = knnsearch(noche.Time, MaxsR.Time(1));
MaxsR.Time(1) = noche.Time(nearestIdx);
MaxsR.Value(1) = noche.Red(nearestIdx);
for i = 2:length(MaxsR.Value)
    if ~isnan(MaxsR.Time(i)) 
        aux = knnsearch(noche.Time(nearestIdx:nearestIdx+256*2), MaxsR.Time(i));
        nearestIdx = aux + nearestIdx;
        lft = max(1, nearestIdx-windowSize);
        rgt = min(length(noche.Time), nearestIdx+windowSize);
        [~, max_index] = max(noche.Red(lft:rgt));
        nearestIdx=lft + max_index;
        MaxsR.Time(i) = noche.Time(nearestIdx);
        MaxsR.Value(i) = noche.Red(nearestIdx);     
    end
end
toc

tic
% Mínimos de R
nearestIdx = knnsearch(noche.Time, MinsR.Time(1));
MinsR.Time(1) = noche.Time(nearestIdx);
MinsR.Value(1) = noche.Red(nearestIdx);
for i = 2:length(MaxsR.Value)
    if ~isnan(MinsR.Time(i)) 
        aux = knnsearch(noche.Time(nearestIdx:nearestIdx+256*2), MinsR.Time(i));
        nearestIdx = aux + nearestIdx;
        lft = max(1, nearestIdx-windowSize);
        rgt = min(length(noche.Time), nearestIdx+windowSize);
        [~, max_index] = min(noche.Red(lft:rgt));
        nearestIdx=lft + max_index;
        MinsR.Time(i) = noche.Time(nearestIdx);
        MinsR.Value(i) = noche.Red(nearestIdx);     
    end
end
toc


%Version antigua
%(Añadir una busqueda de maximos y minimos en torno al dato (ventana de +-5))
% windowSize=5;
% %PROBAR SI FUNCIONA EL TAMAÑO DE VENTANA
% 
% tic
% %Maximos de IR
% nearestIdx = knnsearch(noche.Time, MaxsIr.Time(1));
% MaxsIr.Time(i) = noche.Time(nearestIdx);
% MaxsIr.Value(i) = noche.IR(nearestIdx);
% for i = 2:length(MaxsIr.Value)
%     if ~isnan(MaxsIr.Time(i)) 
%         aux = knnsearch(noche.Time(nearestIdx:nearestIdx+300), MaxsIr.Time(i));
%         nearestIdx = aux + nearestIdx;
%         % Asegurarse de que no salimos de los límites de la señal
%         lft = max(1, nearestIdx-windowSize);  % Índice inicial(>=1)
%         rgt = min(length(MaxsIr.Value), nearestIdx+windowSize);  % Índice final (<= que el largo de la señal)
% 
%         % Extraer la ventana de la señal
%         window = noche.IR(lft:rgt);
% 
%         % Encontrar el valor máximo en esa ventana
%         [max_value, max_index] = max(window);
% 
%         MaxsIr.Time(i) = noche.Time(nearestIdx);
%         MaxsIr.Value(i) = noche.IR(nearestIdx);     
%     end
% end
% toc
% 
% tic
% % Mínimos de IR
% nearestIdx = knnsearch(noche.Time, MinsIr.Time(1));
% MinsIr.Time(i) = noche.Time(nearestIdx);
% MinsIr.Value(i) = noche.IR(nearestIdx);
% for i = 2:length(MaxsIr.Value)
%     if ~isnan(MinsIr.Time(i)) 
%         aux = knnsearch(noche.Time(nearestIdx:nearestIdx+300), MinsIr.Time(i));
%         nearestIdx = aux + nearestIdx;
%         MinsIr.Time(i) = noche.Time(nearestIdx);
%         MinsIr.Value(i) = noche.IR(nearestIdx);     
%     end
% end
% toc
% 
% tic
% %Maximos de R
% nearestIdx = knnsearch(noche.Time, MaxsR.Time(1));
% MaxsR.Time(i) = noche.Time(nearestIdx);
% MaxsR.Value(i) = noche.Red(nearestIdx);
% for i = 2:length(MaxsR.Value)
%     if ~isnan(MaxsR.Time(i)) 
%         aux = knnsearch(noche.Time(nearestIdx:nearestIdx+300), MaxsR.Time(i));
%         nearestIdx = aux + nearestIdx;
%         MaxsR.Time(i) = noche.Time(nearestIdx);
%         MaxsR.Value(i) = noche.Red(nearestIdx);     
%     end
% end
% toc
% 
% tic
% % Mínimos de R
% nearestIdx = knnsearch(noche.Time, MinsR.Time(1));
% MinsR.Time(i) = noche.Time(nearestIdx);
% MinsR.Value(i) = noche.Red(nearestIdx);
% for i = 2:length(MaxsR.Value)
%     if ~isnan(MinsR.Time(i)) 
%         aux = knnsearch(noche.Time(nearestIdx:nearestIdx+300), MinsR.Time(i));
%         nearestIdx = aux + nearestIdx;
%         MinsR.Time(i) = noche.Time(nearestIdx);
%         MinsR.Value(i) = noche.Red(nearestIdx);     
%     end
% end
% toc

%% Dibujo de máximos y mínimos
figure
hold on
plot(noche.Time,noche.IR, "k", 'DisplayName', 'IR');
plot(MaxsIr.Time, MaxsIr.Value, "r*", 'DisplayName', 'MAX'); 
plot(MinsIr.Time, MinsIr.Value, "bo", 'DisplayName', 'MINS'); 
legend

figure
hold on
plot(noche.Time,noche.Red, "k", 'DisplayName', 'RED');
plot(MaxsR.Time, MaxsR.Value, "r*", 'DisplayName', 'MAX');
plot(MinsR.Time, MinsR.Value, "bo", 'DisplayName', 'MINS');
legend


%% OBSOLETO
% Toda la noche (asignacion de valores a las tablas de tiempos)
% %Este apartado cuesta un cojon de ejecutar xdddd
% interv = 20; %Buscamos coincidencias en los 3 siguientes segundos
% 
% prev = 1;
% 
% %Probar a optimizar con la función knnsearch que busca los k valores mas
% %cercanos en un espacio dado
% tic
% %Maximos de IR
% for i=1:length(MaxsIr.Value)
%     if ~isnan(MaxsIr.Time(i))
%         ini=prev;
%          % índice final de la búsqueda, asegurando no exceder el tamaño del array
%         fin= min(ini + interv - 1, length(noche.Time));
%         idx = find(noche.Time(ini:fin) == MaxsIr.Time(i), 1);
%         prev=idx;
%         if isempty(idx) %aproximar al valor correcto más cercano
%             %Optimizacion para la busqueda, sumamos los elementos true del
%             %vector para hayar el indice que buscamos
%             ant=sum(noche.Time<MaxsIr.Time(i));
%             sig = ant + 1;
%             % asegurar que sig no exceda el tamaño del array
%             if sig > length(noche.Time)
%                 sig = length(noche.Time); 
%             end
%             if noche.IR(i)>noche.IR(sig)
%                 MaxsIr.Time(i)=noche.Time(ant);
%                 MaxsIr.Value(i)=noche.IR(ant);
%             else
%                 MaxsIr.Time(i)=noche.Time(sig);
%                 MaxsIr.Value(i)=noche.IR(sig);
%             end
%         else %tenemos valor para el tiempo dado
%             MaxsIr.Value(i)=noche.IR(idx);
%         end
%     end
% end
% toc
% 
% tic
% prev = 1;
% %Mínimos de IR
% for i=1:length(MinsIr.Value)
%     if ~isnan(MinsIr.Time(i))
%         ini=prev;
%         fin= min(ini + interv - 1, length(noche.Time));
%         idx = find(noche.Time(ini:fin) == MinsIr.Time(i), 1);
%         prev=idx;
%         if isempty(idx) 
%             ant=sum(noche.Time<MinsIr.Time(i));
%             sig=ant+1;
%             if sig > length(noche.Time)
%                 sig = length(noche.Time); 
%             end
%             if noche.IR(ant)<noche.IR(sig)
%                 MinsIr.Time(i)=noche.Time(ant);
%                 MinsIr.Value(i)=noche.IR(ant);
%             else
%                 MinsIr.Time(i)=noche.Time(sig);
%                 MinsIr.Value(i)=noche.IR(sig);
%             end
%         else
%             MinsIr.Value(i)=noche.IR(idx);
%         end
%     end
% end
% toc
% 
% tic
% prev = 1;
% %Maximos de R
% for i=1:length(MaxsR.Value)
%     if ~isnan(MaxsR.Time(i))
%         ini=prev;
%         fin= min(ini + interv - 1, length(noche.Time));
%         idx = find(noche.Time(ini:fin) == MaxsR.Time(i), 1);
%         prev=idx;
%         if isempty(idx) 
%             ant=sum(noche.Time<MaxsR.Time(i));
%             sig = ant + 1;
%             if sig > length(noche.Time)
%                 sig = length(noche.Time); 
%             end
%             if noche.Red(ant)>noche.Red(sig)%
%                 MaxsR.Time(i)=noche.Time(ant);
%                 MaxsR.Value(i)=noche.Red(ant);
%             else
%                 MaxsR.Time(i)=noche.Time(sig);
%                 MaxsR.Value(i)=noche.Red(sig);
%             end
%         else 
%             MaxsR.Value(i)=noche.Red(idx);
%         end
%     end
% end
% toc
% 
% tic
% prev = 1;
% %Mínimos de R
% for i=1:length(MinsR.Value)
%     if ~isnan(MinsR.Time(i))
%         ini=prev;
%         fin= min(ini + interv - 1, length(noche.Time));
%         idx = find(noche.Time(ini:fin) == MinsR.Time(i), 1);
%         prev=idx;
%         if isempty(idx) 
%             ant=sum(noche.Time<MinsR.Time(i));
%             sig=ant+1;
%             if sig > length(noche.Time)
%                     sig = length(noche.Time); 
%             end
%             if noche.Red(ant)<noche.Red(sig)
%                 MinsR.Time(i)=noche.Time(ant);
%                 MinsR.Value(i)=noche.Red(ant);
%             else
%                 MinsR.Time(i)=noche.Time(sig);
%                 MinsR.Value(i)=noche.Red(sig);
%             end
%         else
%             MinsR.Value(i)=noche.Red(idx);
%         end
%     end
% end
% toc


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
% idxStartR = knnsearch(noche.Time,MinsR.Time(1));    
% idxEndR = knnsearch(noche.Time, MinsR.Time(2));
% red_dc(1) = mean(noche.Red(idxStartR:idxEndR));
% 
% idxStartIR = knnsearch(noche.Time, MinsIr.Time(1));
% idxEndIR = knnsearch(noche.Time, MinsIr.Time(2));
% ir_dc(1) = mean(noche.IR(idxStartIR:idxEndIR));

for i=1:latidos-1
    %Cogemos el minimo de un latido y el siguiente y sacamos la media del
    %valor de la señal entre esos puntos

    %Prueba 
    % red_dc(i) = mean(MinsR.Value(i):MinsR.Value(i+1));
    % ir_dc(i) = mean(MinsIr.Value(i):MinsIr.Value(i+1));

    % red_dc(i) = mean(noche.Red(find(noche.Time==MinsR.Time(i)):find(noche.Time==MinsR.Time(i+1))));
    % ir_dc(i) = mean(noche.IR(find(noche.Time==MinsIr.Time(i)):find(noche.Time==MinsIr.Time(i+1))));

    %Media 2
    % idxStartR = knnsearch(noche.Time(idxStartR:idxStartR+2560),MinsR.Time(i));
    % idxEndR = knnsearch(noche.Time(idxEndR:idxEndR+2560), MinsR.Time(i+1));
    % red_dc(i) = mean(noche.Red(idxStartR:idxEndR));
    % 
    % idxStartIR = knnsearch(noche.Time(idxStartIR:idxStartIR+2560), MinsIr.Time(i));
    % idxEndIR = knnsearch(noche.Time(idxEndIR:idxEndIR+2560), MinsIr.Time(i+1));
    % ir_dc(i) = mean(noche.IR(idxStartIR:idxEndIR));
    
    if ~isnan(MinsIr.Time(i)) && ~isnan(MinsR.Time(i)) && ~isnan(MinsIr.Time(i+1)) && ~isnan(MinsR.Time(i+1))
        [~,izda] = ismember(MinsIr.Time(i),noche.Time);
        [~,dcha] = ismember(MinsIr.Time(i+1),noche.Time);
        ir_dc(i)=mean(noche.IR(izda:dcha));
    
        [~,izda] = ismember(MinsR.Time(i),noche.Time);
        [~,dcha] = ismember(MinsR.Time(i+1),noche.Time);
        red_dc(i)=mean(noche.Red(izda:dcha));
    end
    
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
%referencia_noche=referencia(inicio_estudio:fin_estudio-1,1);
referencia_noche=referencia(1:end);
segundos = 1:1:length(referencia);
segundos=segundos'+crop*60;
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
% plot(MaxsIr.Time,MaxsIr.Value,'DisplayName', 'SPO2');
% plot(MaxsIr.Time,MaxsR.Value,'DisplayName', 'SPO2');
% plot(MaxsIr.Time,MinsR.Value,'DisplayName', 'SPO2');
% plot(MaxsIr.Time,MinsIr.Value,'DisplayName', 'SPO2');
%plot(prueba1(1:end-1),SPO2,'DisplayName', 'SPO2');
%plot(MaxsR.Time(1:end-1),-aR,'DisplayName', 'R');
% ylim([111 112])
%plot(MaxsIr.Time(1:end-1)-600+40,-aRFilt,'DisplayName', 'R apnea low-15');%El -560 esta puesto a ojo
plot(MaxsR.Time(1:end-1),-aRFilt,'DisplayName', 'R apnea');
legend

%% XCORR

% realizamos el calculo de la media de saturación por seg creando primero
% una tabla para facilitar los cálculos
mediaR=table(-aRFilt(1:end), MaxsIr.Time(41:end-1));



%Normalizar
maxSat=max(ref_noche.referencia_noche);
maxR=max(-aRFilt);
newR=-(aRFilt)/maxR;
newSatRef=ref_noche.referencia_noche/maxSat;

%hacer la media de la saturacion de los latidos por segundo para reducir el
%numero de puntos que tenemos y que sea factible hacer la correlación y
%facilitarlo para ver. Es decir, para los segundos 44-45 coger los latidos
%en ese tiempo y obtener la media de saturación indicada por la R

figure
hold on
plot(ref_noche.segundos(1:10717), newSatRef(1:10717), 'DisplayName', 'Sat Referencia');
yyaxis right
plot(MaxsIr.Time(1:10717)-600+40,newR(1:10717),'DisplayName', 'R apnea low-15');
% 
% [correlacion, lags] = xcorr(newR, newSatRef);
% stem(lags,correlacion)


%% Borrador: preprocesado: filtrado y deteccion de artefactos

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

