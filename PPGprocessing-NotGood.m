
clear 
close all force
clc

addpath(genpath('C:\Users\Pablo\Desktop\Calculos\biomedical-signal-processing'));
%cargaBio;

%% carga y representacion de los datos
filename = 'Paciente1';
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


Setup.wdwVariance = 25;
Setup.C = 3;
Setup.plotflag = false;
isArtifactG      =   energyArtifacts ( PPGGreen , fo , Setup );
isArtifactIR     =   energyArtifacts ( PPGIR , fo , Setup );
isArtifactRed    =   energyArtifacts ( PPGRed , fo , Setup );

% clear Setup
% Setup.plotflag  =	true;
% isArtifactG      =   energyArtifacts ( nanzscore(PPGGreen) , fo , Setup );

PPGIR(isArtifactIR)=nan; 
PPGGreen(isArtifactG)=nan;
PPGRed(isArtifactRed)=nan;

%Obtenemos el instante de tiempo de 
% los maximos y los minimos en la señal filtrada
Setup.plotflag = false;

%Instantes de tiempo de: D->Latido, A->Max, B->Min 
%(partiendo desde el instante de tiempo de medicion inicial)
%Es decir, que son isntantes desde que inicia la tabla válida
Setup.plotflag = false;
[ gD , gA , gB ]    =   pulseDelineation ( PPGGreen , fo , Setup );
[ iD , iA , iB ]    =   pulseDelineation ( PPGIR , fo , Setup );
[ nD , nA , nB ]    =   pulseDelineation ( PPGRed , fo , Setup ); 

clear Setup
Setup.tol        =	1.5;
% Setup.plotflag	=	true ;
[ ~, ~ , ~ , tn_iA , cn ]	=	incidences ( iA , Setup );
tn_iA(cn=='X')=[];
[ ~, ~ , ~ , tn_iB , cn ]	=	incidences ( iB , Setup );
tn_iB(cn=='X')=[];
[ ~, ~ , ~ , tn_nA, cn ]	=	incidences ( nA , Setup );
tn_nA(cn=='X')=[];
[ ~, ~ , ~ , tn_nB , cn ]	=	incidences ( nB , Setup );
tn_nB(cn=='X')=[];
% pruebaTN=ismember(nB,tn_nB);
% tn_nB = 
% disp([ '% correct detections: ' num2str(round(mean(cn=='c')*100)) '%'])


% % % % tn_gA(cn=='X')=NaN;
% % % open computeHRVsignals.m
% % % [mdt , tnCorrected] =	correctModulatingSignal ( mdt , fa , tn , cn , tRR , T );
% % % T = 2;
% % % tn_inventado = diff(gA)<= 1/T;


%creamos un eje temporal que se ajuste a nuestra señal
t = 0:1/fo:(length(PPGRed)-1)/fo;

% a_iA=iA; a_iB=iB; %para parte de apnea
% a_nA=nA; a_nB=nB;
% 
% n_iA=iA; n_iB=iB; %para parte noche normal
% n_nA=nA; n_nB=nB;  

%Con tiempos arreglados:
% irA=iA+remuestreado.Time(1); irB=iB+remuestreado.Time(1);
% rA=nA+remuestreado.Time(1); rB=nB+remuestreado.Time(1);
% grA=gA+remuestreado.Time(1); grB=gB+remuestreado.Time(1);

irA=tn_iA+remuestreado.Time(1); irB=tn_iB+remuestreado.Time(1);
rA=tn_nA+remuestreado.Time(1); rB=tn_nB+remuestreado.Time(1);
grA=gA+remuestreado.Time(1); grB=gB+remuestreado.Time(1);

clear k p z sos filterGain

%% Referencia Infra para obtencion maximos rojo (PROBAR al reves a ver cual da mejores resultados)
MaxRojo=zeros(size(rA)); MaxInfra=zeros(size(irA));
MinRojo=zeros(size(rB)); MinInfra=zeros(size(irB));

%Prueba modificación
figure
hold on
plot(remuestreado.Time, remuestreado.RC, 'DisplayName', 'Señal');
yyaxis right
plot(rA, 'DisplayName', 'rA');
legend

recorrer=min(length(rA), length(irA));
rA=rA(1:recorrer);
irA=irA(1:recorrer);

i=1;
umbral=0.1;
%referencia irA
% while  i <= recorrer
%     inf=irA(i)-umbral;
%     sup=irA(i)+umbral;
%     [idx,]=find(rA(1:end)>inf & rA(1:end)<sup,1);
%     if ~isempty(idx)
%         MaxRojo(i)=rA(idx);
%         MaxInfra(i)=irA(i);
%     else
%         MaxRojo(i)=nan;
%         MaxInfra(i)=nan;
%     end
%     i=i+1;
% end

%Prueba referencia rA
while  i <= recorrer
    inf=rA(i)-umbral;
    sup=rA(i)+umbral;
    [idx,]=find(irA(1:end)>inf & irA(1:end)<sup,1);
    if ~isempty(idx)
        MaxRojo(i)=rA(i);
        MaxInfra(i)=irA(idx);
    else
        MaxRojo(i)=nan;
        MaxInfra(i)=nan;
    end
    i=i+1;
end

plot(MaxRojo, 'DisplayName', 'Recalculado');
legend

figure
hold on
% plot(remuestreado.Time, remuestreado.RC, 'DisplayName', 'Señal');
% yyaxis right
plot(MaxRojo, 'DisplayName', 'Recalculado');
legend




recorrer=min(length(rB), length(irB));
rB=rB(1:recorrer);
irB=irB(1:recorrer);

i=1;
%Referencia irB
% while i <= recorrer
%     inf=irB(i)-umbral;
%     sup=irB(i)+umbral;
%     [idx,]=find(rB(1:end)>inf & rB(1:end)<sup,1);
%     if ~isempty(idx)
%         MinRojo(i)=rB(idx);
%         MinInfra(i)=irB(i);
%     else
%         MinRojo(i)=nan;
%         MinInfra(i)=nan;
%     end
%     i=i+1;
% end

%Referencia rB
while i <= recorrer
    inf=rB(i)-umbral;
    sup=rB(i)+umbral;
    [idx,]=find(irB(1:end)>inf & irB(1:end)<sup,1);
    if ~isempty(idx)
        MinRojo(i)=rB(i);
        MinInfra(i)=irB(idx);
    else
        MinRojo(i)=nan;
        MinInfra(i)=nan;
    end
    i=i+1;
end

% % Quitar 0 al final para dejar ambos indices del mismo tamaño
% last_nonzero_idx = find(MaxRojo ~= 0, 1, 'last');
% if ~isempty(last_nonzero_idx)
%     MaxRojo = MaxRojo(1:last_nonzero_idx);
% end
% last_nonzero_idx = find(MaxInfra ~= 0, 1, 'last');
% if ~isempty(last_nonzero_idx)
%     MaxInfra = MaxInfra(1:last_nonzero_idx);
% end
% last_nonzero_idx = find(MinRojo ~= 0, 1, 'last');
% if ~isempty(last_nonzero_idx)
%     MinRojo = MinRojo(1:last_nonzero_idx);
% end
% last_nonzero_idx = find(MinInfra ~= 0, 1, 'last');
% if ~isempty(last_nonzero_idx)
%     MinInfra = MinInfra(1:last_nonzero_idx);
% end

irA=MaxInfra;rA=MaxRojo;irB=MinInfra;rB=MinRojo;

%% Creamos nuevo filtro para mantener la respiración
%IR
fc = 15/(fo/2);
[z,p,k]         =   cheby2( ord , 20 , fc , 'low' ); 
[sos,filterGain] = zp2sos(z,p,k);
%Este nuevo filtro solo quta el ruido, deja la respiración
Setup.plotflag = false;
FiltroPBajoIR = nanfiltfilt( sos , filterGain , -remuestreado.IC(:) , Setup );


%Roja
fc = 15/(fo/2);
[z,p,k]         =   cheby2( ord , 20 , fc , 'low' ); 
[sos,filterGain] = zp2sos(z,p,k);
Setup.plotflag = false;
FiltroPBajoRoja = nanfiltfilt( sos , filterGain , -remuestreado.RC(:) , Setup ); 

clear Setup
Setup.wdwVariance = 25;
Setup.C = 3;
Setup.plotflag = true;
isArtifactIR     =   energyArtifacts ( FiltroPBajoIR , fo , Setup );
isArtifactRed    =   energyArtifacts ( FiltroPBajoRoja , fo , Setup );

FiltroPBajoIR(isArtifactIR)=nan; 
FiltroPBajoRoja(isArtifactRed)=nan;

% clear Setup
% Setup.tol        =	1.5;
% Setup.plotflag	=	true ;
% [ ~, ~ , ~ , AFiltroPBajoIR , cn ]	=	incidences ( FiltroPBajoIR , Setup );
% AFiltroPBajoIR(cn=='X')=[];
% [ ~, ~ , ~ , AFiltroPBajoRoja , cn ]	=	incidences ( FiltroPBajoRoja , Setup );
% AFiltroPBajoRoja(cn=='X')=[];

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

figure
hold on
% plot(MinsR.Time(8750:8760), 1+zeros(11, 1), "ro", 'DisplayName', 'MinR')
plot(MaxsR.Time(8650:8850), 1.1+zeros(201, 1), "r+", 'DisplayName', 'MaxsR') %Referencia en el fallo de i=8755
% plot(MinsIr.Time, 1.2+zeros(height(MinsIr),1), "bo", 'DisplayName', 'MinsIr')
% plot(MaxsIr.Time, 1.3+zeros(height(MaxsIr),1), "b+", 'DisplayName', 'MaxsIr')
legend

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

figure
hold on
plot(MinsR.Time, 1+zeros(height(MinsR),1), "ro", 'DisplayName', 'MinR')
plot(MaxsR.Time, 1.1+zeros(height(MaxsR),1), "r+", 'DisplayName', 'MaxsR')
plot(MinsIr.Time, 1.2+zeros(height(MinsIr),1), "bo", 'DisplayName', 'MinsIr')
plot(MaxsIr.Time, 1.3+zeros(height(MaxsIr),1), "b+", 'DisplayName', 'MaxsIr')
legend
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
% plot(tn_nB+crop*60,zeros(height(tn_nB),1)+1, "og", 'DisplayName', 'tn_nB');
plot(MaxsR.Time, MaxsR.Value, "r*", 'DisplayName', 'MAX');
plot(MinsR.Time, MinsR.Value, "bo", 'DisplayName', 'MINS');
legend

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
hold on
plot(ref_noche.segundos, ref_noche.referencia_noche, 'DisplayName', 'Sat Referencia');
yyaxis right
% plot(MaxsIr.Time,MaxsIr.Value,'DisplayName', 'SPO2');
% plot(MaxsIr.Time,MaxsR.Value,'DisplayName', 'SPO2');
% plot(MaxsIr.Time,MinsR.Value,'DisplayName', 'SPO2');
% plot(MaxsIr.Time,MinsIr.Value,'DisplayName', 'SPO2');
% plot(prueba1(1:end-1),SPO2','DisplayName', 'SPO2');
plot(MaxsR.Time(1:end-1),-aR,'DisplayName', 'R');
% ylim([111 112])
% plot(MaxsR.Time(1:end-1), MedSatCalc,'DisplayName', 'R apnea low-15');
% plot(MaxsR.Time(1:end-1),-aRFilt,'DisplayName', 'R apnea');
% plot(crop*60:size(mediaLat)+crop*60-1, -mediaLat','DisplayName', 'R apnea');
% plot(crop*60:size(mediaLat)+crop*60-1, -mediaLat','DisplayName', 'R apnea');
legend

%% XCORR
% realizamos el calculo de la media de saturación por seg
lista=floor(MinsR.Time);
mediaLat=mean(aRFilt(1));
for i=2:floor(MinsR.Time(end-3)-crop*60)
    ventana=lista>=i+(crop*60) & lista<(i+1+(crop*60));
    if any(ventana) 
            mediaLat=[mediaLat,mean(aRFilt(ventana))];
    else
             mediaLat=[mediaLat,nan(1,1)];
    end
end

mediaLat=mediaLat';

%Limites puestos a ojo
if filename=="Paciente0"
    for i=1:height(mediaLat)
        if mediaLat(i)>-1.1 || mediaLat(i)<=-2.5
            mediaLat(i)=nan;
        end
    end
elseif filename=="Paciente1"
    for i=1:height(mediaLat)
        if mediaLat(i)>-1.1 || mediaLat(i)<=-2.5
            mediaLat(i)=nan;
        end
    end
elseif filename=="Paciente3"
    for i=1:height(mediaLat)
        if mediaLat(i)>-0.75 || mediaLat(i)<=-3.25
            mediaLat(i)=nan;
        end
    end
end
    
mediaLat=fillmissing(mediaLat,'linear');
nonanSatRef=fillmissing(ref_noche.referencia_noche,'linear');

[correl, lag] = xcorr(nonanSatRef, -mediaLat,6000);

figure
stem(lag,correl)

x=nonanSatRef;
y=mediaLat;

%Si el desplazamiento es positivo se realiza sobre el gráfico x, si es 
% negativo sobre el gáfico y
if filename=="Paciente0"
    %Desplazamiento de 0    delay = 1;
    x=x(delay:end);
    x=circshift(nonanSatRef ,-delay); nonanSatRef(end-delay:end) = [];
elseif filename=="Paciente1"
    delay = 88;
    y=y(delay:end);
    y=circshift(mediaLat ,-delay); y(end-delay:end) = [];

elseif filename=="Paciente3"
    delay = 115;
    y=y(delay:end);
    y=circshift(mediaLat ,-delay); y(end-delay:end) = [];
end

min_length = min(length(x),length(y));
x=x(1:min_length);
y=y(1:min_length);
   
mediaLat=medfilt1(y,25);

figure
plot(x);yyaxis right
% plot(-y);
hold on
plot(-mediaLat);
    
[R,p] = corrcoef(x,mediaLat)


%% Pseudo calculo del SPO2 con Latidos calculados con la mediana
c1= -16.666666;
c2= 8.333333;
c3= 100;
PseudoSPO2= c1*mediaLat.^2+c2*mediaLat+c3; 

for i=1:height(PseudoSPO2)
    if PseudoSPO2(i)<40
        PseudoSPO2(i)=nan;
    end
    if PseudoSPO2(i)>70
        PseudoSPO2(i)=nan;
    end
end

PseudoSPO2=fillmissing(-PseudoSPO2,'linear');

PseudoSPO2Norm = (PseudoSPO2-min(PseudoSPO2))/(max(PseudoSPO2)-min(PseudoSPO2))*(max(x)-min(x));
PseudoSPO2Norm = PseudoSPO2Norm-max(PseudoSPO2Norm)+max(x);

%El añadido es la media de saturación sin nan
if filename=="Paciente0"
    ;
elseif filename=="Paciente1"
    x2 = [zeros(88,1)+mean(x);x];
elseif filename=="Paciente3"
    x2 = [zeros(115,1)+mean(x);x];
end

figure
plot(x2);
hold on
plot(round(PseudoSPO2Norm));

size=min(height(x2), height(PseudoSPO2Norm));

[R, p]=corrcoef(x2(1:size), round(PseudoSPO2Norm(1:size)))

%% Normalización
% for i=1:height(mediaLat)
%     if mediaLat(i)<-1.8
%         mediaLat(i)=nan;
%     end
% end
% mediaLat=fillmissing(-mediaLat,'linear');

for i=1:height(x)
    if x(i)<85
        x(i)=nan;
    end
    if x(i)>100
        x(i)=nan;
    end
end
x2=fillmissing(x,'linear');
if filename=="Paciente0"
    ;
elseif filename=="Paciente1"
    x2 = [zeros(88,1)+mean(x2);x2];
elseif filename=="Paciente3"
    x2 = [zeros(115,1)+mean(x2);x2];
end

LatNorm=(mediaLat-min(mediaLat))/(max(mediaLat)-min(mediaLat))*(max(x2)-min(x2));
LatNorm=LatNorm+mean(x2)-mean(LatNorm);

% errorR=abs(latNorm(1:height(x))-x)/x;
errAbs = LatNorm(1:height(x2))-x2;
disp("ABS: " + mean(abs(errAbs)))

errRel = abs(errAbs)./x2;
disp("REL: " + mean(errRel))

fc=0.15;
ord=1;

[Bb,Ba] = butter(ord, fc/(1/2), 'low');
FilteredRef=filtfilt(Bb,Ba, LatNorm);

[R,p]=corrcoef(x2, FilteredRef(1:height(x2)))


figure
plot(LatNorm,'DisplayName', 'Media de latidos por segundo normalizada');
hold on
plot(x2, 'DisplayName', 'Sat Referencia');
yyaxis right
plot(errAbs, 'DisplayName', 'Error Absoluto')
plot(errRel, 'DisplayName', 'Error Relativo')
legend

%% DTW test
% Compute DTW alignment between the two sequences
[dist, ix, iy] = dtw(nonanSatRef, -mediaLat);

% ix and iy give the indices that match HR_A to HR_B
aligned_A = nonanSatRef(ix);
aligned_B = -mediaLat(iy);

figure
plot(aligned_A);
hold on;
yyaxis right
plot(aligned_B);
legend('Aligned REF', 'Aligned CALC-SAT');
title('Aligned SAT TO REF');

corrcoef(aligned_B, aligned_A)
[A,p] = corrcoef(aligned_B, aligned_A)

%¿APLICAR?
% maxSat=max(ref_noche.referencia_noche);
% maxR=max(-aRFilt);
% newR=-(aRFilt)/maxR;
% newSatRef=ref_noche.referencia_noche/maxSat;

% [correlacion, lags] = xcorr(newR, newSatRef);
% stem(lags,correlacion)

%% Atenuacion de saturación de referencia
for i=1:height(nonanSatRef)
    if nonanSatRef(i)<85
        nonanSatRef(i)=nan;
    end
    if nonanSatRef(i)>100
        nonanSatRef(i)=nan;
    end
end
nonanSatRef=fillmissing(nonanSatRef,'linear');

fc=0.15;
ord=1;

[Bb,Ba] = butter(ord, fc/(1/2), 'low');
FilteredRef=filtfilt(Bb,Ba, nonanSatRef);
% for i=1:height(FilteredRef)
%     FilteredRef(i)=round(FilteredRef(i));
% end

figure
hold on
plot(nonanSatRef, 'DisplayName', 'Señal Original');
plot(FilteredRef,'DisplayName', 'Señal Filtrada');
yyaxis right;
plot(-y2,'DisplayName', 'SPO2 calculada mediana');
legend;

[A,p] = corrcoef(y2(1:height(FilteredRef)), FilteredRef)

clear ord; Ba; Bb; t; 



