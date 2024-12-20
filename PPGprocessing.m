
clear 
close all force
clc

addpath(genpath('C:\Users\Pablo\Desktop\Calculos\biomedical-signal-processing'));
%cargaBio;

%% Carga y representacion de los datos
filename = 'Paciente7';
crop = 10; %Cuantos minutos ignoraremos del inicio y el final
signals = parquetread(['C:\Users\Pablo\Desktop\Calculos\' filename 'LH.gzip']);signals = signals(3:end,:);
referencia = parquetread(['C:\Users\Pablo\Desktop\Calculos\' filename '_Pulso_SpO2.parquet']);

%Limpiamos puntos incoherentes (60<SPO2<100)

for i=1:length(referencia.Saturacion)
    if (referencia.Saturacion(i)> 100)
        referencia.Saturacion(i)=nan;
    end
end

%Eliminamos los CROP primeros minutos y los CROP últimos minutos
referencia=referencia.Saturacion(crop*60:end-crop*60);
%sumamos los lugares donde sea menor que el tiempo buscado para hallar el
%índice rápido
signals=signals(sum(signals.TimeTrial<crop*60)+1:sum(signals.TimeTrial<signals.TimeTrial(end)-crop*60)-1,:);

clear i

%% Resamplear
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
    Time=signals.TimeTrial(huecos(i)):1/fo:signals.TimeTrial(huecos(i)+1);
    Time=Time';
    [n,~]=size(Time);
    GC=nan(n,1);
    RC=nan(n,1);
    IC=nan(n,1);
    resampled = [resampled;table(Time,GC,RC,IC)];

    if (height(signals)-huecos(i))>120*256%*
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

Setup.plotflag = false;
isArtifactG      =   energyArtifacts ( PPGGreen , fo , Setup );
isArtifactIR     =   energyArtifacts ( PPGIR , fo , Setup );
isArtifactRed    =   energyArtifacts ( PPGRed , fo , Setup );

PPGIR(isArtifactIR)=nan; 
PPGGreen(isArtifactG)=nan;
PPGRed(isArtifactRed)=nan;

%Obtenemos el instante de tiempo de 
% los maximos y los minimos en la señal filtrada
Setup.plotflag = false;

%Instantes de tiempo de: D->Latido, A->Max, B->Min 
%(partiendo desde el instante de tiempo de medicion inicial)
[ gD , gA , gB ]    =   pulseDelineation ( PPGGreen , fo , Setup );
[ iD , iA , iB ]    =   pulseDelineation ( PPGIR , fo , Setup );
[ nD , nA , nB ]    =   pulseDelineation ( PPGRed , fo , Setup ); 


%creamos un eje temporal que se ajuste a nuestra señal
t = 0:1/fo:(length(PPGRed)-1)/fo;

%Con tiempos arreglados:
irA=iA+remuestreado.Time(1); irB=iB+remuestreado.Time(1);
rA=nA+remuestreado.Time(1); rB=nB+remuestreado.Time(1);
grA=gA+remuestreado.Time(1); grB=gB+remuestreado.Time(1);

clear k p z sos filterGain

%% Referencia Infra para obtencion maximos rojo 
MaxRojo=zeros(size(rA)); MaxInfra=zeros(size(irA));
MinRojo=zeros(size(rB)); MinInfra=zeros(size(irB));

umbral=0.1;

i=1;
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

%Referencia rA
% while  i <= recorrer
%     inf=rA(i)-umbral;
%     sup=rA(i)+umbral;
%     [idx,]=find(irA(1:end)>inf & irA(1:end)<sup,1);
%     if ~isempty(idx)
%         MaxRojo(i)=rA(i);
%         MaxInfra(i)=irA(idx);
%     else
%         MaxRojo(i)=nan;
%         MaxInfra(i)=nan;
%     end
%     i=i+1;
% end

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

%Referencia rB
% while i <= recorrer
%     inf=rB(i)-umbral;
%     sup=rB(i)+umbral;
%     [idx,]=find(irB(1:end)>inf & irB(1:end)<sup,1);
%     if ~isempty(idx)
%         MinRojo(i)=rB(i);
%         MinInfra(i)=irB(idx);
%     else
%         MinRojo(i)=nan;
%         MinInfra(i)=nan;
%     end
%     i=i+1;
% end

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
% fc = 15/(fo/2);
% [z,p,k]         =   cheby2( ord , 20 , fc , 'low' ); 
% [sos,filterGain] = zp2sos(z,p,k);
FiltroPBajoRoja = nanfiltfilt( sos , filterGain , -remuestreado.RC(:) , Setup ); 

clear k p z sos filterGain



%% Toda la noche (creacion de las tablas)
noche=table(remuestreado.Time, FiltroPBajoRoja, FiltroPBajoIR,'VariableNames', {'Time', 'Red', 'IR'});

MaxsIr=table(irA, nan(length(irA),1), 'VariableNames', {'Time', 'Value'});
MinsIr=table(irB, nan(length(irB),1), 'VariableNames', {'Time', 'Value'});
MaxsR=table(rA, nan(length(rA),1), 'VariableNames', {'Time', 'Value'});
MinsR=table(rB, nan(length(rB),1), 'VariableNames', {'Time', 'Value'});

%% Obtención de valores con busqueda knn:

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

%% Cálculo de R y SPO2

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
end
toc
%% Cálculo FFT y limpieza de la señal
% fourierApn=fft(aR);
% aL=length(aR);
% 
% A2=abs(fourierApn/aL);
% A1=A2(1:aL/2+1);
% A1(2:end-1)= 2*A1(2:end-1);
% fApn=fo/aL*(0:(aL/2));

%En vistas de las graficas nos interesa establecer el corte en 20Hz, ya que
%es el punto suficientemente alejado como para no perder precisión, pero
%suficientemente cercano para limpiar bien la señal.
fc = 20/(fo/2); 
[z,p,k]         =   cheby2( ord , 20 , fc , 'low' ); 
[sos,filterGain] = zp2sos(z,p,k);
aRFilt = nanfiltfilt( sos, filterGain , -aR, Setup );

SPO2=zeros(1,length(aR));
SPO2Filt= zeros(1,length(aR));

for i=1:length(SPO2)
    SPO2Filt(i)=(110-aRFilt(i));
    SPO2(i)=(110-aR(i));
end

%% Dibujo valores de aR
inicio_estudio=floor(noche.Time(1));
fin_estudio=ceil(noche.Time(end));
referencia_noche=referencia(1:end);
segundos = 1:1:length(referencia);
segundos=segundos'+crop*60;
ref_noche=table(referencia_noche, segundos, 'VariableNames', {'referencia_noche', 'segundos'});

figure
hold on
plot(seconds(ref_noche.segundos), ref_noche.referencia_noche, 'DisplayName', 'Sat Referencia');
ylabel('Oxygen Saturation (%)');
yyaxis right
plot(seconds(MaxsR.Time(1:end-1)),-aRFilt,'DisplayName', 'aR apnea');
ylabel('aR');
xtickformat('hh:mm:ss')
xlabel('Tiempo (hh:mm:ss)');
legend

%% Busqueda valores con los que sincronizar
% realizamos el calculo de la media de saturación por seg
k=1;
condicion=1;
lista=floor(MinsR.Time);
mediaLat=mean(aRFilt(1));
while condicion
    if isnan(MinsR.Time(end-k))
        k=k+1;
    else
        condicion=0;
    end
end

for i=2:floor(MinsR.Time(end-k)-crop*60)
    ventana=lista>=i+(crop*60) & lista<(i+1+(crop*60));
    if any(ventana) 
            mediaLat=[mediaLat,mean(aRFilt(ventana))];
    else
            mediaLat=[mediaLat,nan(1,1)];
    end
end

mediaLat=mediaLat';
% figure;boxplot(mediaLat);
 
%Limites puestos a ojo
if filename=="Paciente1"
    for i=1:height(mediaLat)
        if mediaLat(i)>-1.1 || mediaLat(i)<=-2
            mediaLat(i)=nan;
        end
    end
elseif filename=="Paciente3"
    for i=1:height(mediaLat)
        if mediaLat(i)>-0.75 || mediaLat(i)<=-3.25
            mediaLat(i)=nan;
        end
    end
elseif filename=="Paciente6"
    for i=1:height(mediaLat)
        if mediaLat(i)>1.5 || mediaLat(i)<=-3
            mediaLat(i)=nan;
        end
    end
elseif filename=="Paciente7"
    for i=1:height(mediaLat)
        if mediaLat(i)>1.5 || mediaLat(i)<=-2.50
            mediaLat(i)=nan;
        end
    end
elseif filename=="Paciente9"
    for i=1:height(mediaLat)
        if mediaLat(i)>0.5 || mediaLat(i)<=-2
            mediaLat(i)=nan;
        end
    end
elseif filename=="Paciente10"
    for i=1:height(mediaLat)
        if mediaLat(i)>6 || mediaLat(i)<=-13
            mediaLat(i)=nan;
        end
    end
elseif filename=="Paciente11"
    for i=1:height(mediaLat)
        if mediaLat(i)>10|| mediaLat(i)<=-10
            mediaLat(i)=nan;
        end
    end
end

mediaLat=fillmissing(mediaLat,'linear');
nonanSatRef=fillmissing(ref_noche.referencia_noche,'linear');

if filename=="Paciente1"
    for i=1:height(nonanSatRef)
        if nonanSatRef(i)<70
            nonanSatRef(i)=nan;
        end
        if nonanSatRef(i)>100
            nonanSatRef(i)=nan;
        end
    end
    nonanSatRef=fillmissing(nonanSatRef,'linear');
end

[correl, lag] = xcorr(nonanSatRef, -mediaLat,6000);
% figure
% stem(lag,correl)

x=nonanSatRef;
y=mediaLat;

%Añadidos para sincronizar las señales
if filename=="Paciente1"
    padd = 88; 
elseif filename=="Paciente3"
    padd = 230;
elseif filename=="Paciente6"
    padd = 0;
elseif filename=="Paciente7"
    padd = 188;
elseif filename=="Paciente9"
    padd = 880;
elseif filename=="Paciente10"
    padd = 451;
elseif filename=="Paciente11"
    padd = 200;
end

x = [zeros(padd,1)+mean(x);x];

min_length = min(length(x),length(y));
x=x(1:min_length);
y=y(1:min_length);

nonanSatRef=x;
mediaLat=y;

clear x y

%% Cambio signo mediaLat para ponerlo todo en orden
mediaLat=-mediaLat;

%% Dibujo de la señal sincronizada y cálculo de la correlación

refTemporal=seconds(1:1:size(nonanSatRef));

figure
hold on
plot(refTemporal, nonanSatRef, 'DisplayName', 'Sat Referencia');
ylabel('Oxygen Saturation (%)');
yyaxis right
plot(refTemporal, mediaLat,'DisplayName', 'aR apnea');
ylabel('aR');
xtickformat('hh:mm:ss')
xlabel('Tiempo (hh:mm:ss)');
legend;

[R,p] = corrcoef(nonanSatRef,mediaLat)

%% Aplicación del filtro de medianas
mediaLatMedFilt=medfilt1(mediaLat,10);

figure
hold on
plot(refTemporal, nonanSatRef, 'DisplayName', 'Sat Referencia');
ylabel('Oxygen Saturation (%)');
yyaxis right
plot(refTemporal, mediaLatMedFilt, 'DisplayName', 'aR apnea');
ylabel('aR');
xtickformat('hh:mm:ss')
xlabel('Tiempo (hh:mm:ss)');
legend

% [R,p] = corrcoef(nonanSatRef,mediaLatMedFilt)


%% Pseudo calculo del SPO2 con Latidos calculados con la mediana
c1= -16.666666;
c2= 8.333333;
c3= 100;
pseudo=-mediaLatMedFilt;    %Devolvemos al valor original que es el  
                            %esperado para realizar el cálculo
PseudoSPO2= -(c1*pseudo.^2+c2*pseudo+c3); 

PseudoSPO2Norm = (PseudoSPO2-min(PseudoSPO2))/(max(PseudoSPO2)-min(PseudoSPO2))*(max(nonanSatRef)-min(nonanSatRef));
PseudoSPO2Norm = PseudoSPO2Norm-mean(PseudoSPO2Norm)+mean(nonanSatRef);

figure
hold on
plot(refTemporal, nonanSatRef, 'DisplayName', 'Referencia');
ylabel('Oxygen Saturation (%)');
ylim([86 100])
yticks(86:1:100);
yyaxis right 
plot(refTemporal, PseudoSPO2Norm, 'DisplayName', 'Estimada');
ylim([86 100])
yticks(86:1:100);
ylabel('Oxygen Saturation (%)');
xtickformat('hh:mm:ss')
xlabel('Tiempo (hh:mm:ss)');
legend

size=min(height(nonanSatRef), height(PseudoSPO2Norm));

[R, p]=corrcoef(nonanSatRef(1:size), PseudoSPO2Norm(1:size))


errAbs = (PseudoSPO2Norm(1:height(nonanSatRef)))-nonanSatRef;
disp("ABS: " + mean(abs(errAbs)))

errRel = abs(errAbs)./nonanSatRef;
disp("REL: " + mean(errRel))


%% Normalización (Hecho con la mediana)
LatNorm=(mediaLatMedFilt-min(mediaLatMedFilt))/(max(mediaLatMedFilt)-min(mediaLatMedFilt))*(max(nonanSatRef)-min(nonanSatRef));
LatNorm=LatNorm+mean(nonanSatRef)-mean(LatNorm);

% errAbs = LatNorm(1:height(nonanSatRef))-nonanSatRef;
% disp("ABS: " + mean(abs(errAbs)))
% 
% errRel = abs(errAbs)./nonanSatRef;
% disp("REL: " + mean(errRel))
% 
% [R,p]=corrcoef(nonanSatRef, LatNorm(1:height(nonanSatRef)))

figure
hold on
plot(refTemporal, nonanSatRef, 'DisplayName', 'Referencia');
ylabel('Oxygen Saturation (%)');
ylim([86 100])
yticks(86:1:100);
yyaxis right
plot(refTemporal, LatNorm,'DisplayName', 'Estimada');
ylim([86 100])
yticks(86:1:100);
ylabel('Oxygen Saturation (%)');
xtickformat('hh:mm:ss')
xlabel('Tiempo (hh:mm:ss)');
legend

%% Filtrado de la normalización
fc=0.05;
ord=1;
[Bb,Ba] = butter(ord, fc/(1/2), 'low');
FilteredCalc=filtfilt(Bb,Ba, LatNorm);

[R,p]=corrcoef(nonanSatRef, FilteredCalc(1:height(nonanSatRef)))

errAbs = FilteredCalc(1:height(nonanSatRef))-nonanSatRef;
disp("ABS: " + mean(abs(errAbs)))

errRel = abs(errAbs)./nonanSatRef;
disp("REL: " + mean(errRel))

figure
hold on
plot(refTemporal, nonanSatRef, 'DisplayName', 'Referencia');
ylabel('Oxygen Saturation (%)');
ylim([82 100])
yticks(82:1:100);
plot(refTemporal, FilteredCalc,'DisplayName', 'Estimada');
xtickformat('hh:mm:ss')
xlabel('Tiempo (hh:mm:ss)');
legend
