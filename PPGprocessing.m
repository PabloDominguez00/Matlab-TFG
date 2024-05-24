
clear 
close all force
clc

%addpath(genpath('C:\Users\Pablo\Desktop\Calculos\biomedical-signal-processing'));
%Encintramos los saltos: signals.TimeTrial(diff(signals.TimeTrial)>0.5)
%(Estan redondeados a precisión de segundos)
%Tarea: interpolar la frecuencia de muestreo a 256, desde la frecuencia
%original (calculada con los huecos dejados con la formula anterior)
%% carga y representacion de los datos
filename = 'Paciente2';
signals = parquetread(['C:\Users\Pablo\Desktop\Calculos\' filename 'LH.gzip']);
referencia = struct2table(load('C:\Users\Pablo\Desktop\Calculos\Referencia.mat'));

%copy = signals(:,1:4);

%Limpiamos puntos incoherentes (80<SPO2<100)
for i=1:length(referencia.ans)
    if (referencia.ans(i)> 100)
        referencia.ans(i)=nan;
    elseif (referencia.ans(i) < 80)
        referencia.ans(i)=nan;
    end
end

clear i

%% remuestreo a 256Hz de los diferentes apartados
fo=256;
%Encontramos los saltos: 
fs=1:1:7;
huecos = find(diff(signals.TimeTrial)>0.5);
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

resampled = {};

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
        GC = resample(double(g),ffin*100,fix(round(fi*100,2)));
        RC = resample(double(r),ffin*100,fix(round(fi*100,2)));
        IC = resample(double(ir),ffin*100,fix(round(fi*100,2)));
        t = ((huecos(i-1)+1)/fo):1/fo:(huecos(i-1)+length(GC))/256;

        resampled{end+1}=table(t',GC,RC,IC, 'VariableNames', {'Time', 'GC', 'RC', 'IC'});
        
        %figure
        %hold on
        %plot(t,GC,'m')
        %plot(t,RC,'c')
        %plot(t,IC, 'y')
        %plot(signals.TimeTrial(huecos(i-1)+1:huecos(i))+17.5,signals.Green_Count(huecos(i-1)+1:huecos(i)),'g')
        %plot(signals.TimeTrial(huecos(i-1)+1:huecos(i))+17.5,signals.Red_Count(huecos(i-1)+1:huecos(i)),'r')
        %plot(signals.TimeTrial(huecos(i-1)+1:huecos(i))+17.5,signals.IR_Count(huecos(i-1)+1:huecos(i)),'b')

    end
end

remuestreado =  resampled{1}; %Los intervalos a analizar se encuentran en
%segundos= 1:1:length(referencia.ans); %Referencia sin decimales

clear g r ir fi ffin i GC RC IC t fs huecos %Borramos las variable que no vamos a usar mas adelante

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

% remuestreado.GC = PPGGreen;
% remuestreado.IC = PPGIR;
% remuestreado.RC = PPGRed;

%Obtenemos el instante de tiempo de 
% los maximos y los minimos en la señal filtrada
Setup.plotflag = false;

%instantes de tiempo de: D->Latido, A->Max, B->Min
[ iD , iA , iB ]    =   pulseDelineation ( PPGIR , fo , Setup );
[ nD , nA , nB ]    =   pulseDelineation ( PPGRed , fo , Setup ); 

%creamos un eje temporal que se ajuste a nuestra señal
t = 0:1/fo:(length(PPGRed)-1)/fo;

a_iA=iA; a_iB=iB; 
a_nA=nA; a_nB=nB;

n_iA=iA; n_iB=iB; 
n_nA=nA; n_nB=nB;  

clear k p z sos filterGain
%% Infrarroja
%Creamos nuevo filtro para mantener la respiración
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
plot(iA(~isnan(iA)), FiltroPBajoIR(1+round(iA(~isnan(iA))*fo)), 'b*','LineWidth',1);
% plot(nM(~isnan(nM)), remuestreado.RC(1+round(nM(~isnan(nM))*fs)), 'k*','LineWidth',1);
plot(iB(~isnan(iB)), FiltroPBajoIR(1+round(iB(~isnan(iB))*fo)), 'b*','LineWidth',1);
%xline(iD,'color',[0.35 0.35 0.35],'LineWidth',0.5,'HandleVisibility','off');
title('PPGIR');
%Claculamos la diferencia entre el maximo y el minimo de cada latido:
%diffMaxMinIR=FiltroPBajoIR(1+round(iA(~isnan(iA))*fo))- FiltroPBajoIR(1+round(iB(~isnan(iB))*fo));

%Roja
fc = 15/(fo/2);
[z,p,k]         =   cheby2( ord , 20 , fc , 'low' ); 
[sos,filterGain] = zp2sos(z,p,k);
FiltroPBajoRoja = nanfiltfilt( sos , filterGain , -remuestreado.RC(:) , Setup ); %Este nuevo filtro solo quta el ruido, deja l respiración

figure;
hold on;
plot(t,  FiltroPBajoRoja , 'k','LineWidth',1);
plot(nD(~isnan(nD)), FiltroPBajoRoja(1+round(nD(~isnan(nD))*fo)), 'ro','LineWidth',1);
plot(nA(~isnan(nA)), FiltroPBajoRoja(1+round(nA(~isnan(nA))*fo)), 'b*','LineWidth',1);
plot(nB(~isnan(nB)), FiltroPBajoRoja(1+round(nB(~isnan(nB))*fo)), 'b*','LineWidth',1);
title('PPGRED');

clear k p z sos filterGain

%% Apnea: Obtencion de los instantes temporales donde se producen máximos y minimos en Red e IR
%Obtenemos un intervalo representativo de apnea
Red = FiltroPBajoRoja(40*60*256:60*60*256);
IR = FiltroPBajoIR(40*60*256:60*60*256);
Time=remuestreado.Time(find(remuestreado.Time==(40*60)):find(remuestreado.Time==(60*60)));

apnea=table(Time, Red, IR);
clear Red IR Time

%Claculamos el intervalo donde estan los indices con maximos y minimos de
%los latidos
%IR
maxIR = (a_iA>(40*60) & a_iA<(60*60));
maxIR = find(maxIR==1);
a_iA=a_iA(maxIR(1):maxIR(end));

minIR = (a_iB>(40*60) & a_iB<(60*60));
minIR = find(minIR==1);
a_iB=a_iB(minIR(1):minIR(end));

%Red
maxRed = (a_nA>(40*60) & a_nA<(60*60));
maxRed = find(maxRed==1);
a_nA=a_nA(maxRed(1):maxRed(end));

minRed = (a_nB>(40*60) & a_nB<(60*60));
minRed = find(minRed==1);
a_nB=a_nB(minRed(1):minRed(end));

clear minRed minIR maxIR maxRed

latidos = min(length([a_nA, a_nB, a_iA, a_iB]));
Value=nan(latidos,1);

MaxsIr=table(a_iA, Value, 'VariableNames', {'Time', 'Value'});
MinsIr=table(a_iB, Value, 'VariableNames', {'Time', 'Value'});
MaxsR=table(a_nA, Value, 'VariableNames', {'Time', 'Value'});
MinsR=table(a_nB, Value, 'VariableNames', {'Time', 'Value'});
clear Value a_nA a_nB a_iA a_iB

%% Asignacion correcta del valor de la onda en los valores de tiempo dados 
for i=1:latidos
    %Maximos de IR
    if isempty(find(apnea.Time==MaxsIr.Time(i), 1)) %aproximar al valor correcto más cercano
        %Optimizacion para la busqueda, sumamos los elementos true del
        %vector para hayar el indice que buscamos
        ant=sum(apnea.Time<MaxsIr.Time(i));
        sig=ant+1;
        if apnea.IR(ant)>apnea.IR(sig)
            MaxsIr.Time(i)=apnea.Time(ant);
            MaxsIr.Value(i)=apnea.IR(ant);
        else
            MaxsIr.Time(i)=apnea.Time(sig);
            MaxsIr.Value(i)=apnea.IR(sig);
        end
    else%tenemos valor para el tiempo dado
        MaxsIr.Value(i)=apnea.IR(find(apnea.Time==MaxsIr.Time(i)));
    end

    %Mínimos de IR
    if isempty(find(apnea.Time==MinsIr.Time(i), 1)) 
        ant=sum(apnea.Time<MinsIr.Time(i));
        sig=ant+1;
        if apnea.IR(ant)<apnea.IR(sig)
            MinsIr.Time(i)=apnea.Time(ant);
            MinsIr.Value(i)=apnea.IR(ant);
        else
            MinsIr.Time(i)=apnea.Time(sig);
            MinsIr.Value(i)=apnea.IR(sig);
        end
    else
        MinsIr.Value(i)=apnea.IR(find(apnea.Time==MinsIr.Time(i)));
    end
    
    %Maximos de Red
    if isempty(find(apnea.Time==MaxsR.Time(i), 1))
        ant=sum(apnea.Time<MaxsR.Time(i));
        sig=ant+1;
        if apnea.Red(ant)>apnea.Red(sig)
            MaxsR.Time(i)=apnea.Time(ant);
            MaxsR.Value(i)=apnea.Red(ant);
        else
            MaxsR.Time(i)=apnea.Time(sig);
            MaxsR.Value(i)=apnea.Red(sig);
        end
    else
        MaxsR.Value(i)=apnea.Red(find(apnea.Time==MaxsR.Time(i)));
    end

    %Mínimos de Red
    if isempty(find(apnea.Time==MinsR.Time(i), 1)) 
        ant=sum(apnea.Time<MinsR.Time(i));
        sig=ant+1;
        if apnea.Red(ant)<apnea.Red(sig)
            MinsR.Time(i)=apnea.Time(ant);
            MinsR.Value(i)=apnea.Red(ant);
        else
            MinsR.Time(i)=apnea.Time(sig);
            MinsR.Value(i)=apnea.Red(sig);
        end
    else
        MinsR.Value(i)=apnea.Red(find(apnea.Time==MinsR.Time(i)));
    end
end

%% Calculo de R y SPO2
%Declaramos por optimización
red_dc=zeros(latidos-1,1);
ir_dc=zeros(latidos-1,1);
red_ac=zeros(latidos-1,1);
ir_ac=zeros(latidos-1,1);
aR=zeros(latidos-1,1);

for i=1:latidos-1
    %Cogemos el minimo de un latido y el siguiente y sacamos la media del
    %valor de la señal entre esos puntos
    red_dc(i) = mean(apnea.Red(find(apnea.Time==MinsR.Time(i)):find(apnea.Time==MinsR.Time(i+1))));
    ir_dc(i) = mean(apnea.IR(find(apnea.Time==MinsIr.Time(i)):find(apnea.Time==MinsIr.Time(i+1))));

    red_ac(i) = ((MaxsR.Value(i)-MinsR.Value(i)+(MaxsR.Value(i+1)-MinsR.Value(i+1))))/2;
    ir_ac(i) = ((MaxsIr.Value(i)-MinsIr.Value(i)+(MaxsIr.Value(i+1)-MinsIr.Value(i+1))))/2;

    a=(red_ac(i)/red_dc(i));
    b=(ir_ac(i)/ir_dc(i));
    c=a/b;
    aR(i)=c;
    %j=j+1;
end

SPO2=zeros(1,length(aR));
for i=1:length(SPO2)
    %SPO2(i)=(110-25*aR(i));
    SPO2(i)=(110-25*aRFilt(i));
end

%% Dibujo saturación con apnea:
%Juntamos tiempo y saturación de una apnea
inicio_estudio=floor(apnea.Time(1));
fin_estudio=ceil(apnea.Time(end));
referencia_apnea=referencia(inicio_estudio:fin_estudio-1,1);
segundos_apnea = inicio_estudio:1:fin_estudio-1;
segundos_apnea=segundos_apnea';
r_apnea=addvars(referencia_apnea, segundos_apnea);

figure
subplot(2,1,1)
hold on
plot(MaxsR.Time(1:end-1), red_dc(:),'DisplayName', 'Red DC');
plot(MaxsR.Time(1:end-1), red_ac,'DisplayName', 'Red AC');
plot(MaxsR.Time(1:end-1), ir_dc,'DisplayName', 'Infra DC');
plot(MaxsR.Time(1:end-1), ir_ac,'DisplayName', 'Infra AC');
title('Grafico apnea con respiración media 2')
legend

subplot(2,1,2)
hold on
plot(r_apnea.segundos_apnea, r_apnea.ans, 'DisplayName', 'Sat Referencia');
plot(MaxsR.Time(1:end-1),SPO2-60,'DisplayName', 'SPO2');
%plot(MaxsR.Time(1:end-1),aR+91,'DisplayName', 'R+91');
plot(MaxsR.Time(1:end-1),aRFilt+90,'DisplayName', 'R apnea low-15');
legend

%% Intervalo de sueño normal
%Obtencion de los instantes temporales donde se producen máximos y minimos en Red e IR
%Obtenemos un intervalo representativo de apnea
Red_n = FiltroPBajoRoja(176*60*256:196*60*256); %176'-196'
IR_n = FiltroPBajoIR(176*60*256:196*60*256);
Time_n=remuestreado.Time(find(remuestreado.Time==(176*60)):find(remuestreado.Time==(196*60)));

normal=table(Time_n, Red_n, IR_n, 'VariableNames', {'Time', 'Red', 'IR'});
clear Red_n IR_n Time_n

%Claculamos el intervalo donde estan los indices con maximos y minimos de
%los latidos
%IR
maxIR = (n_iA>(176*60) & n_iA<(196*60));
maxIR = find(maxIR==1);
n_iA=n_iA(maxIR(1):maxIR(end));

minIR = (n_iB>(176*60) & n_iB<(196*60));
minIR = find(minIR==1);
n_iB=n_iB(minIR(1):minIR(end));

%Red
maxRed = (n_nA>(176*60) & n_nA<(196*60));
maxRed = find(maxRed==1);
n_nA=n_nA(maxRed(1):maxRed(end));

minRed = (n_nB>(176*60) & n_nB<(196*60));
minRed = find(minRed==1);
n_nB=n_nB(minRed(1):minRed(end));

clear minRed minIR maxIR maxRed

latidos = min(length([n_nA, n_nB, n_iA, n_iB]));
Value=nan(latidos,1);

MaxsIr=table(n_iA, Value, 'VariableNames', {'Time', 'Value'});
MinsIr=table(n_iB, Value, 'VariableNames', {'Time', 'Value'});
MaxsR=table(n_nA, Value, 'VariableNames', {'Time', 'Value'});
MinsR=table(n_nB, Value, 'VariableNames', {'Time', 'Value'});
clear Value  %n_iA n_iB n_nA n_nB

%% Asignacion correcta del valor de la onda en los valores de tiempo dados 
for i=1:latidos
    %Maximos de IR
    if isempty(find(normal.Time==MaxsIr.Time(i), 1)) %aproximar al valor correcto más cercano
        %Optimizacion para la busqueda, sumamos los elementos true del
        %vector para hayar el indice que buscamos
        ant=sum(normal.Time<MaxsIr.Time(i));
        sig=ant+1;
        if normal.IR(ant)>normal.IR(sig)
            MaxsIr.Time(i)=normal.Time(ant);
            MaxsIr.Value(i)=normal.IR(ant);
        else
            MaxsIr.Time(i)=normal.Time(sig);
            MaxsIr.Value(i)=normal.IR(sig);
        end
    else%tenemos valor para el tiempo dado
        MaxsIr.Value(i)=normal.IR(find(normal.Time==MaxsIr.Time(i)));
    end

    %Mínimos de IR
    if isempty(find(normal.Time==MinsIr.Time(i), 1)) 
        ant=sum(normal.Time<MinsIr.Time(i));
        sig=ant+1;
        if normal.IR(ant)<normal.IR(sig)
            MinsIr.Time(i)=normal.Time(ant);
            MinsIr.Value(i)=normal.IR(ant);
        else
            MinsIr.Time(i)=normal.Time(sig);
            MinsIr.Value(i)=normal.IR(sig);
        end
    else
        MinsIr.Value(i)=normal.IR(find(normal.Time==MinsIr.Time(i)));
    end
    
    %Maximos de Red
    if isempty(find(normal.Time==MaxsR.Time(i), 1))
        ant=sum(normal.Time<MaxsR.Time(i));
        sig=ant+1;
        if normal.Red(ant)>normal.Red(sig)
            MaxsR.Time(i)=normal.Time(ant);
            MaxsR.Value(i)=normal.Red(ant);
        else
            MaxsR.Time(i)=normal.Time(sig);
            MaxsR.Value(i)=normal.Red(sig);
        end
    else
        MaxsR.Value(i)=normal.Red(find(normal.Time==MaxsR.Time(i)));
    end

    %Mínimos de Red
    if isempty(find(normal.Time==MinsR.Time(i), 1)) 
        ant=sum(normal.Time<MinsR.Time(i));
        sig=ant+1;
        if normal.Red(ant)<normal.Red(sig)
            MinsR.Time(i)=normal.Time(ant);
            MinsR.Value(i)=normal.Red(ant);
        else
            MinsR.Time(i)=normal.Time(sig);
            MinsR.Value(i)=normal.Red(sig);
        end
    else
        MinsR.Value(i)=normal.Red(find(normal.Time==MinsR.Time(i)));
    end
end

%% Calculo de R y SPO2 normales
%Declaramos por optimización
red_dc=zeros(latidos-1,1);
ir_dc=zeros(latidos-1,1);
red_ac=zeros(latidos-1,1);
ir_ac=zeros(latidos-1,1);
nR=zeros(latidos-1,1);

for i=1:latidos-1
    %Cogemos el minimo de un latido y el siguiente y sacamos la media del
    %valor de la señal entre esos puntos
    red_dc(i) = mean(normal.Red(find(normal.Time==MinsR.Time(i)):find(normal.Time==MinsR.Time(i+1))));
    ir_dc(i) = mean(normal.IR(find(normal.Time==MinsIr.Time(i)):find(normal.Time==MinsIr.Time(i+1))));

    red_ac(i) = ((MaxsR.Value(i)-MinsR.Value(i)+(MaxsR.Value(i+1)-MinsR.Value(i+1))))/2;
    ir_ac(i) = ((MaxsIr.Value(i)-MinsIr.Value(i)+(MaxsIr.Value(i+1)-MinsIr.Value(i+1))))/2;

    a=(red_ac(i)/red_dc(i));
    b=(ir_ac(i)/ir_dc(i));
    c=a/b;
    nR(i)=c;
    %j=j+1;
end

SPO2=zeros(1,length(nR));
for i=1:length(SPO2)
    %SPO2(i)=(110-25*nR(i));
    SPO2(i)=(110-25*nRFilt(i));
end

%% Dibujo saturación normal:
%Juntamos tiempo y saturación de una apnea
inicio_estudio=floor(normal.Time(1));
fin_estudio=ceil(normal.Time(end));
referencia_normal=referencia(inicio_estudio:fin_estudio-1,1);
segundos_normal = inicio_estudio:1:fin_estudio-1;
segundos_normal=segundos_normal';
r_norm=addvars(referencia_normal, segundos_normal);

figure
subplot(2,1,1)
hold on
plot(MaxsR.Time(1:end-1), red_dc(:),'DisplayName', 'Red DC');
plot(MaxsR.Time(1:end-1), red_ac,'DisplayName', 'Red AC');
plot(MaxsR.Time(1:end-1), ir_dc,'DisplayName', 'Infra DC');
plot(MaxsR.Time(1:end-1), ir_ac,'DisplayName', 'Infra AC');
title('Grafico normal con respiración media 2')
legend

subplot(2,1,2)
hold on
plot(r_norm.segundos_normal, r_norm.ans, 'DisplayName', 'Sat Referencia');
plot(MaxsR.Time(1:end-1),SPO2-53,'DisplayName', 'SPO2');
%plot(MaxsR.Time(1:end-1),nR+91,'DisplayName', 'R+91');
plot(MaxsR.Time(1:end-1),nRFilt+91,'DisplayName', 'R normal low-15');
legend

%% Limpiamos variables
clear a b c ant fc sig i latidos 
%% Plot FFT red normal 
%Separar el codigo y meterlo cada cacho en su sección (El de la apnea a la
%apnea y el normal al normal, para que no haya fallos a la hora de
%ejecutar, que si no el orden de ejecucion del programa es muy rarete
%Basado en doc fft
fourierApn=fft(aR);
fourierNor=fft(nR);
aL=length(aR);
nL=length(nR);

A2=abs(fourierApn/aL);
N2=abs(fourierNor/nL);
A1=A2(1:aL/2+1);
A1(2:end-1)= 2*A1(2:end-1);
N1=N2(1:nL/2+1);
N1(2:end-1)= 2*N1(2:end-1);

fApn=fo/aL*(0:(aL/2));
fNorm=fo/nL*(0:(nL/2));

figure;
subplot(2,1,1)
hold on
%Nos cargamos la componente constante en f=0
plot(fApn(2:end),A1(2:end),'DisplayName', ' R durante apnea');
xlabel("f (Hz)")
ylabel("|fft(R apnea)|")
legend

subplot(2,1,2)
%Nos cargamos la componente constante en f=0
plot(fNorm(2:end),N1(2:end),'DisplayName', 'R tramo normal');
ylim([0,0.08])
xlabel("f (Hz)")
ylabel("|fft(R normal)|")
legend

%En vistas de las graficas nos interesa establecer el corte en 15Hz, ya que
%es el punto suficientemente alejado como para no perder precisión, pero
%suficientemente cercano para limpiar bien la señal.
fc = 20/(fo/2); 
[z,p,k]         =   cheby2( ord , 20 , fc , 'low' ); 
[sos,filterGain] = zp2sos(z,p,k);
aRFilt = nanfiltfilt( sos , filterGain , -aR, Setup );
nRFilt = nanfiltfilt( sos , filterGain , -nR, Setup );









%% De aqui para abajo borrador
%Búsqueda correlación de la señal SPO2
%Calculos Apnea
%Obtenemos un intervalo representativo de apnea
intervalo = (remuestreado.Time>(40*60)) & (remuestreado.Time<(60*60));%40'-60'
intervalo = find(intervalo==1);
apnea = remuestreado(intervalo(1):intervalo(end)+1, 1:4); %+1 --> %256==0

%Juntamos tiempo y saturación de una apnea
inicio_estudio=floor(normal.Time(1));
fin_estudio=ceil(normal.Time(end));
referencia_apnea=referencia(inicio_estudio:fin_estudio-1,1);
segundos_apnea = inicio_estudio:1:fin_estudio-1;
segundos_apnea=segundos_apnea';
r_apnea=addvars(referencia_apnea, segundos_apnea);

%Intentamos aproximar el calculo del SPO2 para la apnea
j=1;
[PVInfr, picosInfra] = findpeaks(apnea.IC, apnea.Time,'MinPeakDistance',0.5);
[PVRed, picosRed] = findpeaks(apnea.RC, apnea.Time,'MinPeakDistance',0.5);

%Media con 2
for i=1:length(picosRed)-1
    %Buscamos el pico y hacemos una media del valor de la señal hasta el
    %siguiente
    red_dc(j) = mean(apnea.RC(find(apnea.Time==picosRed(i)):find(apnea.Time==picosRed(i+1))));
    ir_dc(j) = mean(apnea.IC(find(apnea.Time==picosInfra(i)):find(apnea.Time==picosInfra(i+1))));

    red_ac(j) = max(apnea.RC(find(apnea.Time==picosRed(i)):find(apnea.Time==picosRed(i+1))))-min(apnea.RC(find(apnea.Time==picosRed(i)):find(apnea.Time==picosRed(i+1))));
    ir_ac(j) = max(apnea.IC(find(apnea.Time==picosInfra(i)):find(apnea.Time==picosInfra(i+1))))-min(apnea.IC(find(apnea.Time==picosInfra(i)):find(apnea.Time==picosInfra(i+1))));

    a=(red_ac(j)/red_dc(j));
    b=(ir_ac(j)/ir_dc(j));
    c=a/b;
    R(j)=c;
    j=j+1;
end

SPO2=zeros(1,length(R));
for i=1:length(SPO2)
    SPO2(i)=(110-25*R(i));
end
%% Ploteo saturacion apnea
figure
subplot(2,1,1)
hold on
plot(picosRed(1:end-1), red_dc(1:end),'DisplayName', 'Red DC');
plot(picosRed(1:end-1), red_ac(1:end),'DisplayName', 'Red AC');
plot(picosInfra(1:end-2), ir_dc(1:end),'DisplayName', 'Infra DC');
plot(picosInfra(1:end-2), ir_ac(1:end),'DisplayName', 'Infra AC');
title("Grafico apnea unfiltered media 2")
legend

subplot(2,1,2)
hold on
plot(r_apnea.segundos_apnea, r_apnea.ans, 'DisplayName', 'Sat Referencia');
plot(picosRed(1:end-1),R+90,'DisplayName', 'R+90');
legend

%linkaxes([picosRed(1:end-1), r_apnea.segundos_apnea], 'x')

%% Plot FFT red apnea
%Basado en doc fft
fourierADC=fft(red_dc);
fourierAR=fft(red_ac);
L=length(red_ac);

P2=abs(fourierADC/L);
K2=abs(fourierAR/L);
P1=P2(1:L/2+1);
P1(2:end-1)= 2*P1(2:end-1);
K1=K2(1:L/2+1);
K1(2:end-1)= 2*K1(2:end-1);

frecuenciaApnea=fo/L*(0:(L/2));

figure;
subplot(2,1,1)
hold on
%Nos cargamos la componente constante en f=0
plot(frecuenciaApnea(2:end),P1(2:end),'DisplayName', ' Apnea Red DC');
xlabel("f (Hz)")
ylim([0,8])
ylabel("|fft(Red DC)|")
legend

%Nos cargamos la componente constante en f=0
subplot(2,1,2)
plot(frecuenciaApnea(2:end),K1(2:end),'DisplayName', ' Apnea R'); 
%ylim([0,0.08])
xlabel("f (Hz)")
ylabel("|fft(R)|")
legend

%% Calculos de sueño normal

%Obtenemos un intervalo representativo de sueño normal
intervalo = (remuestreado.Time>(176*60)) & (remuestreado.Time<(196*60));%176'-196'
intervalo = find(intervalo==1);
normal = remuestreado(intervalo(1):intervalo(end)+1, 1:4);%+1 --> %256==0

%Juntamos tiempo y saturación de sueño normal
inicio_estudio=floor(normal.Time(1));
fin_estudio=ceil(normal.Time(end));
referencia_normal=referencia(inicio_estudio:fin_estudio-1,1);
segundos_normal = inicio_estudio:1:fin_estudio-1;
segundos_normal=segundos_normal';
r_normal=addvars(referencia_normal, segundos_normal);

%Intentamos aproximar el calculo del SPO2 para tramo normal
k=1;
[nInfr, nPicosInfra] = findpeaks(normal.IC, normal.Time, 'MinPeakDistance', 0.5);
[nRed, nPicosRed] = findpeaks(normal.RC, normal.Time,'MinPeakDistance',0.5);

%media con 2
for i=1:length(nPicosRed)-1
    %Buscamos el pico y hacemos una media del valor de la señal hasta el
    %siguiente
    n_red_dc(k) = mean(normal.RC(find(normal.Time==nPicosRed(i)):find(normal.Time==nPicosRed(i+1))));
    n_ir_dc(k) = mean(normal.IC(find(normal.Time==nPicosInfra(i)):find(normal.Time==nPicosInfra(i+1))));

    n_red_ac(k) = max(normal.RC(find(normal.Time==nPicosRed(i)):find(normal.Time==nPicosRed(i+1))))-min(normal.RC(find(normal.Time==nPicosRed(i)):find(normal.Time==nPicosRed(i+1))));
    n_ir_ac(k) = max(normal.IC(find(normal.Time==nPicosInfra(i)):find(normal.Time==nPicosInfra(i+1))))-min(normal.IC(find(normal.Time==nPicosInfra(i)):find(normal.Time==nPicosInfra(i+1))));

    n_a=(n_red_ac(k)/n_red_dc(k));
    n_b=(n_ir_ac(k)/n_ir_dc(k));
    n_c=n_a/n_b;
    nR(k)=n_c;
    k=k+1;
end

SPO2=zeros(1,length(nR));
for i=1:length(SPO2)
    SPO2(i)=(110-25*nR(i));
end

%% Plot Saturacion noche normal (media de 2, para 3 restar 2 en vez de 1)
figure
subplot(2,1,1)
hold on
plot(nPicosRed(1:end-1), n_red_dc,'DisplayName', 'Red DC');
plot(nPicosRed(1:end-1), n_red_ac,'DisplayName', 'Red AC');
plot(nPicosInfra(1:end-1), n_ir_dc,'DisplayName', 'Infra DC');
plot(nPicosInfra(1:end-1), n_ir_ac,'DisplayName', 'Infra AC');
title("Grafico tranquilo")
legend

subplot(2,1,2)
hold on
plot(r_normal.segundos_normal, r_normal.ans, 'DisplayName', 'Sat Referencia');
plot(nPicosRed(1:end-1),nR+90,'DisplayName', 'R+90');
legend

%% Plot FFT red normal
%Basado en doc fft
fourierNDC=fft(n_red_dc);
fourierNR=fft(n_red_ac);
nL=length(nR);

R2=abs(fourierNDC/nL);
M2=abs(fourierNR/nL);
R1=R2(1:nL/2+1);
R1(2:end-1)= 2*R1(2:end-1);
M1=M2(1:nL/2+1);
M1(2:end-1)= 2*M1(2:end-1);

frecuenciaNormal=fo/nL*(0:(nL/2));

figure;
subplot(2,1,1)
hold on
%Nos cargamos la componente constante en f=0
plot(frecuenciaNormal(2:end),R1(2:end),'DisplayName', ' Normal Red DC');
ylim([0,8])
xlabel("f (Hz)")
ylabel("|fft(Red DC)|")
legend

subplot(2,1,2)
%Nos cargamos la componente constante en f=0
plot(frecuenciaNormal(2:end),M1(2:end),'DisplayName', ' Normal R');
%ylim([0,0.08])
xlabel("f (Hz)")
ylabel("|fft(R)|")
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
%Setup.plotflag  =	false;
%isArtifact      =   energyArtifacts ( PPGIR , fs , Setup );
%isArtifact      =   energyArtifacts ( PPGRed , fs , Setup );


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

