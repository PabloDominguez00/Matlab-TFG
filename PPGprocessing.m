
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
        t = ((huecos(i-1)+1)/fo):1/fo:(huecos(i-1)+length(GC))/fo;

        resampled{end+1}=table(t',GC,RC,IC, 'VariableNames', {'Time', 'GC', 'RC', 'IC'});
        
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
end

remuestreado =  resampled{1}; %Los intervalos a analizar se encuentran en 1
%segundos= 1:1:length(referencia.ans); %Referencia sin decimales

%Borramos las variable que no vamos a usar mas adelante
clear g r ir fi ffin i GC RC IC t fs huecos 
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

%Instantes de tiempo de: D->Latido, A->Max, B->Min 
%(partiendo desde el instante de tiempo de medicion inicial)
%Es decir, que son isntantes desde que inicia la tabla válida
[ iD , iA , iB ]    =   pulseDelineation ( PPGIR , fo , Setup );
[ nD , nA , nB ]    =   pulseDelineation ( PPGRed , fo , Setup ); 

%creamos un eje temporal que se ajuste a nuestra señal
t = 0:1/fo:(length(PPGRed)-1)/fo;

a_iA=iA; a_iB=iB; %para parte de apnea
a_nA=nA; a_nB=nB;

n_iA=iA; n_iB=iB; %para parte noche normal
n_nA=nA; n_nB=nB;  

%Con tiempos arreglados:
irA=iA+remuestreado.Time(1); irB=iB+remuestreado.Time(1);
rA=nA+remuestreado.Time(1); rB=nB+remuestreado.Time(1);

clear k p z sos filterGain
%% Creamos nuevo filtro para mantener la respiración
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
plot(iB(~isnan(iB)), FiltroPBajoIR(1+round(iB(~isnan(iB))*fo)), 'b*','LineWidth',1);
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
plot(nA(~isnan(nA)), FiltroPBajoRoja(1+round(nA(~isnan(nA))*fo)), 'b*','LineWidth',1);
plot(nB(~isnan(nB)), FiltroPBajoRoja(1+round(nB(~isnan(nB))*fo)), 'b*','LineWidth',1);
title('PPGRED');

clear k p z sos filterGain


%% Toda la noche (creacion de las tablas)
noche=table(remuestreado.Time, FiltroPBajoRoja, FiltroPBajoIR,'VariableNames', {'Time', 'Red', 'IR'});

MaxsIr=table(irA, nan(length(iA),1), 'VariableNames', {'Time', 'Value'});
MinsIr=table(irB, nan(length(iB),1), 'VariableNames', {'Time', 'Value'});
MaxsR=table(rA, nan(length(nA),1), 'VariableNames', {'Time', 'Value'});
MinsR=table(rB, nan(length(nB),1), 'VariableNames', {'Time', 'Value'});

find(MaxsIr.Time<=noche.Time(1), 1)
%clear nA nB iA iB

%% Toda la noche (asignacion de valores a las tablas de tiempos)
% Falta aunar en una sola tabla las muestras resampleadas
%Este apartado cuesta un cojon de ejecutar xdddd
tic
for i=1:length(MaxsIr.Value)
    %Maximos de IR
    idx = find(noche.Time==MaxsIr.Time(i), 1);
    if isempty(idx) %aproximar al valor correcto más cercano
        %Optimizacion para la busqueda, sumamos los elementos true del
        %vector para hayar el indice que buscamos
        ant=sum(noche.Time<MaxsIr.Time(i));
        sig=ant+1;
        if noche.IR(ant)>noche.IR(sig)%
            MaxsIr.Time(i)=noche.Time(ant);
            MaxsIr.Value(i)=noche.IR(ant);
        else
            MaxsIr.Time(i)=noche.Time(sig);
            MaxsIr.Value(i)=noche.IR(sig);
        end
    else%tenemos valor para el tiempo dado
        MaxsIr.Value(i)=noche.IR(idx);
    end

    %Mínimos de IR
    idx = find(noche.Time==MinsIr.Time(i), 1);
    if isempty(idx) 
        ant=sum(noche.Time<MinsIr.Time(i));
        sig=ant+1;
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
toc
for i=1:length(MaxsR.Value)
    %Maximos de Red
    idx = find(noche.Time==MaxsR.Time(i), 1);
    if isempty(idx)
        ant=sum(noche.Time<MaxsR.Time(i));
        sig=ant+1;
        if noche.Red(ant)>noche.Red(sig)
            MaxsR.Time(i)=noche.Time(ant);
            MaxsR.Value(i)=noche.Red(ant);
        else
            MaxsR.Time(i)=noche.Time(sig);
            MaxsR.Value(i)=noche.Red(sig);
        end
    else
        MaxsR.Value(i)=noche.Red(idx);
    end

    %Mínimos de Red
    idx = find(noche.Time==MinsR.Time(i), 1);
    if isempty(idx) 
        ant=sum(noche.Time<MinsR.Time(i));
        sig=ant+1;
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


%% Toda la noche (calculo de R y SPO2)

%Revisar como podriasmos solucionar el tema de los latido totales:
latidos=min(height(MaxsR),height(MaxsIr));
%Declaramos por optimización
red_dc=zeros(latidos-1,1);
ir_dc=zeros(latidos-1,1);
red_ac=zeros(latidos-1,1);
ir_ac=zeros(latidos-1,1);
aR=zeros(latidos-1,1);

for i=1:latidos-1
    %Cogemos el minimo de un latido y el siguiente y sacamos la media del
    %valor de la señal entre esos puntos
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
inicio_estudio=floor(noche.Time(1));
fin_estudio=ceil(noche.Time(end));
referencia_noche=referencia(inicio_estudio:fin_estudio-1,1);
segundos = inicio_estudio:1:fin_estudio-1;
segundos=segundos';
ref_noche=addvars(referencia_noche, segundos);

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
plot(ref_noche.segundos, ref_noche.ans, 'DisplayName', 'Sat Referencia');
yyaxis right
%plot(MaxsIr.Time(1:end-1),SPO2,'DisplayName', 'SPO2');
%plot(prueba1(1:end-1),SPO2,'DisplayName', 'SPO2');
%plot(MaxsR.Time(1:end-1),aR+91,'DisplayName', 'R+91');
% ylim([111 112])
plot(MaxsIr.Time(1:end-1),-aRFilt,'DisplayName', 'R apnea low-15');
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

