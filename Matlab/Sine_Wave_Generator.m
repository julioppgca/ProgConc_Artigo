% Create sine waves for test pourpouses.
% Write to a file in C language matrix with size variable of uint16_t.
% Values up to 12bits.
% Julio Santos - julio.ppgca@gmail.com
% Nov. 2016
% Ver. 1

format('shortEng');
t=linspace(0,1/60,2^8); % time vector, 256 values from 0 to 1/60.
% Sine Wave, pure 60Hz.
y1 = 0.8*sin(2*pi*60*t);
sinal=ceil(y1*2048)+2047;  % should be integer
plot(t,y1) % Just to see how beauty is the Matlab.
title('Generated Sine waves');
xlabel('Time (ms)');
ylabel('Voltage (V)');
hold on
grid

% Write to a file, then just copy and paste ;P
fileID = fopen('SineWave.h','w');
fprintf(fileID,'#ifdef SINE_WAVE_CLEAN \n');
fprintf(fileID,'// t=linspace(0,1/60,2^8);  \n');
fprintf(fileID,'// y1 = 0.8*sin(2*pi*60*t);  \n');
fprintf(fileID,'// sinal=ceil(y1*2048)+2047;  \n');
fprintf(fileID,'// RMS: %.3f Volts  \n', rms(y1));
fprintf(fileID,'uint16_t  SineWave[2][SAMPLES_NUM] = \n{');
fprintf(fileID,'\n\t{\n\t');
j=1;
for i=1:8
  fprintf(fileID,'%d,\t',sinal(j:j+15));
  fprintf(fileID,'\n\t');
  j=j+16;
end
fprintf(fileID,'}, \n\t{\n\t');
j=129;
for i=1:8
  fprintf(fileID,'%d,\t',sinal(j:j+15));
  fprintf(fileID,'\n\t');
  j=j+16;
end
fprintf(fileID,'}');
fprintf(fileID,'\n};\n');
fprintf(fileID,'#endif \n');

% Sine wave with some harmonics, 60Hz, 120Hz and 180Hz
y2 = 0.5*sin(2*pi*60*t)+0.3*sin(2*pi*120*t)+0.05*sin(2*pi*180*t);
sinal=ceil(y2*2048)+2047; 
plot(t,y2)

fprintf(fileID,'\n#ifdef SINE_WAVE_N1 \n');
fprintf(fileID,'// t=linspace(0,1/60,2^8);  \n');
fprintf(fileID,'// y2 = 0.5*sin(2*pi*60*t)+0.3*sin(2*pi*120*t)+0.05*sin(2*pi*180*t);  \n');
fprintf(fileID,'// sinal=ceil(y2*2048)+2047;  \n');
fprintf(fileID,'// RMS: %.3f Volts  \n', rms(y2));
fprintf(fileID,'uint16_t  SineWave[2][SAMPLES_NUM] = \n{');
fprintf(fileID,'\n\t{\n\t');
j=1;
for i=1:8
  fprintf(fileID,'%d,\t',(sinal(j:j+15)));
  fprintf(fileID,'\n\t');
  j=j+16;
end
fprintf(fileID,'}, \n\t{\n\t');
j=129;
for i=1:8
  fprintf(fileID,'%d,\t',sinal(j:j+15));
  fprintf(fileID,'\n\t');
  j=j+16;
end
fprintf(fileID,'}');
fprintf(fileID,'\n};\n');
fprintf(fileID,'#endif\n');
% Sine wave with some harmonics, 60Hzh and 720Hz
y3 = 0.5*sin(2*pi*60*t)+0.1*sin(2*pi*720*t);
sinal=ceil(y3*2048)+2047; 
plot(t,y3)
legend((['S1: ' num2str(rms(y1)) 'Vrms - 60Hz.' ]), ...
       (['S2: ' num2str(rms(y3)) 'Vrms - 60, 120 and 180Hz.' ]), ...
       (['S3: ' num2str(rms(y2)) 'Vrms - 60 and 720Hz.' ]));
hold off


fprintf(fileID,'\n#ifdef SINE_WAVE_N2 \n');
fprintf(fileID,'// t=linspace(0,1/60,2^8);  \n');
fprintf(fileID,'// y3 = 0.5*sin(2*pi*60*t)+0.1*sin(2*pi*720*t);  \n');
fprintf(fileID,'// sinal=ceil(y3*2048)+2047;  \n');
fprintf(fileID,'// RMS: %.3f Volts  \n', rms(y3));
fprintf(fileID,'uint16_t  SineWave[2][SAMPLES_NUM] = \n{');
fprintf(fileID,'\n\t{\n\t');
j=1;
for i=1:8
  fprintf(fileID,'%d,\t',(sinal(j:j+15)));
  fprintf(fileID,'\n\t');
  j=j+16;
end
fprintf(fileID,'}, \n\t{\n\t');
j=129;
for i=1:8
  fprintf(fileID,'%d,\t',sinal(j:j+15));
  fprintf(fileID,'\n\t');
  j=j+16;
end
fprintf(fileID,'}');
fprintf(fileID,'\n};\n');
fprintf(fileID,'#endif \n');
fclose(fileID);


