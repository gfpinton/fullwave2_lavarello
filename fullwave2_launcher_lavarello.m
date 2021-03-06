%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2018-06-21
% LAST MODIFIED: 2021-03-01
% Launch Fullwave 2 code, easy matlab wrapper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /celerina/gfp/mfs/dumbmat/
%nohup matlab -nodisplay -nosplash < fullwave2_launcher_lavarello.m > output.txt &
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;         % speed of sound (m/s)
omega0 = 2*pi*5.5e6; % center radian frequency of transmitted wave
wX = 4e-2;         % width of simulation field (m)
wY = 6e-2;         % depth of simulation field (m)
duration = wY*2.3/c0;  % duration of simulation (s)
p0 = 1e5; % pressure in Pa
%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = 12;           % number of points per spatial wavelength
cfl = 0.4;         % Courant-Friedrichs-Levi condition
%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c0/omega0*2*pi;
nX = round(wX/lambda*ppw);  % number of lateral elements
nY = round(wY/lambda*ppw);  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl);
dX = c0/omega0*2*pi/ppw
dY= c0/omega0*2*pi/ppw
dT = dX/c0*cfl;

ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nX,nY); 
inmap(:,1:8)=ones(nX,8);
imagesc(inmap'), axis equal, axis tight
incoords = mapToCoords(inmap); % note zero indexing for compiled code
plot(incoords(:,1),incoords(:,2),'.')
%%% Generate initial conditions based on input coordinates %%%%%%
ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
fcen=[round(nX/2) round(nY/1.3)]; % center of focus
t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec)
icmat=repmat(icvec,size(incoords,1)/8,1);
for k=2:8
  t=t-dX/c0; icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
  icmat=[icmat' repmat(icvec,size(incoords,1)/8,1)']';
end
imagesc(icmat)

foc=round(3e-2/dX);
fcen=[round(nX/2) foc]; % center of focus

cmap = ones(nX,nY)*1540;   % speed of sound map (m/s)
cmap(fcen(1)-1:fcen(1)+1,fcen(2)-1:fcen(2)+1)=0.5*c0; % scatterer
imagesc(cmap'), axis equal, axis tight

t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
[icmat dd] = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/8,:),icvec,cfl);
for k=2:8
  t=t-dX/c0; icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
  icmat=[icmat' focusCoords(fcen(1),fcen(2),incoords((k-1)*size(incoords,1)/8+1:(k)*size(incoords,1)/8,:),icvec,cfl)']';
end
imagesc(icmat)


%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
modX=3; modY=3;
outcoords = coordsMatrix(nX,nY,modX,modY);
plot(outcoords(:,1),outcoords(:,2),'.')
%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = ones(nX,nY)*1540;   % speed of sound map (m/s)
rhomap = ones(nX,nY)*1000; % density map (kg/m^3)
Amap = ones(nX,nY)*0.5;    % attenuation map (dB/MHz/cm)
betamap = ones(nX,nY)*0.0;    % nonlinearity map 
imagesc(cmap'), axis equal, axis tight

rhos=rand(nX,nY); rhos(find(rhos>0.05))=0; rhos=rhos/0.05;
rhomap=rhomap-rhos*50;
imagesc((1:nX)*dX,(1:nY)*dY,rhomap'), axis equal, axis tight
xlabel('m'), ylabel('m'), cbar=colorbar; title(cbar,'kg/m^3')
saveFig(gcf,'rhomap')
imagesc((1:nX)*dX,(1:nY)*dY,Amap'), axis equal, axis tight
xlabel('m'), ylabel('m'), cbar=colorbar; title(cbar,'dB/MHz/cm')
saveFig(gcf,'Amap1') 

%%% write files  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prep_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
tic
eval('!./fullwave2_try6_nln_relaxing')
toc

pause(60)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
nRun=sizeOfFile('genout.dat')/4/ncoordsout

nRun=sizeOfFile('genout.dat')/4/ncoordsout
while(nRun<nT-3)
  pause(60)
  nRun=sizeOfFile('genout.dat')/4/ncoordsout
end

genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

tx1.pxducer=squeeze(p(:,ceil(9/modX),:));
tx1.dT=dT;
tx1.dX=dX*modX;
tx1.cmap=cmap;
tx1.Amap=Amap;
tx1.rhomap=rhomap;
tx1.betamap=betamap;

for i=1:20:size(p,1)
  imagesc(powcompress(squeeze(p(i,:,:)),1/3)), title(num2str(i)), drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now with attenuation %%
[idx]=circleIdx([nX nY],fcen,1.5e-2/dX/2);
Amap(idx)=1;
imagesc((1:nX)*dX,(1:nY)*dY,Amap'), axis equal, axis tight
xlabel('m'), ylabel('m'), cbar=colorbar; title(cbar,'dB/MHz/cm')
saveFig(gcf,'Amap2') 
%%% write files  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prep_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
tic
eval('!./fullwave2_try6_nln_relaxing & ')
toc

pause(60)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
nRun=sizeOfFile('genout.dat')/4/ncoordsout


nRun=sizeOfFile('genout.dat')/4/ncoordsout
while(nRun<nT-3)
  pause(60)
  nRun=sizeOfFile('genout.dat')/4/ncoordsout
end


genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

tx2.pxducer=squeeze(p(:,ceil(9/modX),:));
tx2.dT=dT;
tx2.dX=dX*modX;
tx2.cmap=cmap;
tx2.Amap=Amap;
tx2.rhomap=rhomap;
tx2.betamap=betamap;


for i=1:20:size(p,1)
  imagesc(powcompress(squeeze(p(i,:,:)),1/3)), title(num2str(i)), drawnow
end


save fullwave_attenuation_lesion.mat tx1 tx2

imagesc(powcompress(tx1.pxducer,1/4))
saveFig(gcf,'pxducer1') 
imagesc(powcompress(tx1.pxducer,1/4))
saveFig(gcf,'pxducer2') 






if(0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's add a surface %%%%%%%%%%%%%%%%%%%%%%%%%%
cmap(:,round(nY/2):end)=0.95*c0; % surface
imagesc(cmap')

prep_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
tic
eval('!./fullwave2_try6_nln_relaxing')
toc
ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

for i=1:20:size(p,1)
  imagesc(squeeze(p(i,:,:)), [-1 1]*p0), title(num2str(i)), drawnow
end
%% use powcompress to view small amplitude fields %%
for i=1:20:size(p,1)
  imagesc(powcompress(squeeze(p(i,:,:)),1/3)), title(num2str(i)), drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's add the scatterer back in and let's focus to the scatterer%%
foc=round(nY/1.3);
fcen=[round(nX/2) foc]; % center of focus

cmap = ones(nX,nY)*1540;   % speed of sound map (m/s)
cmap(fcen(1)-1:fcen(1)+1,fcen(2)-1:fcen(2)+1)=0.5*c0; % scatterer
imagesc(cmap'), axis equal, axis tight

t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
[icmat dd] = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/8,:),icvec,cfl);
for k=2:8
  t=t-dX/c0; icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
  icmat=[icmat' focusCoords(fcen(1),fcen(2),incoords((k-1)*size(incoords,1)/8+1:(k)*size(incoords,1)/8,:),icvec,cfl)']';
end
imagesc(icmat)

prep_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
tic
eval('!./fullwave2_try6_nln_relaxing')
toc
ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

for i=1:20:size(p,1)
  imagesc(squeeze(p(i,:,:)), [-1 1]*p0), title(num2str(i)), drawnow
end
%% use powcompress to view small amplitude fields %%
for i=1:20:size(p,1)
  imagesc(powcompress(squeeze(p(i,:,:)),1/3)), title(num2str(i)), drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's add another scatterer off axis and let's focus to the scatterer%%
foc=round(nY/1.3);
fcen=[round(nX/4) foc]; % center of focus

cmap = ones(nX,nY)*1540;   % speed of sound map (m/s)
cmap(fcen(1)-1:fcen(1)+1,fcen(2)-1:fcen(2)+1)=0.5*c0; % scatterer
imagesc(cmap'), axis equal, axis tight

t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
[icmat dd] = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/8,:),icvec,cfl);
for k=2:8
  t=t-dX/c0; icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
  icmat=[icmat' focusCoords(fcen(1),fcen(2),incoords((k-1)*size(incoords,1)/8+1:(k)*size(incoords,1)/8,:),icvec,cfl)']';
end
imagesc(icmat)

prep_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
tic
eval('!./fullwave2_try6_nln_relaxing')
toc
ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

for i=1:20:size(p,1)
  imagesc(squeeze(p(i,:,:)), [-1 1]*p0), title(num2str(i)), drawnow
end
%% use powcompress to view small amplitude fields %%
for i=1:20:size(p,1)
  imagesc(powcompress(squeeze(p(i,:,:)),1/3)), title(num2str(i)), drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's add many scatterers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foc=round(nY/1.3);
fcen=[round(nX/2) foc]; % center of focus

cmap = ones(nX,nY)*1540;   % speed of sound map (m/s)
rhomap = ones(nX,nY)*1000;
rhos=rand(nX,nY); rhos(find(rhos>0.05))=0; rhos=rhos/0.05;
rhomap=rhomap-rhos*50;
imagesc(rhomap'), axis equal, axis tight

t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
[icmat dd] = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/8,:),icvec,cfl);
for k=2:8
  t=t-dX/c0; icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
  icmat=[icmat' focusCoords(fcen(1),fcen(2),incoords((k-1)*size(incoords,1)/8+1:(k)*size(incoords,1)/8,:),icvec,cfl)']';
end
imagesc(icmat)

prep_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
tic
eval('!./fullwave2_try6_nln_relaxing')
toc
ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

%% use powcompress to view small amplitude fields %%
for i=1:20:size(p,1)
  imagesc(powcompress(squeeze(p(i,:,:)),1/3)), title(num2str(i)), drawnow
end


end
