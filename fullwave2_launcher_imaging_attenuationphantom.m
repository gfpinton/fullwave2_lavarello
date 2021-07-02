%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2018-06-21
% LAST MODIFIED: 2021-06-24
% Launch Fullwave 2 code, easy matlab wrapper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /celerina/gfp/mfs/dumbmat/
%nohup matlab -nodisplay -nosplash < fullwave2_launcher_imaging_homog.m > output.txt &
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;         % speed of sound (m/s)
rho0=1000;
omega0 = 2*pi*5.5e6; % center radian frequency of transmitted wave
wX = 4e-2;         % width of simulation field (m)
wY = 6e-2;         % depth of simulation field (m)
duration = wY*2.5/c0;  % duration of simulation (s)
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
rhosr=0.05;
basedir='/kulm/scratch/imaging_attenuationphantom/'
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
outmap = zeros(nX,nY);  outmap(:,9) = ones(nX,1);
outcoords = mapToCoords(outmap);
%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nlines=128;
beamwidth=round(ppw/2);
nXextend=nX+(nlines-1)*beamwidth;
cmap = ones(nXextend,nY)*1540;   % speed of sound map (m/s)
rhomap = ones(nXextend,nY)*1000; % density map (kg/m^3)
Amap = ones(nXextend,nY)*0.5;    % attenuation map (dB/MHz/cm)
betamap = ones(nXextend,nY)*0.0;    % nonlinearity map 
imagesc(cmap'), axis equal, axis tight
rhos=rand(nXextend,nY); rhos(find(rhos>0.05))=0; rhos=rhos/0.05;
rhos(:,1:10)=0;
rhomap=rhomap-rhos*rho0*rhosr;
imagesc((1:nX)*dX,(1:nY)*dY,rhomap'), axis equal, axis tight
xlabel('m'), ylabel('m'), cbar=colorbar; title(cbar,'kg/m^3')
%saveFig(gcf,'rhomap')
imagesc((1:nX)*dX,(1:nY)*dY,Amap'), axis equal, axis tight
xlabel('m'), ylabel('m'), cbar=colorbar; title(cbar,'dB/MHz/cm')
%saveFig(gcf,'Amap1') 

lld=0;

for n=1:nlines
  idxvec=(1:nX)+(n-1)*beamwidth;

  c=cmap(idxvec,:);
  rho=rhomap(idxvec,:);
  A=Amap(idxvec,:);
  beta=betamap(idxvec,:);

%%% Launch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  outdir=[ basedir '/txrx' num2str(n) '/'];
  eval(['!mkdir -p ' outdir]);
  eval(['!cp fullwave2_try6_nln_relaxing ' outdir]);
  cwd=pwd; addpath(cwd); cd(outdir)
  
  prep_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,c,rho,A,beta,incoords,outcoords,icmat)

  !nohup ./fullwave2_try6_nln_relaxing > output.txt &
  cd(cwd)

   [tmp lld]=system('uptime | awk  ''{print $11}'''); lld=str2num(lld)
   while(lld>30)
    [tmp lld]=system('uptime | awk  ''{print $11}'''); lld=str2num(lld)
    pause(60)
  end
  
end
%%% GENERATE IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncoordsout=size(outcoords,1);
pause(60)
nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
while(nRun<nT-10)
  pause(60)
  nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
end
  

deps = 10e-3:lambda/20:nY*dY;
lats = 0;
fnumber=1;

xducercoords = outcoords;
bm=zeros(length(lats),length(deps),nlines);
idps=cell(length(lats),length(deps));
for ii=1:length(lats)
  lat=lats(ii);
  for jj=1:length(deps)
    dep=deps(jj);
    fceni=round([lat/dY+mean(xducercoords(:,1)) dep/dY ]);
    idx=find(abs(xducercoords(:,1)-fceni(1))<=fceni(2)/fnumber);
    dd=focusProfile(fceni,xducercoords(idx,:),dT/dY*c0);
    idt=round(2*dep/double(c0)/(dT));
    idp=double((nT*(idx-1))+double(idt)+dd);
    idps{ii,jj}=idp;
  end
end

ncoordsout=size(outcoords,1);
idc=find(outcoords(:,1)>-1);
outdir=[ basedir '/txrx' num2str(n) '/'];
pxducer=readpx(outdir,nT,ncoordsout,idc);

for n=1:nlines
  outdir=[ basedir '/txrx' num2str(n) '/'];
  pxducer=readpx(outdir,nT,ncoordsout,idc);
  imagesc(powcompress(pxducer,1/4))
  if(n==1)
    px=pxducer(:,round(size(pxducer,2)/2));
    [val idt0]=max(abs(hilbert(px)))
  end
  for ii=1:length(lats)
    for jj=1:length(deps)
      bm(ii,jj,n)=sum(round(pxducer(idps{ii,jj}+idt0)));
    end
  end
end
%% PLOT THE BMODE IMAGE %%
figure(1)
n=1:nlines; bws=((n-(nlines+1)/2)*beamwidth)*dY;
imagesc(bws*1e3,deps*1e3,dbzero(abs(hilbert(squeeze(bm)))),[-40 0])
colormap gray
xlabel('mm'), ylabel('mm')
axis equal, axis tight

bm_ref=squeeze(bm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOW WITH ATTENUATION% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[idx]=circleIdx([nXextend nY],[nXextend/2 fcen(2)],1.5e-2/dX/2);
Amap(idx)=1;
%rhomap(idx)=(rhomap(idx)-rho0)/10+rho0;

for n=1:nlines
  idxvec=(1:nX)+(n-1)*beamwidth;

  c=cmap(idxvec,:);
  rho=rhomap(idxvec,:);
  A=Amap(idxvec,:);
  beta=betamap(idxvec,:);

%%% Launch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  outdir=[ basedir '/txrxA' num2str(n) '/'];
  eval(['!mkdir -p ' outdir]);
  eval(['!cp fullwave2_try6_nln_relaxing ' outdir]);
  cwd=pwd; addpath(cwd); cd(outdir)
  
  prep_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,c,rho,A,beta,incoords,outcoords,icmat)

  !nohup ./fullwave2_try6_nln_relaxing > output.txt &
  cd(cwd)

  [tmp lld]=system('uptime | awk  ''{print $11}'''); lld=str2num(lld)
   while(lld>30)
    [tmp lld]=system('uptime | awk  ''{print $11}'''); lld=str2num(lld)
    pause(60)
  end
  
end

pause(60)
nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
while(nRun<nT-10)
  pause(60)
  nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
end
 

for n=1:nlines
  outdir=[ basedir '/txrxA' num2str(n) '/'];
  pxducer=readpx(outdir,nT,ncoordsout,idc);
  imagesc(powcompress(pxducer,1/4))
  if(n==1)
    px=pxducer(:,round(size(pxducer,2)/2));
    [val idt0]=max(abs(hilbert(px)))
  end
  for ii=1:length(lats)
    for jj=1:length(deps)
      bm(ii,jj,n)=sum(round(pxducer(idps{ii,jj}+idt0)));
    end
  end
end
%% PLOT THE BMODE IMAGE %%
figure(1)
n=1:nlines; bws=((n-(nlines+1)/2)*beamwidth)*dY;
imagesc(bws*1e3,deps*1e3,dbzero(abs(hilbert(squeeze(bm)))),[-40 0])
colormap gray
xlabel('mm'), ylabel('mm')
axis equal, axis tight


imagesc(bws*1e3,deps*1e3,squeeze(bm))
xlabel('mm'), ylabel('mm')
axis equal, axis tight

bm_A=squeeze(bm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOW WITH ATTENUATION + REFLECTIVITY REDUCTION %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[idx]=circleIdx([nXextend nY],[nXextend/2 fcen(2)],1.5e-2/dX/2);
Amap(idx)=1;
rhomap(idx)=(rhomap(idx)-rho0)/10+rho0;

for n=1:nlines
  idxvec=(1:nX)+(n-1)*beamwidth;

  c=cmap(idxvec,:);
  rho=rhomap(idxvec,:);
  A=Amap(idxvec,:);
  beta=betamap(idxvec,:);

%%% Launch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  outdir=[ basedir '/txrxAR' num2str(n) '/'];
  eval(['!mkdir -p ' outdir]);
  eval(['!cp fullwave2_try6_nln_relaxing ' outdir]);
  cwd=pwd; addpath(cwd); cd(outdir)
  
  prep_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,c,rho,A,beta,incoords,outcoords,icmat)

  !nohup ./fullwave2_try6_nln_relaxing > output.txt &
  cd(cwd)

  [tmp lld]=system('uptime | awk  ''{print $11}'''); lld=str2num(lld)
   while(lld>30)
    [tmp lld]=system('uptime | awk  ''{print $11}'''); lld=str2num(lld)
    pause(60)
  end

end

pause(60)
nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
while(nRun<nT-10)
  pause(60)
  nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
end
 

for n=1:nlines
  outdir=[ basedir '/txrxAR' num2str(n) '/'];
  pxducer=readpx(outdir,nT,ncoordsout,idc);
  imagesc(powcompress(pxducer,1/4))
  if(n==1)
    px=pxducer(:,round(size(pxducer,2)/2));
    [val idt0]=max(abs(hilbert(px)))
  end
  for ii=1:length(lats)
    for jj=1:length(deps)
      bm(ii,jj,n)=sum(round(pxducer(idps{ii,jj}+idt0)));
    end
  end
end
%% PLOT THE BMODE IMAGE %%
figure(1)
n=1:nlines; bws=((n-(nlines+1)/2)*beamwidth)*dY;
imagesc(bws*1e3,deps*1e3,dbzero(abs(hilbert(squeeze(bm)))),[-40 0])
colormap gray
xlabel('mm'), ylabel('mm')
axis equal, axis tight


imagesc(bws*1e3,deps*1e3,squeeze(bm))
xlabel('mm'), ylabel('mm')
axis equal, axis tight

bm_AR=squeeze(bm);

save imaging_attenuationphantom_rf bm_ref bm_A bm_AR bws deps lats
