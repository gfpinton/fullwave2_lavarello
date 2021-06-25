function [pxducer] = readpx(outdir,nT,ncoordsout,varargin)

idc=1:ncoordsout;
optargin=size(varargin,2);
if(optargin>=1)
  idc=varargin{1};
end

nRun=sizeOfFile([outdir '/genout.dat'])/4/ncoordsout
if(nRun<nT)
  disp(['WARNING: nT<nRun, will zero pad from ' num2str(nRun) ' to ' num2str(nT)])
end
if(nRun>0 & nRun<=nT)
  pxducer=readGenoutSlice([outdir 'genout.dat'],0:nRun-1,ncoordsout,idc);
  pxducer(end:nT,:)=0;
end
if(nRun>nT)
  disp(['Too many time points in pxducer, truncating to ' num2str(nT)])
  pxducer=readGenoutSlice([outdir 'genout.dat'],0:nT-1,ncoordsout,idc);
end

