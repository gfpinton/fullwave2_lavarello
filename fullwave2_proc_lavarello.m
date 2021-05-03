%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2021-04-10
% LAST MODIFIED: 2021-04-10
% Launch Fullwave 2 code, easy matlab wrapper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load fullwave_attenuation_lesion.mat 

taxis=(1:size(tx1.pxducer,1))*tx1.dT;
xaxis=(1:size(tx1.pxducer,2))*tx1.dX;
imagesc(xaxis,taxis,tx1.pxducer)
cbar=colorbar, title(cbar,'Pressure (Pa)'), xlabel('distance (m)'), ylabel('time (s)')
saveFig(gcf,'figures/pxducer1')

imagesc(xaxis,taxis,tx2.pxducer)
cbar=colorbar, title(cbar,'Pressure (Pa)'), xlabel('distance (m)'), ylabel('time (s)')
saveFig(gcf,'figures/pxducer2')

plot(taxis,tx1.pxducer(:,round(end/2))), hold on
plot(taxis,tx2.pxducer(:,round(end/2))), hold off
grid on, label('Pressure (Pa)'), xlabel('time (s)')
legend('Case 1','Case 2')
saveFig(gcf,'figures/paxial')

plot(taxis,tx1.pxducer(:,round(end/2))-tx2.pxducer(:,round(end/2)))
grid on, label('Pressure (Pa)'), xlabel('time (s)')
saveFig(gcf,'figures/paxial_diff')


