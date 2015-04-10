m = 0:0.01:30;

% noise
mu1 = 10;
s1 = 1;
L1 = (1/sqrt(2*pi*s1^2)).*exp(-(m-mu1).^2./(2*s1^2));
P1 = 0.5;


% signal
mu2 = 20;
s2 = 3;
L2 = (1/sqrt(2*pi*s2^2)).*exp(-(m-mu2).^2./(2*s2^2));
P2 = 1-P1;


L11 = 1;  % correct reject
L12 = 2;  % false alarm
L21 = 4; % miss
L22 = 1;  % hit


EL1 = L1.*P1.*L11 + L2.*P2.*L21; % Estimated loss for saying "no"
EL2 = L1.*P1.*L12 + L2.*P2.*L22; % Estimated loss for saying "yes"




% plot
figure(1);
clf;
plot(m,L1,'b-');
hold on;

plot(m,L2,'r-');
plot(m,EL1,'c--');
plot(m,EL2,'m-');

q = find(EL1>EL2);

if size(q,2)~=0
plot([m(q(1)) m(q(1))],[0 2],'g-','linewidth',1);
end
plot(m,zeros(size(m)),'k:');
xlabel('sensory signal m','fontsize',20);
ylabel('probdensity/loss-gain','fontsize',20);
%axis([xlim 0 1]);
