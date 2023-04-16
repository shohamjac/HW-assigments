t = 0.02;
K_pd = 5; %V
K_VLC = 3350; %Hz/V
k = K_pd * K_VLC;

ws = k/t;
zis = 1 / (4 * k * t);

G = tf([k],[t,1,0]);
H = tf([k],[t,1,k]);
%H = tf([ws],[1,2*sqrt(zis)*sqrt(ws),ws])

bode(G)
figure(2)
bode(H)


%h = sigmaplot(G, 'r');
%legend('Reference H0','location','southwest')
%setoptions(h,'YlimMode','manual','Ylim',{[-60 0]})
%H2 = G/(1+G);

w = logspace(2,5.1,100);
