t = 0.02;
K_pd = 5; %V
K_VLC = 3300; %Hz/V
k = K_pd * K_VLC;

G = tf([k],[t,1,k])