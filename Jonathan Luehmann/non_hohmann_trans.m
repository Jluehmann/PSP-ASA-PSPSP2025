mu_Sun = 132712440017.99; 

AU = 149597870.69;

mu_tat = 320000;
mu_end = 1395700;

r_tat = 4.5; % [AU]
r_end = 3.9;
r_a = 4.8;
a = 4.35;

vp_tat = sqrt(mu_Sun/(r_tat*AU));
e_ti = 1-(r_end/a);
P_ti = (a*AU)*(1-e_ti.^2);
h_ti = sqrt(mu_Sun*P_ti);
v_ti = sqrt(mu_Sun * (2/(r_tat*AU) - 1/(a*AU)));
gam_i = acos(h_ti/(r_tat*AU*v_ti));
v_infi = sqrt(v_ti.^2*vp_tat-2*v_ti*vp_tat*cos(gam_i));

vp_end = sqrt(mu_Sun/(r_end*AU));
v_tf = sqrt(mu_Sun * (2/(r_end*AU) - 1/(a*AU)));
gam_f = acos(h_ti/(r_end*AU*v_tf));
v_inff = sqrt(v_tf.^2*vp_end-2*v_tf*vp_end*cos(gam_f));

e_f = 1 + (100+24000)*(v_inff.^2/mu_end);

ta = abs(acos(P_ti/(r_tat*AU*e_ti)-(1/e_ti)));

figure
orbitplot2D(a*AU, e_ti, [0:0.01:2], 0, "Transfer", r_scale=AU); hold on
orbitplot2D(r_tat*AU, 0, [0:0.01:2*pi], 0, "Tatoine", r_scale=AU); hold on
orbitplot2D(r_end*AU, 0, [0:0.01:2*pi], 0, "Endor", r_scale=AU);
axis equal 
hold off