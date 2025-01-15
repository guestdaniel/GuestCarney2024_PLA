x = cosine_ramp(scale_dbspl(pure_tone(1e3, 0.0, 0.25, 100e3), 50.0), 0.01, 100e3);
ihc = sim_ihc_zbc2014(x, 1e3);
an1 = sim_an_zbc2025(ihc, 1e3, implnt=1);
an2 = sim_an_zbc2025(ihc, 1e3, implnt=2);
plot(an1); hold on;
plot(an2);