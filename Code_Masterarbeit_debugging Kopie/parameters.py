ars1 = 1.1576
ars2 = 1.0389

prs1 = 0.0512
prs2 = 4.6659
prs3 = -7.8387
prs4 = 5.4571
prs5 = 0.1098
prs6 = -0.0044
prs7 = 0.4021

ko = 1.0546
k1w = 3.5421
k2w = -0.2786
k1b = 2.2658
k2b = 0.0577

sy= 0.011
lambda_0 = 400

band_wavelengths_sen = {
    'B1':443,
    'B2':490 ,
    'B3': 560,
    'B4': 665,
    'B5': 705,
    'B6': 740,
    'B7': 783,
    'B11': 1610
}
#a*p von https://elib.dlr.de/127859/1/Diss_Gege_1994.pdf
#    a*x=0 https://www.researchgate.net/profile/Thomas-Heege/publication/224794670_Flugzeuggestutzte_Fernerkundung_von_Wasserinhaltsstoffen_am_Bodensee/links/5fcfad44299bf188d403dbcf/Flugzeuggestuetzte-Fernerkundung-von-Wasserinhaltsstoffen-am-Bodensee.pdf
#empirical_coefficients_sen = {
 #   'B1': {'aw':0.0107 , 'bbw': 0.00175, 'Rb_sand':0.33,'Rb_LC':0.0525,'Rb_art': 0.1, 'Rb_ga':0.05, 'Rb_ba':0.025, 'Rb_ra':0.075, 'a*_P':0.014, 'a*_X':0,'b*b_X':0.0086},
    
  #  'B2': {'aw': 0.0181, 'bbw': 0.00115, 'Rb_sand': 0.36,'Rb_LC':0.1,'Rb_art': 0.1, 'Rb_ga':0.09, 'Rb_ba':0.03, 'Rb_ra':0.07, 'a*_P':0.015, 'a*_X':0,'b*b_X':0.0086},
    
  #  'B3': {'aw': 0.0621, 'bbw': 0.00065, 'Rb_sand': 0.42,'Rb_LC':0.15,'Rb_art': 0.1, 'Rb_ga': 0.23, 'Rb_ba': 0.07, 'Rb_ra': 0.07, 'a*_P':0.0136, 'a*_X':0,'b*b_X':0.0086},
    
    #'B4': {'aw': 0.4295, 'bbw': 0.0003, 'Rb_sand': 0.53,'Rb_LC':0.12,'Rb_art': 0.1, 'Rb_ga': 0.12, 'Rb_ba': 0.05, 'Rb_ra': 0.17, 'a*_P':0.011,'a*_X':0,'b*b_X':0.0086},
    
 #   'B5': {'aw':0.70675, 'bbw': 0.00025, 'Rb_sand':0.65 ,'Rb_LC':0.1,'Rb_art': 0.1, 'Rb_ga': 0.2, 'Rb_ba':0.1, 'Rb_ra': 0.22, 'a*_P':0.0075,'a*_X':0,'b*b_X':0.0086},
    
   #  'B6': {'aw':2.5319, 'bbw':0.0002, 'Rb_sand':0.66 ,'Rb_LC':0.11,'Rb_art': 0.1, 'Rb_ga': 0.45, 'Rb_ba':0.3, 'Rb_ra':0.5 , 'a*_P':0.0025,'a*_X':0,'b*b_X':0.0086},
    
   # 'B7': {'aw':2.6254 , 'bbw': 0.00015, 'Rb_sand': 0.67, 'Rb_LC':0.11,'Rb_art': 0.1,'Rb_ga':0.45, 'Rb_ba':0.3, 'Rb_ra':0.5, 'a*_P':0.001, 'a*_X':0.0,'b*b_X':0.0086},
    
   # 'B11': {'aw': 726.99, 'bbw': 0.000001, 'Rb_sand': 0.67, 'Rb_art': 0.1,'Rb_LC':0.11,'Rb_ga': 0.45, 'Rb_ba': 0.3, 'Rb_ra': 0.5, 'a*_P':0.001, 'a*_X':0.0,'b*b_X':0.0086}
#}

#mit werten aus dem patent: 
#https://www.researchgate.net/publication/248976578_Three-component_model_of_ocean_colour_and_its_application_to_remote_sensing_of_phytoplankton_pigments_in_coastal_waters bbx
#phytoplankton, manual wasi 
#gelbstoff wasi 2d
#bottom wasi 2d
empirical_coefficients_sen = {
    'B1': {'aw':0.0107 , 'bbw': 0.00175, 'Rb_sand':0.33,'Rb_LC':0.0525,'Rb_art': 0.1, 'Rb_ga':0.05, 'Rb_ba':0.025, 'Rb_ra':0.075, 'a*_P':0.0331, 'a*_Y':0.9661,'b*b_X':0.015},
    
    'B2': {'aw': 0.0181, 'bbw': 0.00115, 'Rb_sand': 0.36,'Rb_LC':0.1,'Rb_art': 0.1, 'Rb_ga':0.09, 'Rb_ba':0.03, 'Rb_ra':0.07, 'a*_P':0.0254, 'a*_Y':0.5978,'b*b_X':0.015},
    
    'B3': {'aw': 0.0621, 'bbw': 0.00065, 'Rb_sand': 0.42,'Rb_LC':0.15,'Rb_art': 0.1, 'Rb_ga': 0.23, 'Rb_ba': 0.07, 'Rb_ra': 0.07,'a*_P':0.0136, 'a*_Y':0.3411,'b*b_X':0.015},
    
    'B4': {'aw': 0.4295, 'bbw': 0.0003,'Rb_sand': 0.53,'Rb_LC':0.12,'Rb_art': 0.1, 'Rb_ga': 0.12, 'Rb_ba': 0.05, 'Rb_ra': 0.17, 'a*_P':0.0166,'a*_Y':0.1808,'b*b_X':0.015},
    
    'B5': {'aw':0.70675, 'bbw': 0.00025, 'Rb_sand':0.65 ,'Rb_LC':0.1,'Rb_art': 0.1, 'Rb_ga': 0.2, 'Rb_ba':0.1, 'Rb_ra': 0.22, 'a*_P':0.00713,'a*_Y':0.1482,'b*b_X':0.015},
    
    'B6': {'aw':2.5319, 'bbw':0.0002, 'Rb_sand':0.66 ,'Rb_LC':0.11,'Rb_art': 0.1, 'Rb_ga': 0.45, 'Rb_ba':0.3, 'Rb_ra':0.5 , 'a*_P':0.000549,'a*_Y':0.1264,'b*b_X':0.015},
    
    'B7': {'aw':2.6254 , 'bbw': 0.00015, 'Rb_sand': 0.67, 'Rb_LC':0.11,'Rb_art': 0.1,'Rb_ga':0.45, 'Rb_ba':0.3, 'Rb_ra':0.5, 'a*_P':0, 'a*_Y':0.1056,'b*b_X':0.015},
    
    'B11': {'aw': 726.99, 'bbw': 0.000001, 'Rb_sand': 0.67, 'Rb_art': 0.1,'Rb_LC':0.11,'Rb_ga': 0.45, 'Rb_ba': 0.3, 'Rb_ra': 0.5, 'a*_P':0, 'a*_Y':0.000001,'b*b_X':0.015}
}

band_wavelengths_land = {
    'SR_B1':443,
    'SR_B2':482 ,
    'SR_B3': 560,
    'SR_B4': 655,
    'SR_B5': 842,
    'SR_B6': 1610
}
empirical_coefficients_land = {
    'SR_B1': {'aw':0.0107 , 'bbw': 0.00175, 'Rb_sand':0.33,'Rb_LC':0.0525,'Rb_art': 0.1, 'Rb_ga':0.05, 'Rb_ba':0.025, 'Rb_ra':0.075, 'a*_P':0.0331, 'a*_Y':0.9661,'b*b_X':0.015},
    
    'SR_B2': {'aw': 0.0181, 'bbw': 0.00115, 'Rb_sand': 0.36,'Rb_LC':0.1,'Rb_art': 0.1, 'Rb_ga':0.09, 'Rb_ba':0.03, 'Rb_ra':0.07, 'a*_P':0.0267, 'a*_Y':0.6565,'b*b_X':0.015},
    
    'SR_B3': {'aw': 0.0621, 'bbw': 0.00065, 'Rb_sand': 0.42,'Rb_LC':0.15,'Rb_art': 0.1, 'Rb_ga': 0.23, 'Rb_ba': 0.07, 'Rb_ra': 0.07, 'a*_P':0.0136, 'a*_Y':0.3411,'b*b_X':0.015},
    
    'SR_B4': {'aw': 0.4295, 'bbw': 0.0003, 'Rb_sand': 0.53,'Rb_LC':0.12,'Rb_art': 0.1, 'Rb_ga': 0.12, 'Rb_ba': 0.05, 'Rb_ra': 0.17, 'a*_P':0.0121,'a*_Y':0.1907,'b*b_X':0.015},
    
    'SR_B5': {'aw':4.604, 'bbw':0.000104 , 'Rb_sand': 0.67, 'Rb_LC':0.11,'Rb_art': 0.1,'Rb_ga':0.45, 'Rb_ba':0.3, 'Rb_ra':0.5,'a*_P':0,'a*_Y':0.0847,'b*b_X':0.015},
    
   'SR_B6': {'aw': 726.99, 'bbw': 0.00000001,'Rb_sand': 0.67, 'Rb_art': 0.1,'Rb_LC':0.11,'Rb_ga': 0.45, 'Rb_ba': 0.3, 'Rb_ra': 0.5, 'a*_P':0, 'a*_Y':0.000001,'b*b_X':0.015}
}