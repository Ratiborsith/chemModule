import re

# Sample text containing IUPAC names
text = """
3-phenyl-1-(thiophen-2-yl)prop-2-en-1-one was prepared from benzaldehyde and 1-(thiophen-2-yl)ethanone via aldol condensation using procedure described by Azam (Parveen, H.; Iqbal, P. F.; Azam, A. Synth. Commu., 2008, 38, 3973). 1H NMR (400 MHz, CDCl3) δ 7.88-7.80 (m, 2H), 7.67 (dd, J=4.9, 1.1 Hz, 1H), 7.66-7.59 (m, 2H), 7.47-7.34 (m, 4H), 7.18 (dd, J=5.0, 3.8 Hz, 1H). ESI-MS (m/z): 215 [M+H]+.
Figure US20230136730A1-20230504-C00154
[0414]
4-phenyl-6-(thiophen-2-yl)-2-thioxo-1,2-dihydropyridine-3-carbonitrile. To a solution of 3-phenyl-1-(thiophen-2-yl)prop-2-en-1-one (2.34 mmol, 500 mg) and cyanothioacetamide (7.0 mmol, 717 mg, 3.0 equiv.) in ethanol (7 mL), a few drops of piperidine were added. The reaction was refluxed for 3 h. The solid that formed was collected and recrystallized from acetic acid to give designed product in 46% isolated yield. 1H NMR (400 MHz, DMSO-d6) δ 8.17 (d, J=3.8 Hz, 1H), 7.96 (d, J=5.0 Hz, 1H), 7.74-7.62 (m, 2H), 7.54 (dd, J=5.1, 2.0 Hz, 3H), 7.31-7.19 (m, 1H), 7.01 (s, 1H). ESI-MS (m/z): 295 [M+H]+.
Figure US20230136730A1-20230504-C00155
[0415]
2-(((butylthio)methyl)sulfinyl)-4-phenyl-6-(thiophen-2-yl)nicotinonitrile. Acetic Acid (900 μL) and hydrogen peroxide (0.57 mmol, 1.5 equiv., 30% solution in water) were added to the solution of 2-(((butylthio)methyl)sulfinyl)-4-phenyl-6-(thiophen-2-yl)nicotinonitrile (0.38 mmol, 150 mg) in chloroform (900 μL). The reaction mixture was stirring at 32° C. for 45 min. The reaction was then diluted with EtOAc and washed with saturated NaHCO3 solution, dried over magnesium sulfate, filtered and concentrated under reduced pressure to give 153 mg of designed product (98%). 1H NMR (400 MHz, CDCl3) δ 7.75 (dd, J=3.8, 1.1 Hz, 1H), 7.66-7.57 (m, 2H), 7.58-7.51 (m, 4H), 7.47 (s, 1H), 7.16 (dd, J=5.0, 3.8 Hz, 1H), 4.74 (d, J=13.0 Hz, 1H), 4.41 (d, J=13.0 Hz, 1H), 2.97 (dt, J=13.0, 8.2 Hz, 1H), 2.81 (dt, J=12.9, 7.3 Hz, 1H), 1.94-1.76 (m, 2H), 1.53-1.38 (m, 2H), 0.94 (t, J=7.4 Hz, 3H). ESI-MS (m/z): 413 [M+H]+.
Figure US20230136730A1-20230504-C00156
[0416]
SW033291 2-(butylsulfinyl)-4-phenyl-6-(thiophen-2-yl)thieno[2,3-b]pyridin-3-amine was prepared using procedure describe by Kalugin (Kalugin V. E. Russian. Chem. Bull., Int. Ed., 2006, 55, 529). To the solution of 4-(((butylthio)methyl)sulfinyl)-2,6-diphenylpyrimidine-5-carbonitrile (0.53 mmol, 220 mg) in DMF (0.25 M)/EtOH (0.5 M) was added KOH (0.32 mmol, 18 mg, 0.6 equiv., 0.1 M in water). The reaction mixture was stirred at 35° C. for 40 min. Once complete, the reaction was diluted with EtOAc and washed with 10% aq. solution of acidic acid, the organic phase was separated and aqueous layer was extracted twice with EtOAc, dried over magnesium sulfate, filtered and concentrated under reduced pressure to give 211 mg of SW033291 2-(butylsulfinyl)-4-phenyl-6-(thiophen-2-yl)thieno[2,3-b]pyridin-3-amine (96%). 1H NMR (400 MHz, CDCl3) δ 7.67-7.60 (m, 1H), 7.57-7.35 (m, 7H), 7.10 (dd, J=5.0, 3.7 Hz, 1H), 4.54 (s, 2H), 3.26 (ddd, J=12.8, 9.1, 6.0 Hz, 1H), 3.09 (ddd, J=12.8, 9.1, 6.6 Hz, 1H), 1.83-1.61 (m, 2H), 1.53-1.38 (m, 2H), 0.93 (t, J=7.3 Hz, 3H). ESI-MS (m/z): 413 [M+H]+.
"""

# Regular expression to match IUPAC names
pattern = r'\b(?:\d*-)?[a-zA-Z\d()-]+(?:-(?:[a-zA-Z\d()-]+|\d+-[a-zA-Z\d]+))*\b'

# Extracting IUPAC names using regex
iupac_names = re.findall(pattern, text)

# Print extracted IUPAC names
for name in iupac_names:
    print(name)
