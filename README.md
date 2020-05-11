# gii-minimize

gii (global instability index) minimize is an experimental project attempting use simple concepts of bond strain and electrostatics with generic scipy minimize methods to optimize inorganic structures. Ad-hoc parameters used bond-valence [1], and shannon radii [2]. 

Good results have been found for rotationally distorted perovskites such as LiNbO3 (a-a-a- glazer tilt) and CaTiO3 (a-a-c+ glazer tilt). Minimizer also identifies local distortions such as the Jahn-Teller distortion in LaMnO3 and polar distortion in tetragonal BaTiO3*

Improvements needed:
i. automatically identify point group operations
ii. simple io... this version is a poc

*Additional constraint needed
1. I. D. Brown, The Chemical Bond in Inorganic Chemistry. "The Bond Valence Method" Oxford University Press, Oxford, 2002
2. R. D. Shannon and C. T. Prewitt, "Revised values of effective ionic radii", Acta Cryst. (1970). B26, 1046-1048
