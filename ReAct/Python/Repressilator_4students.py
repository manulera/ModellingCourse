
import numpy as np
from Gilles import *
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

# Initial conditions
user_input = ['TetR_site', 1,
              'TetR_site_b', 0,
              'TetR_mRNA', 0,
              'TetR_Prot', 0,
              'TetR_Prot2', 0,

              'LacI_site', 0,
              'LacI_site_b', 1,
              'LacI_mRNA', 0,
              'LacI_Prot', 0,
              'LacI_Prot2', 0,

              'Gammacl_site', 0,
              'Gammacl_site_b', 1,
              'Gammacl_mRNA', 0,
              'Gammacl_Prot', 0,
              'Gammacl_Prot2', 0,

              'GFP_site', 0,
              'GFP_site_b', 1,
              'GFP_mRNA', 0,
              'GFP_Prot', 0
              ]

# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)

k = ()
reactions = (

    (-1,'TetR_site'), (1, 'TetR_mRNA'), k[0] , # mRNA transcription

    (1,'TetR_mRNA'),(), k[1], # mRNA degradation

    (-1, 'TetR_mRNA'), (1, 'TetR_Prot'), k[2], # Translation

    (1, 'TetR_Prot'), (), k[3], # Protein degradation

    (2, 'TetR_Prot'), (1, 'TetR_Prot2'), k[4],

    (1, 'TetR_Prot2'), (2, 'TetR_Prot'),  k[5],

    (1, 'TetR_Prot2',1,'Gammacl_site'), (1,'Gammacl_site_b'), k[6], # Binding of the repressor

    (1, 'Gammacl_site_b'), (1, 'TetR_Prot2',1,'Gammacl_site'), k[7], # Unbinding of the repressor

    # ------------------------------------------------------------------------------------------------------------------

    (-1, 'Gammacl_site'), (1, 'Gammacl_mRNA'), k[0],  # mRNA transcription

    (1, 'Gammacl_mRNA'), (), k[1],  # mRNA degradation

    (-1, 'Gammacl_mRNA'), (1, 'Gammacl_Prot'), k[2],  # Translation

    (1, 'Gammacl_Prot'), (), k[3],  # Protein degradation

    (2, 'Gammacl_Prot'), (1, 'Gammacl_Prot2'), k[4],

    (1, 'Gammacl_Prot2'), (2, 'Gammacl_Prot'),  k[5],

    (1, 'Gammacl_Prot2', 1, 'LacI_site'), (1, 'LacI_site_b'), k[6],  # Binding of the repressor

    (1, 'LacI_site_b'), (1, 'Gammacl_Prot2', 1, 'LacI_site'), k[7],  # Unbinding of the repressor

    # ------------------------------------------------------------------------------------------------------------------

    (-1, 'LacI_site'), (1, 'LacI_mRNA'), k[0],  # mRNA transcription

    (1, 'LacI_mRNA'), (), k[1],  # mRNA degradation

    (-1, 'LacI_mRNA'), (1, 'LacI_Prot'), k[2],  # Translation

    (1, 'LacI_Prot'), (), k[3],  # Protein degradation

    (2, 'LacI_Prot'), (1, 'LacI_Prot2'), k[4],

    (1, 'LacI_Prot2'), (2, 'LacI_Prot'),  k[5],

    (1, 'LacI_Prot2', 1, 'TetR_site'), (1, 'TetR_site_b'), k[6],  # Binding of the repressor

    (1, 'TetR_site_b'), (1, 'LacI_Prot2', 1, 'TetR_site'), k[7],  # Unbinding of the repressor

    # ------------------------------------------------------------------------------------------------------------------

    (-1, 'GFP_site'), (1, 'GFP_mRNA'), k[0]*3,  # mRNA transcription

    (1, 'GFP_mRNA'), (), k[1],  # mRNA degradation

    (-1, 'GFP_mRNA'), (1, 'GFP_Prot'), k[2],  # Translation

    (1, 'GFP_Prot'), (), k[3],  # Protein degradation

    (1, 'TetR_Prot2',1,'GFP_site'), (1,'GFP_site_b'), k[4], # Binding of the repressor

    (1, 'GFP_site_b'), (1, 'TetR_Prot2',1,'GFP_site'), k[5], # Unbinding of the repressor

)
# dt is used for the deterministic calculation, and the
dt=0.1
t = np.arange(0, 100, dt)

(solution,(tgill, valsgill, _, _),rows,mode)=ReAct(user_input,reactions,t,mode=1)

Gillesplot(solution,t,tgill, valsgill,rows,mode,which2plot=['TetR_Prot','Gammacl_Prot','LacI_Prot','GFP_Prot'])

#plt.plot(valsgill[0][0,:], valsgill[0][1,:], linewidth=2)

plt.show()