# Trying to make everything from the reactions only
# Kinases Model II
# Time


import numpy as np
dt = 0.01
t = np.arange(0, 1, dt)

# Parameters

user_input = ['Ssk1', 1,
              'Ssk2', 5.3,
              'Ssk2P', 0,
              'Pbs2', 42,
              'Pbs2P', 8.3,
              'Pbs2PP', 0,
              'Hog1', 79,
              'Hog1P', 0.9,
              'Hog1PP', 0, ]

k = (1.2, 1.438, 0.011, 1.438, 0.011, 1.438, 0.011, 1.438, 0.011, 1.438, 0.011,
     0.29, 0.053, 0.11, 0.091)

# Empty reaction ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (

    (1, 'Ssk1'), (), k[0],

    (-1, 'Ssk1', 1, 'Ssk2'), (1, 'Ssk2P'), k[1],

    (1, 'Ssk2P'), (1, 'Ssk2'), k[2],

    (-1, 'Ssk2P', 1, 'Pbs2'), (1, 'Pbs2P'), k[3],

    (1, 'Pbs2P'), (1, 'Pbs2'), k[4],

    (-1, 'Ssk2P', 1, 'Pbs2P'), (1, 'Pbs2PP'), k[5],

    (1, 'Pbs2PP'), (1, 'Pbs2P'), k[6],

    (-1, 'Pbs2PP', 1, 'Hog1'), (1, 'Hog1P'), k[7],

    (1, 'Hog1P'), (1, 'Hog1'), k[8],

    (-1, 'Pbs2PP', 1, 'Hog1P'), (1, 'Hog1PP'), k[9],

    (1, 'Hog1PP'), (1, 'Hog1P'), k[10])

(solution, names) = Stochy2(user_input, reactions, t)

fig, ax = plt.subplots()

for i in [0, 5, 8]:
    ax.plot(t, solution[:, i], label=names[i])
legend = ax.legend(loc='center', shadow=True)
# Show over time
plt.show()