# Parameters

user_input = ['A', 100,
              'B', 0]

k = (0.2,0.1)

# Empty reaction ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)
reactions = (
    (1,'A'),(1,'B'),k[0],
    (1,'B'),(1,'A'),k[1],
)

#Function to calculate the possible combinations?

