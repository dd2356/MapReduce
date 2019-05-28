import matplotlib.pyplot as plt
import numpy as np

times = [[None for j in range(2**i)] for i in range(4)]
times[0][0] = {"read": 2.85, "map": 19.28, "shuffle": 0.16, "communicate": 0.02, "reduce": 0.67}

times[1][0] = {"read": 1.47, "map": 10.23, "shuffle": 0.09, "communicate": 0.28, "reduce": 0.31}
times[1][1] = {"read": 1.48, "map": 10.20, "shuffle": 0.09, "communicate": 0.29, "reduce": 0.31}

times[2][0] = {"read": 0.87, "map": 5.60, "shuffle": 0.05, "communicate": 0.46, "reduce": 0.15}
times[2][1] = {"read": 0.86, "map": 5.62, "shuffle": 0.05, "communicate": 0.40, "reduce": 0.15}
times[2][2] = {"read": 0.86, "map": 5.59, "shuffle": 0.05, "communicate": 0.44, "reduce": 0.15}
times[2][3] = {"read": 0.86, "map": 5.52, "shuffle": 0.05, "communicate": 0.51, "reduce": 0.15}

times[3][0] = {"read": 0.49, "map": 4.33, "shuffle": 0.04, "communicate": 0.89, "reduce": 0.09}
times[3][1] = {"read": 0.50, "map": 4.48, "shuffle": 0.04, "communicate": 0.68, "reduce": 0.10}
times[3][2] = {"read": 0.49, "map": 4.44, "shuffle": 0.04, "communicate": 0.85, "reduce": 0.09}
times[3][3] = {"read": 0.51, "map": 4.25, "shuffle": 0.04, "communicate": 0.91, "reduce": 0.09}
times[3][4] = {"read": 0.52, "map": 4.45, "shuffle": 0.03, "communicate": 0.81, "reduce": 0.09}
times[3][5] = {"read": 0.51, "map": 4.52, "shuffle": 0.03, "communicate": 0.65, "reduce": 0.09}
times[3][6] = {"read": 0.51, "map": 4.45, "shuffle": 0.03, "communicate": 0.58, "reduce": 0.09}
times[3][7] = {"read": 0.51, "map": 4.49, "shuffle": 0.04, "communicate": 0.82, "reduce": 0.09}

for key in times[0][0].keys():
	r = np.arange(len(times))
	for i in r:
		for j in range(2**i):
			plt.loglog(2**i, times[i][j][key], "*")
	plt.loglog(2**r, times[0][0][key] / 2**r)

	plt.title(key)
	plt.show()


