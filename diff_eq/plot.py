import matplotlib.pyplot as plt

def read_data_from_file(filename):
    data = {}
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line:
                name, x, y = line.split(',')
                x = int(x)
                y = float(y)
                if name not in data:
                    data[name] = ([], [])
                data[name][0].append(x)
                data[name][1].append(y)
    return data

data1 = read_data_from_file("results_vector.txt")
data2 = read_data_from_file("results_scalar.txt")

methods = set(data1.keys()).union(set(data2.keys()))

fig, axs = plt.subplots(len(methods), 1, figsize=(10, 15))

for i, method in enumerate(methods):
    if method in data1 and method in data2:
        x_values1, y_values1 = data1[method]
        x_values2, y_values2 = data2[method]
        
        axs[i].plot(x_values1, y_values1, marker='o', label=f"{method} (File 1)")
        axs[i].plot(x_values2, y_values2, marker='x', label=f"{method} (File 2)")
        
        axs[i].set_title(f"Сравнение метода {method}", fontsize=14)
        axs[i].set_xlabel("Количество переменных (N)", fontsize=12)
        axs[i].set_ylabel("Время выполнения (секунды)", fontsize=12)
        axs[i].legend()
        axs[i].grid()
        axs[i].set_xscale('log')
        axs[i].set_yscale('linear')

plt.tight_layout(pad=3.0)

plt.savefig("results_comparison_plots.png")

plt.show()
