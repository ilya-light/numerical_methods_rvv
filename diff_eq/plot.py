import matplotlib.pyplot as plt

def read_data_from_file(filename):
    data = {}
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line:  # Проверка на пустую строку
                name, x, y = line.split(',')
                x = int(x)
                y = float(y)
                if name not in data:
                    data[name] = ([], [])  # Инициализация списков для x и y
                data[name][0].append(x)  # Добавление x
                data[name][1].append(y)  # Добавление y
    return data

# Чтение данных из двух файлов
data1 = read_data_from_file("results_vector.txt")
data2 = read_data_from_file("results_scalar.txt")

# Уникальные названия методов
methods = set(data1.keys()).union(set(data2.keys()))

# Создание подграфиков
fig, axs = plt.subplots(len(methods), 1, figsize=(10, 15))

# Построение графиков для каждого метода
for i, method in enumerate(methods):
    if method in data1 and method in data2:
        x_values1, y_values1 = data1[method]
        x_values2, y_values2 = data2[method]
        
        axs[i].plot(x_values1, y_values1, marker='o', label=f"{method} (File 1)")
        axs[i].plot(x_values2, y_values2, marker='x', label=f"{method} (File 2)")
        
        axs[i].set_title(f"Сравнение метода {method}", fontsize=14)  # Уменьшение размера шрифта
        axs[i].set_xlabel("Количество переменных (N)", fontsize=12)
        axs[i].set_ylabel("Время выполнения (секунды)", fontsize=12)
        axs[i].legend()
        axs[i].grid()
        axs[i].set_xscale('log')  # Логарифмическая шкала по оси x
        axs[i].set_yscale('linear')  # Линейная шкала по оси y

# Настройка общего оформления
plt.tight_layout(pad=3.0)  # Увеличение отступов между подграфиками

# Сохранение графиков в файл
plt.savefig("results_comparison_plots.png")

# Отображение графиков
plt.show()