import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import tkinter as tk
import tkinter.ttk as ttk

values = [u"Двухточечная первого",
          u"Двухточечная второго",
          u"Трехточечная второго"]

def LW_5(event):
    l = 3.1415926535897932384626
    error_path = './error_5.csv'
    input_path, output_path = './input_5.txt', './output_5.csv'

    try:
        a, b, c = (float(entry_a.get()), float(entry_b.get()), float(entry_c.get()))
        t, n, k = (float(entry_t.get()), int(entry_n.get()), int(entry_k.get()))
        boundary = values.index(combobox_boundary.get())
    except ValueError:
        label.configure(text = 'a: float, b: float, c: float, t: float, n: int, k: int')
        return
    h, tau = l / n, t / k
    sigma = a * a * tau / (h * h)
    if (sigma > 0.5):
        label.configure(text = 'sigma > 0.5')
        return

    with open(input_path, 'w') as file:
        print(a, b, c, t, n, k, boundary, file = file)
    os.system('./LW_5 {} < {} > {}'.format(error_path, input_path, output_path))

    x, u = [], ([], [], [], [])
    with open(output_path, 'r') as file:
        reader = csv.reader(file)
        for line in reader:
            x.append(float(line[0]))
            for i in range(0, 4):
                u[i].append(float(line[i + 1]))

    fig = plt.figure()
    subplot = fig.add_subplot(111, facecolor = '#FFFFFF')
    subplot.plot(x, u[0], color = 'red', lw = 2, label = 'u')
    subplot.plot(x, u[1], color = 'green', lw = 2, label = 'explicit')
    subplot.plot(x, u[2], color = 'blue', lw = 2, label = 'implicit')
    subplot.plot(x, u[3], color = 'cyan', lw = 2, label = 'crank_nicolson')
    plt.title('sigma = ' + str(sigma))
    plt.legend()

    t, error = [], ([], [], [])
    with open(error_path, 'r') as file:
        reader = csv.reader(file)
        for line in reader:
            t.append(float(line[0]))
            for i in range(0, 3):
                error[i].append(float(line[i + 1]))

    fig = plt.figure()
    subplot = fig.add_subplot(111, facecolor = '#FFFFFF')
    subplot.plot(t, error[0], color = 'green', lw = 2, label = 'explicit')
    subplot.plot(t, error[1], color = 'blue', lw = 2, label = 'implicit')
    subplot.plot(t, error[2], color = 'cyan', lw = 2, label = 'crank_nicolson')
    plt.legend()

    plt.show()

master = tk.Tk()

label = tk.Label(master, text = 'Введите коэффициенты:')
label.grid(row = 0, column = 0, columnspan = 6)

tk.Label(master, text='a = ').grid(row = 1, column = 0)
entry_a = tk.Entry(master)
entry_a.grid(row = 1, column = 1)

tk.Label(master, text='b = ').grid(row = 1, column = 2)
entry_b = tk.Entry(master)
entry_b.grid(row = 1, column = 3)

tk.Label(master, text='c = ').grid(row = 1, column = 4)
entry_c = tk.Entry(master)
entry_c.grid(row = 1, column = 5)

tk.Label(master, text='T = ').grid(row = 2, column = 0)
entry_t = tk.Entry(master)
entry_t.grid(row = 2, column = 1)

tk.Label(master, text='N = ').grid(row = 2, column = 2)
entry_n = tk.Entry(master)
entry_n.grid(row = 2, column = 3)

tk.Label(master, text='K = ').grid(row = 2, column = 4)
entry_k = tk.Entry(master)
entry_k.grid(row = 2, column = 5)

tk.Label(master, text='граничные: ').grid(row = 3, column = 2)
combobox_boundary = ttk.Combobox(master, values = values)
combobox_boundary.set(values[0])
combobox_boundary.grid(row = 3, column = 3)

button_apply = tk.Button(master, text = 'Построить график')
button_apply.grid(row = 4, column = 0, columnspan = 6)

button_apply.bind('<Button-1>', LW_5)
master.bind('<Return>', LW_5)

master.mainloop()
