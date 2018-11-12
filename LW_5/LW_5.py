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
    input_path = './input_5.txt'
    output_path = './output_5.csv'
    try:
        (a, b, c) = (float(entry_a.get()), float(entry_b.get()), float(entry_c.get()))
        (t, n, k) = (float(entry_t.get()), int(entry_n.get()), int(entry_k.get()))
        num_diff = values.index(combobox.get())
    except ValueError:
        label.configure(text = 'a: float, b: float, c: float, n: int, k: int, t: float')
        return

    with open(input_path, 'w') as file:
        print(a, b, c, t, n, k, num_diff, file = file)
    os.system('./LW_5 < ' + input_path + ' > ' + output_path)
    with open(output_path, 'r') as file:
        reader = csv.reader(file)
        x, u = [], ([], [], [], [])
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

combobox = ttk.Combobox(master, values = values)
combobox.set(values[0])
combobox.grid(row = 3, column = 0, columnspan = 6)

button_apply = tk.Button(master, text = 'Построить график')
button_apply.grid(row = 4, column = 0, columnspan = 6)

button_apply.bind('<Button-1>', LW_5)
master.bind('<Return>', LW_5)

master.mainloop()
