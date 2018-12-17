import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import PhotoImage

def LW_7(event):
    error_path = './error_7.csv'
    input_path, output_path = './input_7.txt', './output_7.csv'

    try:
        n_1, n_2 = int(entry_n1.get()), int(entry_n2.get())
        epsilon, y, omega = float(entry_epsilon.get()), float(entry_y.get()), float(entry_omega.get())
    except ValueError:
        label.configure(text = 'n_1: int, n_2: int, epsilon: float, y: float, omega: float')
        return

    with open(input_path, 'w') as file:
        print(n_1, n_2, epsilon, y, omega, file = file)
    os.system('./LW_7 {} < {} > {}'.format(error_path, input_path, output_path))

    x, u = [], ([], [], [])
    with open(output_path, 'r') as file:
        reader = csv.reader(file)
        for line in reader:
            x.append(float(line[0]))
            for i in range(0, 3):
                u[i].append(float(line[i + 1]))

    fig = plt.figure()
    subplot = fig.add_subplot(111, facecolor = '#FFFFFF')
    subplot.plot(x, u[0], color = 'red', lw = 2, label = 'u')
    subplot.plot(x, u[1], color = 'green', lw = 2, label = 'successive')
    subplot.plot(x, u[2], color = 'blue', lw = 2, label = 'seidel')
    plt.legend()

    t, error = [], ([], [])
    with open(error_path, 'r') as file:
        reader = csv.reader(file)
        for line in reader:
            t.append(float(line[0]))
            for i in range(0, 2):
                error[i].append(float(line[i + 1]))

    fig = plt.figure()
    subplot = fig.add_subplot(111, facecolor = '#FFFFFF')
    subplot.plot(t, error[0], color = 'green', lw = 2, label = 'successive')
    subplot.plot(t, error[1], color = 'blue', lw = 2, label = 'seidel')
    plt.legend()

    plt.show()

master = tk.Tk()

photo_image = PhotoImage(file = 'LW_7.png')

label_photo = tk.Label(master, image = photo_image)
label_photo.image = photo_image
label_photo.grid(row = 0, column = 0, columnspan = 4)

label = tk.Label(master, text = 'Введите коэффициенты:')
label.grid(row = 1, column = 0, columnspan = 4)

tk.Label(master, text = 'N1 = ').grid(row = 2, column = 0)
entry_n1 = tk.Entry(master)
entry_n1.grid(row = 2, column = 1)

tk.Label(master, text = 'N2 = ').grid(row = 2, column = 2)
entry_n2 = tk.Entry(master)
entry_n2.grid(row = 2, column = 3)

tk.Label(master, text = 'epsilon = ').grid(row = 3, column = 0)
entry_epsilon = tk.Entry(master)
entry_epsilon.grid(row = 3, column = 1)

tk.Label(master, text = 'y = ').grid(row = 3, column = 2)
entry_y = tk.Entry(master)
entry_y.grid(row = 3, column = 3)

tk.Label(master, text = 'omega = ').grid(row = 3, column = 4)
entry_omega = tk.Entry(master)
entry_omega.grid(row = 3, column = 5)

button_apply = tk.Button(master, text = 'Построить график')
button_apply.grid(row = 4, column = 0, columnspan = 4)

button_apply.bind('<Button-1>', LW_7)
master.bind('<Return>', LW_7)

master.mainloop()
