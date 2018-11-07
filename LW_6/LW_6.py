import csv
import matplotlib.pyplot as plt
import numpy as np
import os
from tkinter import *

def LW_6(event):
    input_path = './input_6.txt'
    output_path = './output_6.csv'
    try:
        (a, b, c, d) = (float(entry_a.get()), float(entry_b.get()), float(entry_c.get()), float(entry_d.get()))
        (t, n, k) = (float(entry_t.get()), int(entry_n.get()), int(entry_k.get()))
    except ValueError:
        label.configure(text = 'a: float, b: float, c: float, d: float, n: int, k: int, t: float')
        return

    with open(input_path, 'w') as file:
        print(a, b, c, d, t, n, k, file = file)
    os.system('./LW_6 < ' + input_path + ' > ' + output_path)
    with open(output_path, 'r') as file:
        reader = csv.reader(file)
        x, u = [], ([], [], [], [])
        for line in reader:
            x.append(float(line[0]))
            for i in range(0, 3):
                u[i].append(float(line[i + 1]))

        fig = plt.figure()
        subplot = fig.add_subplot(111, facecolor = '#FFFFFF')
        subplot.plot(x, u[0], color = 'red', lw = 2, label = 'u')
        subplot.plot(x, u[1], color = 'green', lw = 2, label = 'explicit')
        subplot.plot(x, u[2], color = 'blue', lw = 2, label = 'implicit')
        plt.legend()
        plt.show()

tk = Tk()

label = Label(tk, text = 'Введите коэффициенты:')
label.grid(row = 0, column = 0, columnspan = 8)

Label(tk, text='a = ').grid(row = 1, column = 0)
entry_a = Entry(tk)
entry_a.grid(row = 1, column = 1)

Label(tk, text='b = ').grid(row = 1, column = 2)
entry_b = Entry(tk)
entry_b.grid(row = 1, column = 3)

Label(tk, text='c = ').grid(row = 1, column = 4)
entry_c = Entry(tk)
entry_c.grid(row = 1, column = 5)

Label(tk, text='d = ').grid(row = 1, column = 6)
entry_d = Entry(tk)
entry_d.grid(row = 1, column = 7)

Label(tk, text='T = ').grid(row = 2, column = 0)
entry_t = Entry(tk)
entry_t.grid(row = 2, column = 1)

Label(tk, text='N = ').grid(row = 2, column = 2)
entry_n = Entry(tk)
entry_n.grid(row = 2, column = 3)

Label(tk, text='K = ').grid(row = 2, column = 4)
entry_k = Entry(tk)
entry_k.grid(row = 2, column = 5)

button_apply = Button(tk, text = 'Построить график')
button_apply.grid(row = 3, column = 0, columnspan = 8)

button_apply.bind('<Button-1>', LW_6)
tk.bind('<Return>', LW_6)

tk.mainloop()
