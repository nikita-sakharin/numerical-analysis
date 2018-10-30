#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <cstddef>

#include <iomanip>
#include <iostream>
#include <limits>

#include <QGraphicsScene>
#include <QMainWindow>

#include "../generic/header.hpp"
#include "parabolic_pde.hpp"

namespace Ui
{
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    void resizeEvent(QResizeEvent *) override;
    ~MainWindow() override;

private slots:
    void on_applyArg_clicked();

private:
    static constexpr int INDENT_10 = 10,
        INDENT_20 = 20,
        INDENT_3 = 3,
        INDENT_BOTTOM = 64,
        INDENT_TOP = 40;

    Ui::MainWindow *ui;

    void widgetMove(const QResizeEvent *);
    template<typename T>
    void vectorAddLine(const QPen &,
        const T, const ublas::vector<T> &);
    template<typename T>
    void functionAddLine(const QPen &, const T, const T, std::size_t, const T,
                         T (const T &, const T &, const T &));
};

#endif // MAINWINDOW_H
