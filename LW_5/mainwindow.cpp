#include <cmath>

#include <QResizeEvent>

#include "mainwindow.h"
#include "ui_mainwindow.h"

static constexpr ldbl u_0_t(const ldbl &, const ldbl &) noexcept;
static constexpr ldbl u_l_t(const ldbl &, const ldbl &) noexcept;
static constexpr ldbl u_x_0(const ldbl &, const ldbl &) noexcept;
static constexpr ldbl u_exact(const ldbl &, const ldbl &, const ldbl &) noexcept;

static constexpr ldbl L = PI_LDBL;

MainWindow::MainWindow(QWidget * const parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QGraphicsScene * const scene = new QGraphicsScene(ui->graphicsOfU);
    ui->graphicsOfU->setScene(scene);
}

MainWindow::~MainWindow()
{
    delete ui->graphicsOfU->scene();
    delete ui;
}

void MainWindow::resizeEvent(QResizeEvent * const event)
{
    QMainWindow::resizeEvent(event);
    widgetMove(event);
    ui->graphicsOfU->fitInView(ui->graphicsOfU->scene()->sceneRect());
}

void MainWindow::on_applyArg_clicked()
{
    QGraphicsScene * const scene = ui->graphicsOfU->scene();
    scene->clear();

    const ldbl a = static_cast<ldbl>(ui->fieldA->text().toDouble()),
              t = static_cast<ldbl>(ui->fieldT->text().toDouble());
    const std::size_t
              n = static_cast<std::size_t>(ui->fieldN->text().toULong()),
              k = static_cast<std::size_t>(ui->fieldK->text().toULong());
    const QPen darkGrayPen(QBrush(Qt::darkGray), 0),
            redPen(QBrush(Qt::red), 0),
            greenPen(QBrush(Qt::green), 0),
            bluePen(QBrush(Qt::blue), 0),
            cyanPen(QBrush(Qt::cyan), 0);
    functionAddLine(redPen, a, L, n, t, u_exact);
    vectorAddLine(greenPen, L,
        explicit_fdm<ldbl>(a, L, n, t, k, u_0_t, u_l_t, u_x_0));
    vectorAddLine(bluePen, L,
        implicit_fdm<ldbl>(a, L, n, t, k, u_0_t, u_l_t, u_x_0));
    vectorAddLine(cyanPen, L,
        crank_nicolson<ldbl>(a, L, n, t, k, u_0_t, u_l_t, u_x_0));
    scene->setSceneRect(scene->itemsBoundingRect());
    ui->graphicsOfU->fitInView(scene->sceneRect());

    scene->addLine(scene->sceneRect().x(), 0.0,
                   scene->sceneRect().x() + scene->sceneRect().width(), 0.0, darkGrayPen);
    scene->addLine(0.0, scene->sceneRect().y(),
                   0.0, scene->sceneRect().y() + scene->sceneRect().height(), darkGrayPen);
}

void MainWindow::widgetMove(const QResizeEvent * const event)
{
    const QSize argSize = ui->argA->size(), fieldSize = ui->fieldA->size(),
                applySize = ui->applyArg->size(), eventSize = event->size();
    ui->applyArg->move((eventSize.width() - applySize.width()) / 2,
                       eventSize.height() - INDENT_BOTTOM - applySize.height());
    ui->fieldN->move((eventSize.width() - INDENT_10) / 2 - fieldSize.width(),
                     ui->applyArg->y() - fieldSize.height() - INDENT_20);
    ui->argN->move(ui->fieldN->x() - INDENT_10 - argSize.width(),
                   ui->fieldN->y() + INDENT_3);
    ui->argK->move((eventSize.width() + INDENT_10) / 2,
                   ui->argN->y());
    ui->fieldK->move(ui->argK->x() + argSize.width() + INDENT_10,
                     ui->fieldN->y());
    ui->argA->move(ui->argN->x(),
                   ui->argN->y() - INDENT_10 - argSize.height());
    ui->fieldA->move(ui->fieldN->x(),
                     ui->fieldN->y() - INDENT_10 - fieldSize.height());
    ui->argT->move(ui->argK->x(), ui->argA->y());
    ui->fieldT->move(ui->fieldK->x(), ui->fieldA->y());
    ui->graphicsOfU->resize(eventSize.width() - INDENT_10 - INDENT_10,
                            ui->fieldA->y() - INDENT_TOP - INDENT_20);
}

template<typename T>
void MainWindow::vectorAddLine(const QPen &pen,
    const T l, const ublas::vector<T> &lineVector)
{
    QGraphicsScene * const scene = ui->graphicsOfU->scene();
    const std::size_t n_upper = lineVector.size() - 1;
    const T h = l / n_upper;
    for (std::size_t i = 0; i < n_upper; ++i)
    {
        scene->addLine(static_cast<qreal>(h * i), static_cast<qreal>(lineVector[i]),
                       static_cast<qreal>(h * (i + 1)), static_cast<qreal>(lineVector[i + 1]), pen);
    }
}

template<typename T>
void MainWindow::functionAddLine(const QPen &pen, const T a,
    const T l, const std::size_t n_upper, const T t,
    T u(const T &, const T &, const T &))
{
    QGraphicsScene * const scene = ui->graphicsOfU->scene();
    const T h = l / n_upper;
    for (std::size_t i = 0; i < n_upper; ++i)
    {
        scene->addLine(static_cast<qreal>(h * i),
                       static_cast<qreal>(u(a, h * i, t)),
                       static_cast<qreal>(h * (i + 1)),
                       static_cast<qreal>(u(a, h * (i + 1), t)),
                       pen);
    }
}

static constexpr ldbl u_0_t(const ldbl &a, const ldbl &t) noexcept
{
    return std::exp(-a * t);
}

static constexpr ldbl u_l_t(const ldbl &a, const ldbl &t) noexcept
{
    return -std::exp(-a * t);
}

static constexpr ldbl u_x_0(const ldbl &, const ldbl &x) noexcept
{
    return std::cos(x);
}

static constexpr ldbl u_exact(const ldbl &a, const ldbl &x, const ldbl &t) noexcept
{
    return std::exp(-a * t) * std::cos(x);
}
