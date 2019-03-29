#ifndef __MAINWINDOW__
#define __MAINWINDOW__

#include <QtGui>
#include <QTimer>
#include <QtWidgets>
#include "graphics_widget.h"
#include "body.h"
#include "simulation.h"

class MainWindow : public QMainWindow
{
Q_OBJECT
	QTimer* timer;
	GraphicsWidget* gl;
	Simulation* sim;

	QLineEdit *mass;
	QLineEdit *r, *phi, *x, *y;
	QLineEdit *ur, *uphi, *ux, *uy;

	QLabel *l_angmom;
	QLabel *ext_time, *board_time;
	QLabel *vel, *proper_vel;
	QLabel *force, *thrust_pg, *thrust_rt;

	QCheckBox *useBoardTime,*shipTracking;

	QPushButton *start, *stop, *reset;

	QSlider *F_slider, *angmom;

protected:
	void keyPressEvent(QKeyEvent* event);
	void keyReleaseEvent(QKeyEvent* event);
public:
	MainWindow();
	~MainWindow();

	void createMenus();
	void createComponents();

public slots:
	void massChanged(const QString&);
	void angmomChanged(int);

	void rChanged(const QString&);
	void phiChanged(const QString&);
	void xChanged(const QString&);
	void yChanged(const QString&);
	void urChanged(const QString&);
	void uphiChanged(const QString&);
	void uxChanged(const QString&);
	void uyChanged(const QString&);

	void stateUpdate(QObject* o = NULL);

	void buttonClicked();

	void forceChanged(int);
	void boardTime(int);
	void trackShip(int);
	void timeWarp(double);

	void autoStop();

signals:
	void enable(bool);
	void keyPressed(int);
	void keyReleased(int);
	void forceNew(double);
};

#endif
