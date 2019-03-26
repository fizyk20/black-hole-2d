#ifndef __GRAPHICS_WIDGET__
#define __GRAPHICS_WIDGET__

#include <QGLWidget>
#include <QMouseEvent>
#include "simulation.h"

class GraphicsWidget : public QGLWidget
{
Q_OBJECT
	double center[2];
	double zoom;
	bool trackShip;

	QPoint lastPos;
	Simulation* sim;
protected:
	void initializeGL();
	void resizeGL(int w, int h);
	void paintGL();
	void mousePressEvent(QMouseEvent*);
	void mouseMoveEvent(QMouseEvent*);

public:
	GraphicsWidget(QWidget* parent, Simulation*);
	~GraphicsWidget();

	void setView(double, double, double); //center_x, center_y, zoom

	void setTrackShip(bool);
	void drawBlackHole();
	void drawShip();
};

#endif
