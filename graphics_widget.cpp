#include "graphics_widget.h"
#include "MainWindow.h"
#include <math.h>

#define NUM_SIDES 32

GraphicsWidget::GraphicsWidget(QWidget* parent, Simulation* s) : QGLWidget(parent)
{
	zoom = 1.0;
	center[0] = 0.0;
	center[1] = 0.0;
	trackShip = false;

	sim = s;
}

GraphicsWidget::~GraphicsWidget()
{
}

/*************************************************************************/

void GraphicsWidget::initializeGL()
{
	makeCurrent();
	glDisable(GL_LIGHTING);
	glClearColor(1.0,1.0,1.0,1.0);
}

void GraphicsWidget::resizeGL(int w, int h)
{
	if(w==0) w=1;
	if(h==0) h=1;
	glViewport(0,0,w,h);
	setView(center[0], center[1], zoom);
}

void GraphicsWidget::paintGL()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glDisable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT);

	//draw frame
	glBegin(GL_LINE_LOOP);
		glColor3d(0.0,0.0,0.0);
		glVertex2d(center[0]-3.99/zoom,center[1]+3.99*height()/width()/zoom);
		glVertex2d(center[0]+3.99/zoom,center[1]+3.99*height()/width()/zoom);
		glVertex2d(center[0]+3.99/zoom,center[1]-3.99*height()/width()/zoom);
		glVertex2d(center[0]-3.99/zoom,center[1]-3.99*height()/width()/zoom);
	glEnd();

	CBody* ship = sim->getShip();
	vector4 pos = sim->geometry()->convertTo(ship->getVector(0),COORD_CARTESIAN);

	if(trackShip)
		glTranslated(center[0]-pos[1], center[1]-pos[2], 0.0);
	drawBlackHole();
	drawShip();

	glFlush();
}

void GraphicsWidget::mousePressEvent(QMouseEvent* event)
{
	lastPos = event->pos();
}

void GraphicsWidget::mouseMoveEvent(QMouseEvent* event)
{
	int dx = event->x() - lastPos.x();
	int dy = event->y() - lastPos.y();
	lastPos = event->pos();

	if(event->buttons() & Qt::LeftButton && !trackShip)
	{
		center[0] -= 8.0*(double)dx/(width()*zoom);
		center[1] += 8.0*(double)dy/(width()*zoom);
		setView(center[0], center[1], zoom);
		updateGL();
	}
	if(event->buttons() & Qt::RightButton)
	{
		zoom *= exp(-(double)dy/height());
		setView(center[0], center[1], zoom);
		updateGL();
	}
}

/*************************************************************************/

void GraphicsWidget::setView(double x, double y, double z)
{
	center[0] = x;
	center[1] = y;
	zoom = z;

	double w,h;
	w = 4.0 / zoom;
	h = w*(double)height()/(double)width();

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(x-w, x+w, y-h, y+h, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void GraphicsWidget::setTrackShip(bool a)
{
	trackShip = a;
}

void GraphicsWidget::drawBlackHole()
{
	double M = sim->geometry()->getM();
	double a = sim->geometry()->getA();
	double r = M+sqrt(M*M-a*a);
	int i;

	glBegin(GL_TRIANGLE_FAN);
		glColor3d(0.0,0.0,0.0);
		glVertex2d(0.0,0.0);
		double phi;
		for(i=0; i<NUM_SIDES+1; i++)
		{
			phi = 2.0*M_PI/(double)NUM_SIDES*i;
			glVertex2d(r*cos(phi), r*sin(phi));
		}
	glEnd();
}

void GraphicsWidget::drawShip()
{
	CBody* ship = sim->getShip();
	if(ship==NULL) return;

	vector4 pos;
	vector3 front;
	int i;

	pos = sim->geometry()->convertTo(ship->getVector(0),COORD_CARTESIAN);
	front = ship->getCartesianDir(ship->getVector(3));
	double angle = atan2(front[1],front[0]);

	//ship
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslated(pos[1],pos[2],0.0);
	glRotated(angle*180/M_PI,0.0,0.0,1.0);
	glBegin(GL_TRIANGLE_FAN);
		glColor3d(0.0,1.0,0.0);
		glVertex2d(0.1,0.0);
		glVertex2d(0.0,-0.03);
		glVertex2d(-0.03,0.0);
		glVertex2d(0.0,0.03);
	glEnd();
	glPopMatrix();

	//path
	QVector<vector4> path = sim->getPath();
	glBegin(GL_LINE_STRIP);
		glColor3d(0.0,1.0,0.0);
		for(i=0; i<path.size(); i++)
		{
			vector4 point = sim->geometry()->convertTo(path[i],COORD_CARTESIAN);
			glVertex2d(point[1],point[2]);
		}
		glVertex2d(pos[1],pos[2]);
	glEnd();
}

/*************************************************************************/
