#include "simulation.h"
#include <QtOpenGL>
#include <QTime>

#define MAX_PATH_POINTS 5000
#define SKIP 2
#define STEP 1

Simulation::Simulation()
{
	time_warp = 1.0;
	skip = 0;
	board_time = false;
	F_pg = F_rt = 0.1;
	angVel = 2.5;

	g = new Geometry(M_KERR_EF);
	g -> setA(0.0);

	ship = new CBody(g->convertFrom(vector4(0.0,1.0,M_PI/2,0.0),COORD_SPHERICAL), vector4(1.0,0.0,0.0,sqrt(0.1)), g);
}

Simulation::~Simulation()
{
	if(ship != NULL) delete ship;
}

long Simulation::getTime()
{
	QTime t = QTime::currentTime();
	return t.msec() + 1000*t.second() + 1000*60*t.minute() + 1000*3600*t.hour();
}

long Simulation::deltaT(long t1, long t2)
{
	return (t2-t1)%(24*3600*1000);
}

Geometry* Simulation::geometry()
{
	return g;
}

/********************************************************************************************/

CBody* Simulation::getShip()
{
	return ship;
}

QVector<vector4> Simulation::getPath()
{
	return path;
}

void Simulation::run()
{
	path.clear();
	ship -> setTau(0.0);

	exec();
}

void Simulation::launch()
{
	time = getTime();
	QTimer::singleShot(0, this, SLOT(simulate()));
	running = true;
}

void Simulation::stop()
{
	running = false;
}

void Simulation::reset()
{
	stop();
	path.clear();

	double r = g->getM()*10.0;
	double uphi = sqrt(g->getM()/r/r/r);
	ship -> setVector(0,g->convertFrom(vector4(0.0,r,M_PI/2,0.0),COORD_SPHERICAL));
	ship -> setVector(1,vector4(1.0,0.0,0.0,uphi));
	ship -> setTau(0.0);
	ship -> generateBasis();
	//emit update();

	time_warp = 1.0;
	emit timeWarp(time_warp);
}

void Simulation::simulate()
{
	double dt = (double)deltaT(time,getTime())/1000.0*time_warp;
	time = getTime();

	vector3 f = vector3(F_rt*(keys[Qt::Key_D]-keys[Qt::Key_A]), F_pg*(keys[Qt::Key_W]-keys[Qt::Key_S]), 0.0);
	if(0==keys[Qt::Key_Shift])
		ship->setForce(f);
	else
		ship->applyForce(f);
	ship->setAngVel(vector3(0.0,0.0,(keys[Qt::Key_U]-keys[Qt::Key_O])*angVel));

	ship->propagate(dt, board_time);

	if(skip%SKIP==0) path.push_back(ship->getVector(0));
	skip++;
	skip = skip%SKIP;
	if(path.size()>MAX_PATH_POINTS) path.pop_front();

	//emit update();

	if(ship->getVector(0)[1] <= 0.2*g->getM())
	{
		stop();
		emit autoStop();
	}

	if(running) QTimer::singleShot(0, this, SLOT(simulate()));
}

void Simulation::keyPressed(int k)
{
	keys[k] = 1;
	if(k==Qt::Key_T)
	{
		time_warp*=2.0;
		emit timeWarp(time_warp);
	}
	if(k==Qt::Key_R)
	{
		time_warp/=2.0;
		emit timeWarp(time_warp);
	}
}

void Simulation::keyReleased(int k)
{
	keys[k] = 0;
}

void Simulation::setPgForce(double f)
{
	F_pg = f;
}

void Simulation::setRtForce(double f)
{
	F_rt = f;
}

void Simulation::setAngVel(double v)
{
	angVel = v;
}

void Simulation::setUseBoardTime(bool b)
{
	board_time = b;
}

bool Simulation::useBoardTime()
{
	return board_time;
}
