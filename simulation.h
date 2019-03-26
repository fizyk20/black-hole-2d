#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "body.h"
#include <QVector>
#include <QMap>
#include <QThread>
#include <QTimer>
#include <QPointF>

class Simulation : public QThread
{
Q_OBJECT
	Geometry* g;

	CBody* ship;
	double time_warp;
	long time;
	bool running;

	QVector<vector4> path;
	int skip;

	bool board_time;

	QMap<int,int> keys;
	double F_pg, F_rt, angVel;

	long getTime();
	long deltaT(long, long);
public:
	Simulation();
	~Simulation();

	CBody* getShip();
	QVector<vector4> getPath();
	Geometry* geometry();

	void setUseBoardTime(bool);
	bool useBoardTime();

	void run();

public slots:
	void launch();
	void stop();
	void reset();
	void simulate();
	void keyPressed(int);
	void keyReleased(int);

	void setPgForce(double);
	void setRtForce(double);
	void setAngVel(double);

signals:
	void timeWarp(double);
	void autoStop();
};

#endif
