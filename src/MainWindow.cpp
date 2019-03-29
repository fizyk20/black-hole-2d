#include "MainWindow.h"
#include <math.h>

#define POS_PREC 4
#define VEL_PREC 7
#define ROW_MIN_H 15

MainWindow::MainWindow()
{
	setWindowTitle("Black Hole");
	resize(1110,635);
	setFixedSize(1110,635);
	createComponents();
	createMenus();
	show();
	statusBar()->showMessage("Time warp: 1.0000x");
	timer = NULL;
}

MainWindow::~MainWindow()
{
	if(timer != NULL) delete timer;
	delete gl;
	delete sim;
}

void MainWindow::createComponents()
{
	sim = new Simulation;
	sim -> start();

	gl = new GraphicsWidget(this, sim);
	gl->move(10,10);
	gl->resize(800,600);
	gl->show();

	QWidget* components = new QWidget(this);
	components -> move(810,0);
	components -> resize(300,600);

	QGridLayout* main = new QGridLayout(components);
	components -> setLayout(main);

	mass = new QLineEdit("0.1");
	mass->setAlignment(Qt::AlignRight);
	mass->setValidator(new QDoubleValidator(mass));

	angmom = new QSlider(Qt::Horizontal);
	angmom -> setMinimum(0);
	angmom -> setMaximum(100);
	angmom -> setValue(0);
	angmom -> setTickPosition(QSlider::TicksBelow);
	angmom -> setTickInterval(10);

	l_angmom = new QLabel("Angular momentum / mass: 0.00");

	main -> addWidget(new QLabel("Black hole mass:"), 0, 0, 1, 2);
	main -> addWidget(mass, 0, 2, 1, 2);
	main -> addWidget(l_angmom, 1, 0, 1, 4);
	main -> addWidget(angmom, 2, 0, 1, 4);

	main -> setRowMinimumHeight(3, ROW_MIN_H);

	//position
	r = new QLineEdit("1.0000");
	r->setAlignment(Qt::AlignRight);
	r->setValidator(new QDoubleValidator(r));
	phi = new QLineEdit("0.0000");
	phi->setAlignment(Qt::AlignRight);
	phi->setValidator(new QDoubleValidator(phi));
	x = new QLineEdit("1.0000");
	x->setAlignment(Qt::AlignRight);
	x->setValidator(new QDoubleValidator(x));
	y = new QLineEdit("0.0000");
	y->setAlignment(Qt::AlignRight);
	y->setValidator(new QDoubleValidator(y));

	main->addWidget(new QLabel("R:"), 4, 0);
	main->addWidget(r, 4, 1);
	main->addWidget(new QLabel("Phi:"), 5, 0);
	main->addWidget(phi, 5, 1);
	main->addWidget(new QLabel("X:"), 4, 2);
	main->addWidget(x, 4, 3);
	main->addWidget(new QLabel("Y:"), 5, 2);
	main->addWidget(y, 5, 3);

	main -> setRowMinimumHeight(6, ROW_MIN_H);

	ur = new QLineEdit("0.0000000");
	ur->setAlignment(Qt::AlignRight);
	ur->setValidator(new QDoubleValidator(ur));
	uphi = new QLineEdit("0.3162278");
	uphi->setAlignment(Qt::AlignRight);
	uphi->setValidator(new QDoubleValidator(uphi));
	ux = new QLineEdit("0.0000000");
	ux->setAlignment(Qt::AlignRight);
	ux->setValidator(new QDoubleValidator(ux));
	uy = new QLineEdit("0.3162278");
	uy->setAlignment(Qt::AlignRight);
	uy->setValidator(new QDoubleValidator(uy));

	main->addWidget(new QLabel("VR:"), 7, 0);
	main->addWidget(ur, 7, 1);
	main->addWidget(new QLabel("VPhi:"), 8, 0);
	main->addWidget(uphi, 8, 1);
	main->addWidget(new QLabel("VX:"), 7, 2);
	main->addWidget(ux, 7, 3);
	main->addWidget(new QLabel("VY:"), 8, 2);
	main->addWidget(uy, 8, 3);

	main -> setRowMinimumHeight(9, ROW_MIN_H);

	vel = new QLabel("0.31662278");
	proper_vel = new QLabel("0.3535534");
	board_time = new QLabel("0.0000");
	ext_time = new QLabel("0.0000");

	main -> addWidget(new QLabel("Velocity:"), 10, 0, 1, 2);
	main -> addWidget(vel, 10, 2, 1, 2);
	main -> addWidget(new QLabel("Proper velocity:"), 11, 0, 1, 2);
	main -> addWidget(proper_vel, 11, 2, 1, 2);
	main -> addWidget(new QLabel("Time:"), 12, 0, 1, 2);
	main -> addWidget(ext_time, 12, 2, 1, 2);
	main -> addWidget(new QLabel("On-board time:"), 13, 0, 1, 2);
	main -> addWidget(board_time, 13, 2, 1, 2);

	main -> setRowMinimumHeight(14, ROW_MIN_H);

	useBoardTime = new QCheckBox("Use on-board time for simulation");
	main -> addWidget(useBoardTime, 15, 0, 1, 4);

	main -> setRowMinimumHeight(16, ROW_MIN_H);

	start = new QPushButton("Start", this);
	main -> addWidget(start, 17, 0, 1, 2);
	stop = new QPushButton("Stop", this);
	stop -> setEnabled(false);
	main -> addWidget(stop, 17, 2, 1, 2);
	reset = new QPushButton("Reset", this);
	main -> addWidget(reset, 18, 0, 1, 4);

	main -> setRowMinimumHeight(19, ROW_MIN_H);

	main -> addWidget(new QLabel("Thrust controls"), 20, 0, 1, 4);

	F_slider = new QSlider(Qt::Horizontal);
	F_slider -> setMinimum(-80);
	F_slider -> setMaximum(20);
	F_slider -> setValue(-20);
	F_slider -> setTickPosition(QSlider::TicksBelow);
	main -> addWidget(F_slider, 21, 0, 1, 4);

	force = new QLabel("Thrust: 0.1000");
	main -> addWidget(force, 22, 0, 1, 4);
	thrust_pg = new QLabel("Current forward thrust:  0.0000");
	thrust_rt = new QLabel("Current right thrust:	0.0000");
	main -> addWidget(thrust_pg, 23, 0, 1, 4);
	main -> addWidget(thrust_rt, 24, 0, 1, 4);

	main -> setRowStretch(25, 5);

	shipTracking = new QCheckBox("Ship tracking");
	main -> addWidget(shipTracking, 26, 0, 1, 4);

	//connect(sim, SIGNAL(update()), this, SLOT(stateUpdate()));
	//connect(sim, SIGNAL(update()), gl, SLOT(updateGL()));
	connect(this, SIGNAL(destroyed()), sim, SLOT(quit()));
	connect(sim, SIGNAL(autoStop()), this, SLOT(autoStop()));
	connect(sim, SIGNAL(timeWarp(double)), this, SLOT(timeWarp(double)));

	connect(this, SIGNAL(keyPressed(int)), sim, SLOT(keyPressed(int)));
	connect(this, SIGNAL(keyReleased(int)), sim, SLOT(keyReleased(int)));

	connect(start, SIGNAL(clicked()), sim, SLOT(launch()));
	connect(start, SIGNAL(clicked()), this, SLOT(buttonClicked()));
	connect(stop, SIGNAL(clicked()), sim, SLOT(stop()));
	connect(stop, SIGNAL(clicked()), this, SLOT(buttonClicked()));
	connect(reset, SIGNAL(clicked()), sim, SLOT(reset()));
	connect(reset, SIGNAL(clicked()), this, SLOT(buttonClicked()));

	connect(this, SIGNAL(enable(bool)), useBoardTime, SLOT(setEnabled(bool)));
	connect(this, SIGNAL(enable(bool)), mass, SLOT(setEnabled(bool)));
	connect(this, SIGNAL(enable(bool)), angmom, SLOT(setEnabled(bool)));
	connect(this, SIGNAL(enable(bool)), x, SLOT(setEnabled(bool)));
	connect(this, SIGNAL(enable(bool)), y, SLOT(setEnabled(bool)));
	connect(this, SIGNAL(enable(bool)), r, SLOT(setEnabled(bool)));
	connect(this, SIGNAL(enable(bool)), phi, SLOT(setEnabled(bool)));
	connect(this, SIGNAL(enable(bool)), ux, SLOT(setEnabled(bool)));
	connect(this, SIGNAL(enable(bool)), uy, SLOT(setEnabled(bool)));
	connect(this, SIGNAL(enable(bool)), ur, SLOT(setEnabled(bool)));
	connect(this, SIGNAL(enable(bool)), uphi, SLOT(setEnabled(bool)));

	connect(mass, SIGNAL(textEdited(const QString&)), this, SLOT(massChanged(const QString&)));
	connect(angmom, SIGNAL(valueChanged(int)), this, SLOT(angmomChanged(int)));
	connect(x, SIGNAL(textEdited(const QString&)), this, SLOT(xChanged(const QString&)));
	connect(y, SIGNAL(textEdited(const QString&)), this, SLOT(yChanged(const QString&)));
	connect(r, SIGNAL(textEdited(const QString&)), this, SLOT(rChanged(const QString&)));
	connect(phi, SIGNAL(textEdited(const QString&)), this, SLOT(phiChanged(const QString&)));
	connect(ux, SIGNAL(textEdited(const QString&)), this, SLOT(uxChanged(const QString&)));
	connect(uy, SIGNAL(textEdited(const QString&)), this, SLOT(uyChanged(const QString&)));
	connect(ur, SIGNAL(textEdited(const QString&)), this, SLOT(urChanged(const QString&)));
	connect(uphi, SIGNAL(textEdited(const QString&)), this, SLOT(uphiChanged(const QString&)));

	connect(useBoardTime, SIGNAL(stateChanged(int)), this, SLOT(boardTime(int)));
	connect(shipTracking, SIGNAL(stateChanged(int)), this, SLOT(trackShip(int)));
	connect(F_slider, SIGNAL(valueChanged(int)), this, SLOT(forceChanged(int)));
	connect(this, SIGNAL(forceNew(double)), sim, SLOT(setPgForce(double)));
	connect(this, SIGNAL(forceNew(double)), sim, SLOT(setRtForce(double)));
}

void MainWindow::createMenus()
{
}

void MainWindow::timeWarp(double t)
{
	statusBar()->showMessage("Time warp: "+QString::number(t,'f',4)+"x");
}

void MainWindow::massChanged(const QString& s)
{
	sim->geometry()->setM(s.toDouble());
	gl -> updateGL();
}

void MainWindow::angmomChanged(int val)
{
	double a_per_M = (double) val / 100.0;
	if(a_per_M == 1.0) {
		a_per_M = 0.999;
	}
	double M = sim->geometry()->getM();
	l_angmom->setText("Angular momentum / mass: " + QString::number(a_per_M,'f',2));
	sim->geometry()->setA(a_per_M * M);
	gl -> updateGL();
}

void MainWindow::rChanged(const QString& s)
{
	double rd = s.toDouble();
	double phid = phi->text().toDouble();
	CBody* ship = sim->getShip();

	ship->setVector(0,sim->geometry()->convertFrom(vector4(0.0,rd,M_PI/2,phid),COORD_SPHERICAL));
	stateUpdate(r);
}

void MainWindow::phiChanged(const QString& s)
{
	double phid = s.toDouble();
	double rd = r->text().toDouble();
	CBody* ship = sim->getShip();

	ship->setVector(0,sim->geometry()->convertFrom(vector4(0.0,rd,M_PI/2,phid),COORD_SPHERICAL));
	stateUpdate(phi);
}

void MainWindow::xChanged(const QString& s)
{
	double xd = s.toDouble();
	double yd = y->text().toDouble();
	double rd = sqrt(xd*xd+yd*yd);
	double phid = atan2(yd,xd);

	CBody* ship = sim->getShip();

	ship->setVector(0,sim->geometry()->convertFrom(vector4(0.0,rd,M_PI/2,phid),COORD_SPHERICAL));
	stateUpdate(x);
}

void MainWindow::yChanged(const QString& s)
{
	double yd = s.toDouble();
	double xd = x->text().toDouble();
	double rd = sqrt(xd*xd+yd*yd);
	double phid = atan2(yd,xd);

	CBody* ship = sim->getShip();

	ship->setVector(0,sim->geometry()->convertFrom(vector4(0.0,rd,M_PI/2,phid),COORD_SPHERICAL));
	stateUpdate(y);
}

void MainWindow::urChanged(const QString& s)
{
	CBody* ship = sim->getShip();
	vector4 x = ship->getVector(0);
	Tensor jac = sim->geometry()->Jacobian(COORD_SPHERICAL,x);
	Tensor invjac = sim->geometry()->invJacobian(COORD_SPHERICAL,x);
	vector4 vel = toVector4((invjac*Vector(ship->getVector(1))).contracted(1,2));
	vel[1] = s.toDouble();
	ship->setVector(1,toVector4((jac*Vector(vel)).contracted(1,2)));

	stateUpdate(ur);
}

void MainWindow::uphiChanged(const QString& s)
{
	CBody* ship = sim->getShip();
	vector4 x = ship->getVector(0);
	Tensor jac = sim->geometry()->Jacobian(COORD_SPHERICAL,x);
	Tensor invjac = sim->geometry()->invJacobian(COORD_SPHERICAL,x);
	vector4 vel = toVector4((invjac*Vector(ship->getVector(1))).contracted(1,2));
	vel[3] = s.toDouble();
	ship->setVector(1,toVector4((jac*Vector(vel)).contracted(1,2)));

	stateUpdate(uphi);
}

void MainWindow::uxChanged(const QString& s)
{
	CBody* ship = sim->getShip();
	vector4 x = ship->getVector(0);
	Tensor jac = sim->geometry()->Jacobian(COORD_CARTESIAN,x);
	Tensor invjac = sim->geometry()->invJacobian(COORD_CARTESIAN,x);
	vector4 vel = toVector4((invjac*Vector(ship->getVector(1))).contracted(1,2));
	vel[1] = s.toDouble();
	ship->setVector(1,toVector4((jac*Vector(vel)).contracted(1,2)));

	stateUpdate(ux);
}

void MainWindow::uyChanged(const QString& s)
{
	CBody* ship = sim->getShip();
	vector4 x = ship->getVector(0);
	Tensor jac = sim->geometry()->Jacobian(COORD_CARTESIAN,x);
	Tensor invjac = sim->geometry()->invJacobian(COORD_CARTESIAN,x);
	vector4 vel = toVector4((invjac*Vector(ship->getVector(1))).contracted(1,2));
	vel[2] = s.toDouble();
	ship->setVector(1,toVector4((jac*Vector(vel)).contracted(1,2)));

	stateUpdate(uy);
}

void MainWindow::buttonClicked()
{
	if(sender() == (QObject*) start)
	{
		emit enable(false);
		start -> setEnabled(false);
		stop -> setEnabled(true);
		timer = new QTimer(this);
		connect(timer, SIGNAL(timeout()), this, SLOT(stateUpdate()));
		connect(timer, SIGNAL(timeout()), gl, SLOT(updateGL()));
		timer->start(20);
	}
	if(sender() == (QObject*) reset)
	{
		gl -> setView(0.0, 0.0, 1.0);
	}
	if(sender() == (QObject*) stop || sender() == (QObject*) reset)
	{
		if(timer != NULL) delete timer;
		timer = NULL;
		emit enable(true);
		start -> setEnabled(true);
		stop -> setEnabled(false);
		stateUpdate();
		gl->updateGL();
	}
}

void MainWindow::stateUpdate(QObject* o)
{
	CBody* ship = sim->getShip();
	Geometry* g = sim->geometry();
	vector4 u = ship->getVector(1), xv = ship->getVector(0);
	Tensor jac_bl, invjac_bl;
	jac_bl = g->Jacobian(COORD_SPHERICAL,xv);
	invjac_bl = g->invJacobian(COORD_SPHERICAL,xv);
	//vector4 cart_u = sim->geometry()->convertVel(u,xv), cart_x = sim->geometry()->convertPos(xv);
	vector4 spher_x = g->convertTo(xv,COORD_SPHERICAL);
	vector4 spher_u = toVector4((invjac_bl*Vector(u)).contracted(1,2));
	vector4 cart_x = g->convertTo(xv,COORD_CARTESIAN);
	vector4 cart_u = toVector4((g->invJacobian(COORD_CARTESIAN,xv)*Vector(u)).contracted(1,2));

	ext_time -> setText(QString::number(cart_x[0],'f',POS_PREC));
	board_time -> setText(QString::number(ship->getTau(),'f',POS_PREC));

	if(o != (QObject*) x) x -> setText(QString::number(cart_x[1],'f',POS_PREC));
	if(o != (QObject*) y) y -> setText(QString::number(cart_x[2],'f',POS_PREC));
	if(o != (QObject*) r) r -> setText(QString::number(spher_x[1],'f',POS_PREC));
	if(o != (QObject*) phi) phi -> setText(QString::number(spher_x[3],'f',POS_PREC));

	double vx = cart_u[1]/cart_u[0], vy = cart_u[2]/cart_u[0];
	double vr = spher_u[1]/spher_u[0], vphi = spher_u[3]/spher_u[0];
	if(o != (QObject*) ux) ux -> setText(QString::number(vx,'f',VEL_PREC));
	if(o != (QObject*) uy) uy -> setText(QString::number(vy,'f',VEL_PREC));
	if(o != (QObject*) ur) ur -> setText(QString::number(vr,'f',VEL_PREC));
	if(o != (QObject*) uphi) uphi -> setText(QString::number(vphi,'f',VEL_PREC));
	vel -> setText(QString::number(sqrt(vx*vx+vy*vy),'f',VEL_PREC));

	//Calculating proper velocity
	Tensor g_bl = (g->g(xv)*jac_bl*jac_bl).contracted(0,2).contracted(0,2);
	Tensor tv = Vector(vector4(1,0,0,0));
	double ltv = toDouble((g_bl*tv*tv).contracted(0,2).contracted(0,1));
	vector4 u_perp, u_para;
	u_para = toDouble((g_bl*tv*Vector(spher_u)).contracted(0,2).contracted(0,1))*vector4(1,0,0,0)/ltv;
	u_perp = spher_u - u_para;
	double vel_perp = toDouble((g_bl*Vector(u_perp)*Vector(u_perp)).contracted(0,2).contracted(0,1));
	double vel_para = toDouble((g_bl*Vector(u_para)*Vector(u_para)).contracted(0,2).contracted(0,1));
	proper_vel -> setText(QString::number(sqrt(fabs(vel_perp/vel_para)),'f',VEL_PREC));

	thrust_pg -> setText("Current forward thrust:  "+QString::number(ship->getForce()[1],'f',4));
	thrust_rt -> setText("Current right thrust:	"+QString::number(ship->getForce()[0],'f',4));

	//DEBUG
	//vector3 f = ship->getCartesianDir(ship->getVector(3));
	//statusBar()->showMessage(QString::number(atan2(f[1],f[0]),'f',3)+"   "+QString::number(sim->geometry()->g_v(ship->getVector(1),ship->getVector(1),ship->getVector(0)),'f',4));
}

void MainWindow::boardTime(int state)
{
	if(state == Qt::Checked) sim->setUseBoardTime(true);
	if(state == Qt::Unchecked) sim->setUseBoardTime(false);
}

void MainWindow::trackShip(int state)
{
	if(state == Qt::Checked) gl->setTrackShip(true);
	if(state == Qt::Unchecked) gl->setTrackShip(false);
	gl->updateGL();
}

void MainWindow::forceChanged(int val)
{
	double f = pow(10, (double)val/20.0);
	force -> setText("Thrust: "+QString::number(f,'f',4));
	emit forceNew(f);
}

void MainWindow::autoStop()
{
	if(timer != NULL) delete timer;
	timer = NULL;
	emit enable(true);
	QMessageBox::information(this,"Automatic stop","Simulation was automatically stopped because you reached the singularity. Nobody knows what happens to you now.");
}

void MainWindow::keyPressEvent(QKeyEvent* event)
{
	emit keyPressed(event->key());
}

void MainWindow::keyReleaseEvent(QKeyEvent* event)
{
	emit keyReleased(event->key());
}
