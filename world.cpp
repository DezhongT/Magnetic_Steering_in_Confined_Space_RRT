#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				// boolean
	saveData = m_inputData.GetBoolOpt("saveData");			// boolean
	
	// Physical parameters
	RodLength = m_inputData.GetScalarOpt("RodLength");      // meter
    rodRadius = m_inputData.GetScalarOpt("rodRadius");      // meter

    gVector = m_inputData.GetVecOpt("gVector");             // m/s^2
    maxIter = m_inputData.GetIntOpt("maxIter");             // maximum number of iterations
	numVertices = m_inputData.GetIntOpt("numVertices");     // int_num
	youngM = m_inputData.GetScalarOpt("youngM");            // Pa
	Poisson = m_inputData.GetScalarOpt("Poisson");          // dimensionless
	deltaTime = m_inputData.GetScalarOpt("deltaTime");      // seconds
	totalTime= m_inputData.GetScalarOpt("totalTime");       // seconds
	tol = m_inputData.GetScalarOpt("tol");                  // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");				// small number, e.g. 0.1%
	density = m_inputData.GetScalarOpt("density");          // kg/m^3
	viscosity = m_inputData.GetScalarOpt("viscosity");      // viscosity in Pa-s

	baVector = m_inputData.GetVecOpt("baVector");           // magnetic field
	brVector = m_inputData.GetVecOpt("brVector");           // magnetic field
	muZero = m_inputData.GetScalarOpt("muZero");            // magnetic field

	scaleRendering = m_inputData.GetScalarOpt("scaleRendering");

	thickness =m_inputData.GetScalarOpt("thickness");
	dBar =m_inputData.GetScalarOpt("dBar");
	stiffness =m_inputData.GetScalarOpt("stiffness");

	
	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus
}

world::~world()
{
	delete m_externalContactForce;
	delete m_dampingForce;
	delete m_magneticForce;
	delete m_gravityForce;
	delete m_inertialForce;
	delete m_twistingForce;
	delete m_bendingForce;
	delete m_stretchForce;
	delete stepper;
	delete rod;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
	
	int systemRet = system("mkdir datafiles"); //make the directory

	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	ostringstream name;
	name.precision(6);
	name << fixed;
    name << "datafiles/simDER";
    //name << "_rodLength_" << RodLength;
    //name << "_Ba_" << baVector(0);
    name << ".txt";

    outfile.open(name.str().c_str());
    outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	outfile.close();
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	if (timeStep == Nstep)
	{
		;
	}
}

void world::setRodStepper()
{
	// Set up geometry
	rodGeometry();	

	// Create the rod 
	rod = new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
		youngM, shearM, RodLength);

	// Find out the tolerance, e.g. how small is enough?
	characteristicForce = M_PI * pow(rodRadius ,4)/4.0 * youngM / pow(RodLength, 2);
	forceTol = tol * characteristicForce;
	
	// Set up boundary condition
	rodBoundaryCondition();
	
	// setup the rod so that all the relevant variables are populated
	rod->setup();
	// End of rod setup
	
	// set up the time stepper
	stepper = new timeStepper(*rod);
	totalForce = stepper->getForce();

	// declare the forces
	m_stretchForce = new elasticStretchingForce(*rod, *stepper);
	m_bendingForce = new elasticBendingForce(*rod, *stepper);
	m_twistingForce = new elasticTwistingForce(*rod, *stepper);
	m_inertialForce = new inertialForce(*rod, *stepper);
	m_gravityForce = new externalGravityForce(*rod, *stepper, gVector);
	m_magneticForce = new externalMagneticForce(*rod, *stepper, baVector, brVector, muZero);
	m_dampingForce = new dampingForce(*rod, *stepper, viscosity);
	m_externalContactForce = new externalContactForce(*rod, *stepper, dBar, stiffness, thickness);

	Nstep = totalTime/deltaTime;

	// Allocate every thing to prepare for the first iteration
	rod->updateTimeStep();
	
	timeStep = 0;
	currentTime = 0.0;

	xInitial = rod->getVertex(rod->nv - 1);
}

// Setup geometry
void world::rodGeometry()
{
	/*
	vertices = MatrixXd(numVertices, 3);

    double helixA = helixradius;
    double helixB = helixpitch / (2.0 * M_PI); 

    double helixT = RodLength / sqrt(helixA*helixA + helixB*helixB);
    double delta_t = helixT / (numVertices - 1); // step for t->[0, T]

    int i = 0;
    for (double tt = 0.0; i < numVertices; tt += delta_t)
    {
        vertices(i, 0) = helixB * tt;
        vertices(i, 1) = helixA * cos(tt);
        vertices(i, 2) = helixA * sin(tt);
        i++;
    }

    double midPoint = vertices(numVertices/2, 0);

    for (int i = 0; i < numVertices; i++)
    {
    	vertices(i, 0) = vertices(i, 0) - midPoint;
    }
    */

    vertices = MatrixXd(numVertices, 3);

    double deltaL = RodLength / (numVertices - 1);

    for (int i = 0; i < numVertices; i++)
    {
    	vertices(i, 0) = deltaL * i - RodLength/2;
        vertices(i, 1) = 0.0;
        vertices(i, 2) = 0.0;
    }

}

void world::rodBoundaryCondition()
{
	// Apply boundary condition

	rod->setThetaBoundaryCondition(0, 0);

	Vector3d x0 = rod->getVertex(0);
	rod->setVertexBoundaryCondition(x0, 0);

	Vector3d x1 = rod->getVertex(1);
	rod->setVertexBoundaryCondition(x1, 1);
}

void world::assembleStaticSystem()
{
	stepper->setZero();

	// Internal elastic mechanics.
	m_stretchForce->computeFs();
	m_stretchForce->computeJs();

	m_bendingForce->computeFb();
	m_bendingForce->computeJb();

	m_twistingForce->computeFt();
	m_twistingForce->computeJt();

	// Configuration-dependent and static external loading.
	m_gravityForce->computeFg();
	m_gravityForce->computeJg();

	m_magneticForce->computeFm();
	m_magneticForce->computeJm();

	m_externalContactForce->computeFc();
}

void world::assembleDynamicSystem()
{
	assembleStaticSystem();

	// Dynamic relaxation terms are deliberately kept separate from the
	// physical static residual assembled above.
	m_inertialForce->computeFi();
	m_inertialForce->computeJi();

	m_dampingForce->computeFd();
	m_dampingForce->computeJd();
}

double world::computeResidualNorm() const
{
	double squaredNorm = 0.0;
	for (int i = 0; i < rod->uncons; ++i)
	{
		squaredNorm += totalForce[i] * totalForce[i];
	}
	return sqrt(squaredNorm);
}

double world::getStaticResidualNorm()
{
	rod->prepareForIteration();
	assembleStaticSystem();
	return computeResidualNorm();
}
	

void world::updateTimeStep()
{
	double normf = forceTol * 10.0;	
	double normf0 = 0;
	
	bool solved = false;
	
	iter = 0;

	// Start with a trial solution for our solution x
	rod->updateGuess(); // x = x0 + u * dt
		
	while (solved == false)
	{
		rod->prepareForIteration();
		assembleDynamicSystem();
		
		// Compute norm of the force equations.
		normf = computeResidualNorm();
		if (iter == 0) normf0 = normf;
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}
		
		if (solved == false)
		{
			const int linearSolveInfo = stepper->integrator(); // Solve equations of motion
			if (linearSolveInfo != 0)
			{
				cerr << "Error. LAPACK band solve failed with info="
					 << linearSolveInfo << ". Exiting.\n";
				break;
			}
			rod->updateNewtonX(totalForce); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}
	
	rod->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << " iter=" << iter << endl;
	}

	currentTime += deltaTime;
		
	timeStep++;
	
	if (solved == false)
	{
		timeStep = Nstep; // we are exiting
	}
}

int world::simulationRunning()
{
	if (timeStep < Nstep) 
	{
		return 1;
	}
	else 
	{
		return -1;
	}
}

int world::numPoints()
{
	return rod->nv;
}

double world::getScaledCoordinate(int i)
{
	return rod->x[i] * scaleRendering;
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}

Vector3d world::getScaledCoordinateSurface(int i, int j)
{
	
	Vector3i NodeIndex = rod->v_triangular[i];

	Vector3d xCurrent = rod->v_nodes[NodeIndex(j)];

	return xCurrent * scaleRendering;
}

int world::numTriangle()
{
	return rod->plateTri;
}
