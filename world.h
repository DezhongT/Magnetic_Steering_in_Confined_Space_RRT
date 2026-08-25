#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

#include <time.h>

// include elastic rod class
#include "elasticRod.h"

// include force classes
#include "elasticStretchingForce.h"
#include "elasticBendingForce.h"
#include "elasticTwistingForce.h"
#include "externalGravityForce.h"
#include "externalMagneticForce.h"
#include "inertialForce.h"
#include "externalContactForce.h"

// include external force
#include "dampingForce.h"

// include time stepper
#include "timeStepper.h"

// include input file and option
#include "setInput.h"

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	world(const world&) = delete;
	world& operator=(const world&) = delete;
	void setRodStepper();
	void updateTimeStep();
	double getStaticResidualNorm();
	int simulationRunning();
	int numPoints();
	double getScaledCoordinate(int i);
	double getCurrentTime();
	double getTotalTime();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);

	Vector3d getScaledCoordinateSurface(int i, int j);
	int numTriangle();
		
private:

	// Physical parameters
	double RodLength;
	double rodRadius;
	int numVertices;
	double youngM;
	double Poisson;
	double shearM;
	double deltaTime;
	double totalTime;
	double density;
	Vector3d gVector;
	double viscosity;
	double scaleRendering;

	Vector3d baVector;
    Vector3d brVector;
    double muZero;

    double thickness;
    double dBar;
    double stiffness;
    
	double tol, stol;
	int maxIter; // maximum number of iterations
	double characteristicForce;
	double forceTol;
	
	// Geometry
	MatrixXd vertices;
	
	// Rod
	elasticRod *rod = nullptr;
	
	// set up the time stepper
	timeStepper *stepper = nullptr;
	double *totalForce = nullptr;
	double currentTime;
	
	// declare the forces
	elasticStretchingForce *m_stretchForce = nullptr;
	elasticBendingForce *m_bendingForce = nullptr;
	elasticTwistingForce *m_twistingForce = nullptr;
	inertialForce *m_inertialForce = nullptr;
	externalGravityForce *m_gravityForce = nullptr;
	externalMagneticForce *m_magneticForce = nullptr;
	dampingForce *m_dampingForce = nullptr;
	externalContactForce *m_externalContactForce = nullptr;
	
	int Nstep;
	int timeStep;
	int iter;

	void rodGeometry();
	void rodBoundaryCondition();
	void assembleStaticSystem();
	void assembleDynamicSystem();
	double computeResidualNorm() const;
    
	bool render; // should the OpenGL rendering be included?
	bool saveData; // should data be written to a file?

	Vector3d xInitial;
};

#endif
