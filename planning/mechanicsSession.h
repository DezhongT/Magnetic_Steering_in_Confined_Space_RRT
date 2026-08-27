#ifndef MECHANICSSESSION_H
#define MECHANICSSESSION_H

#include "actuation.h"
#include "continuation/fieldContinuation.h"
#include "plannerState.h"

#include <limits>
#include <memory>
#include <string>

class world;

struct MechanicsConfig
{
	double rodLength = 1.0;
	double rodRadius = 1.0e-2;
	int numVertices = 20;
	double youngModulus = 1.0e7;
	double poissonRatio = 0.5;
	int maximumNewtonIterations = 200;
	Vector3d gravity = Vector3d::Zero();
	std::string domainType = "planar_slab";
	double planeHalfThickness = 1.0e-1;
	Vector3d shellCenter = Vector3d::Zero();
	double shellRadius = 1.0;
	double shellMinusThickness = 0.2;
	double shellPlusThickness = 0.2;
	Vector3d cavityCenter = Vector3d::Zero();
	double cavityRadius = 2.0;
	Vector3d obstacleCenter = Vector3d(0.0, 0.5, 0.0);
	double obstacleRadius = 0.1;
	Vector3d secondObstacleCenter = Vector3d(0.0, -0.5, 0.0);
	double secondObstacleRadius = 0.1;
	double barrierDistance = 1.5e-2;
	double barrierStiffness = 1.0e2;
	double tipSafeDistance = 0.0;
	double tipDipoleMoment = 1.0e-3;
	double insertionStiffness = 1.0e3;
	Vector3d insertionAxis = Vector3d::UnitX();
	Actuation initialActuation;
};

struct ContinuationEdgeResult
{
	bool success = false;
	bool rolledBack = false;
	PlannerState state;
	int attemptedSteps = 0;
	int rejectedSteps = 0;
	int storedPoints = 0;
	int contactsAdded = 0;
	int contactsReleased = 0;
	double reachedPathFraction = 0.0;
	double minimumStabilityMargin =
		std::numeric_limits<double>::infinity();
	std::string failureReason;
};

struct LocalSteeringResult
{
	Vector3d tipPosition = Vector3d::Zero();
	Matrix<double, 3, 4> tipActuationJacobian =
		Matrix<double, 3, 4>::Zero();
	double linearResidualNorm = std::numeric_limits<double>::infinity();
};

class MechanicsSession
{
public:
	explicit MechanicsSession(const MechanicsConfig &config = MechanicsConfig());
	~MechanicsSession();
	MechanicsSession(const MechanicsSession &) = delete;
	MechanicsSession &operator=(const MechanicsSession &) = delete;

	PlannerState solveInitialState();
	LocalSteeringResult evaluateLocalSteering(const PlannerState &state);
	ContinuationEdgeResult attemptContinuation(
		const PlannerState &startState,
		const Actuation &targetActuation,
		const FieldContinuationOptions &options = FieldContinuationOptions());

private:
	std::unique_ptr<world> simulation;
};

#endif
