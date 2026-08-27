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
#include "tipMagneticForce.h"
#include "inertialForce.h"
#include "externalContactForce.h"
#include "contact/planarBarrierContactForce.h"
#include "geometry/planarSlabDomain.h"
#include "geometry/sphericalShellDomain.h"
#include "geometry/sphericalObstacleDomain.h"
#include "geometry/doubleSphericalObstacleDomain.h"
#include "contact/contactKktSystem.h"
#include "contact/contactKktEquilibriumResult.h"
#include "contact/contactConstraintAnalysis.h"
#include "continuation/fieldContinuation.h"
#include "continuation/actuationContinuation.h"
#include "actuation.h"
#include "plannerState.h"
#include "insertion/insertionModel.h"
#include "insertion/proximalGuideInsertionModel.h"

// include external force
#include "dampingForce.h"

// include time stepper
#include "timeStepper.h"

// include input file and option
#include "setInput.h"
#include "staticEvaluation.h"
#include "equilibriumResult.h"

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
	StaticEvaluation evaluateStaticSystem();
	EquilibriumResult solveStaticEquilibrium();
	void setAppliedField(const Vector3d &field);
	Vector3d getAppliedField() const;
	MatrixXd computeConfigurationFieldSensitivity();
	const ContactDetectionResult &getLastContactDetection() const;
	PlanarContactKktSeed buildPlanarContactKktSeed();
	ContactKktEquilibriumResult solvePlanarContactKktEquilibrium();
	ContactKktEquilibriumResult correctPlanarContactKktEquilibrium(
		const ContactKktEquilibriumResult &warmStart);
	ContactStabilityAnalysis analyzePlanarContactStability(
		const ContactKktEquilibriumResult &equilibrium);
	ContactEquilibriumSensitivity computePlanarContactFieldSensitivity(
		const ContactKktEquilibriumResult &equilibrium);
	ActuationEquilibriumSensitivity computePlanarContactActuationSensitivity(
		const ContactKktEquilibriumResult &equilibrium);
	Actuation getActuation() const;
	void setActuation(const Actuation &actuation);
	PlannerState capturePlannerState(
		const ContactKktEquilibriumResult &equilibrium,
		const ContactStabilityAnalysis &stability) const;
	void restorePlannerState(const PlannerState &state);
	FieldContinuationResult continuePlanarContactField(
		const Vector3d &targetField,
		const FieldContinuationOptions &options = FieldContinuationOptions());
	ActuationContinuationResult continuePlanarContactActuation(
		const Actuation &targetActuation,
		const FieldContinuationOptions &options = FieldContinuationOptions());
	ActuationContinuationResult continuePlanarContactActuation(
		const PlannerState &startState,
		const Actuation &targetActuation,
		const FieldContinuationOptions &options = FieldContinuationOptions());
	RodState captureRodState() const;
	void restoreRodState(const RodState &state);
	void applyFreeDofIncrement(const VectorXd &increment);
	int numFreeDofs() const;
	int simulationRunning();
	int numPoints();
	double getScaledCoordinate(int i);
	double getCurrentTime();
	double getTotalTime();
	double getVelocityNorm() const;
	
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
	string magneticModel;
	double tipDipoleMoment;

    double thickness;
	Vector3d shellCenter;
	double shellRadius;
	double shellMinusThickness;
	double shellPlusThickness;
	Vector3d cavityCenter;
	double cavityRadius;
	Vector3d obstacleCenter;
	double obstacleRadius;
	Vector3d secondObstacleCenter;
	double secondObstacleRadius;
    double dBar;
    double stiffness;
	string contactModel;
	double tipSafeDistance;
	int maxLineSearchIter;
	double lineSearchReduction;
	double lineSearchArmijo;
	double kktGapTolerance;
	double kktMultiplierTolerance;
	double kktComplementarityTolerance;
	int kktMaxActiveSetUpdates;
	string insertionModel;
	double insertionCoordinate;
	double insertionStiffness;
	Vector3d insertionAxis;
    
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
	tipMagneticForce *m_tipMagneticForce = nullptr;
	dampingForce *m_dampingForce = nullptr;
	externalContactForce *m_externalContactForce = nullptr;
	ConfinedDomain *m_contactDomain = nullptr;
	PlanarBarrierContactForce *m_planarBarrierContactForce = nullptr;
	InsertionModel *m_insertionModel = nullptr;
	const ContactKktEquilibriumResult *m_kktWarmStart = nullptr;
	
	int Nstep;
	int timeStep;
	int iter;

	void rodGeometry();
	void rodBoundaryCondition();
	ActuationContinuationResult continuePlanarContactActuationImpl(
		const Actuation &targetActuation,
		const FieldContinuationOptions &options,
		const ContactKktEquilibriumResult *warmStart);
	void assembleStaticSystem(bool includeContact = true);
	void assembleDynamicSystem();
	StaticEvaluation evaluatePhysicalStaticSystem();
	ContactKktLinearization evaluatePlanarKktLinearization(
		std::vector<ContactCandidate> &contacts,
		const VectorXd &multipliers);
	double computeResidualNorm() const;
    
	bool render; // should the OpenGL rendering be included?
	bool saveData; // should data be written to a file?

	Vector3d xInitial;
};

#endif
