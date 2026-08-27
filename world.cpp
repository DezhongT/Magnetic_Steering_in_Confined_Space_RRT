#include "world.h"

#include <stdexcept>

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
	magneticModel = m_inputData.GetStringOpt("magneticModel");
	tipDipoleMoment = m_inputData.GetScalarOpt("tipDipoleMoment");

	scaleRendering = m_inputData.GetScalarOpt("scaleRendering");

	thickness =m_inputData.GetScalarOpt("thickness");
	shellCenter = m_inputData.GetVecOpt("shellCenter");
	shellRadius = m_inputData.GetScalarOpt("shellRadius");
	shellMinusThickness = m_inputData.GetScalarOpt("shellMinusThickness");
	shellPlusThickness = m_inputData.GetScalarOpt("shellPlusThickness");
	cavityCenter = m_inputData.GetVecOpt("cavityCenter");
	cavityRadius = m_inputData.GetScalarOpt("cavityRadius");
	obstacleCenter = m_inputData.GetVecOpt("obstacleCenter");
	obstacleRadius = m_inputData.GetScalarOpt("obstacleRadius");
	secondObstacleCenter = m_inputData.GetVecOpt("secondObstacleCenter");
	secondObstacleRadius = m_inputData.GetScalarOpt("secondObstacleRadius");
	dBar =m_inputData.GetScalarOpt("dBar");
	stiffness =m_inputData.GetScalarOpt("stiffness");
	contactModel = m_inputData.GetStringOpt("contactModel");
	tipSafeDistance = m_inputData.GetScalarOpt("tipSafeDistance");
	maxLineSearchIter = m_inputData.GetIntOpt("maxLineSearchIter");
	lineSearchReduction = m_inputData.GetScalarOpt("lineSearchReduction");
	lineSearchArmijo = m_inputData.GetScalarOpt("lineSearchArmijo");
	kktGapTolerance = m_inputData.GetScalarOpt("kktGapTolerance");
	kktMultiplierTolerance = m_inputData.GetScalarOpt("kktMultiplierTolerance");
	kktComplementarityTolerance =
		m_inputData.GetScalarOpt("kktComplementarityTolerance");
	kktMaxActiveSetUpdates = m_inputData.GetIntOpt("kktMaxActiveSetUpdates");
	insertionModel = m_inputData.GetStringOpt("insertionModel");
	insertionCoordinate = m_inputData.GetScalarOpt("insertionCoordinate");
	insertionStiffness = m_inputData.GetScalarOpt("insertionStiffness");
	insertionAxis = m_inputData.GetVecOpt("insertionAxis");
	if (maxLineSearchIter < 0 || !(lineSearchReduction > 0.0) ||
		!(lineSearchReduction < 1.0) || !(lineSearchArmijo > 0.0) ||
		!(lineSearchArmijo < 0.5))
	{
		throw invalid_argument("Invalid static Newton line-search parameters");
	}
	if (!(kktGapTolerance > 0.0) || !(kktMultiplierTolerance > 0.0) ||
		!(kktComplementarityTolerance > 0.0) || kktMaxActiveSetUpdates < 0)
	{
		throw invalid_argument("Invalid planar KKT solver parameters");
	}

	
	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus
}

world::~world()
{
	delete m_insertionModel;
	delete m_planarBarrierContactForce;
	delete m_contactDomain;
	delete m_externalContactForce;
	delete m_dampingForce;
	delete m_tipMagneticForce;
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
	(void)outfile;
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
	if (magneticModel == "legacy")
	{
		m_magneticForce = new externalMagneticForce(
			*rod, *stepper, baVector, brVector, muZero);
	}
	else if (magneticModel == "axial_tip")
	{
		m_tipMagneticForce = new tipMagneticForce(
			*rod, *stepper, baVector, tipDipoleMoment);
	}
	else
	{
		throw invalid_argument("Unknown magneticModel: " + magneticModel);
	}
	m_dampingForce = new dampingForce(*rod, *stepper, viscosity);
	if (insertionModel == "proximal_guide")
	{
		m_insertionModel = new ProximalGuideInsertionModel(
			*rod, *stepper, insertionAxis, insertionStiffness,
			insertionCoordinate);
	}
	else if (insertionModel != "none")
	{
		throw invalid_argument("Unknown insertionModel: " + insertionModel);
	}
	if (contactModel == "legacy")
	{
		m_externalContactForce = new externalContactForce(
			*rod, *stepper, dBar, stiffness, thickness);
	}
	else if (contactModel == "planar_barrier")
	{
		m_contactDomain = new PlanarSlabDomain(
			Vector3d::Zero(), Vector3d::UnitZ(), thickness, thickness);
		m_planarBarrierContactForce = new PlanarBarrierContactForce(
			*rod, *stepper, *m_contactDomain, rodRadius, dBar,
			stiffness, tipSafeDistance);
	}
	else if (contactModel == "spherical_shell_barrier")
	{
		m_contactDomain = new SphericalShellDomain(
			shellCenter, shellRadius, shellMinusThickness,
			shellPlusThickness);
		m_planarBarrierContactForce = new PlanarBarrierContactForce(
			*rod, *stepper, *m_contactDomain, rodRadius, dBar,
			stiffness, tipSafeDistance);
	}
	else if (contactModel == "spherical_obstacle_barrier")
	{
		m_contactDomain = new SphericalObstacleDomain(
			cavityCenter, cavityRadius, obstacleCenter, obstacleRadius);
		m_planarBarrierContactForce = new PlanarBarrierContactForce(
			*rod, *stepper, *m_contactDomain, rodRadius, dBar,
			stiffness, tipSafeDistance);
	}
	else if (contactModel == "double_spherical_obstacle_barrier")
	{
		m_contactDomain = new DoubleSphericalObstacleDomain(
			cavityCenter, cavityRadius, obstacleCenter, obstacleRadius,
			secondObstacleCenter, secondObstacleRadius);
		m_planarBarrierContactForce = new PlanarBarrierContactForce(
			*rod, *stepper, *m_contactDomain, rodRadius, dBar,
			stiffness, tipSafeDistance);
	}
	else if (contactModel != "none")
	{
		throw invalid_argument("Unknown contactModel: " + contactModel);
	}

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

void world::assembleStaticSystem(bool includeContact)
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
	if (m_insertionModel != nullptr)
	{
		m_insertionModel->assemble();
	}

	if (m_magneticForce != nullptr)
	{
		m_magneticForce->computeFm();
		m_magneticForce->computeJm();
	}
	if (m_tipMagneticForce != nullptr)
	{
		m_tipMagneticForce->computeFm();
		m_tipMagneticForce->computeJm();
	}

	if (includeContact && m_externalContactForce != nullptr)
	{
		m_externalContactForce->computeFc();
	}
	if (includeContact && m_planarBarrierContactForce != nullptr)
	{
		m_planarBarrierContactForce->computeFc();
	}
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
	return evaluateStaticSystem().residualNorm();
}

StaticEvaluation world::evaluateStaticSystem()
{
	rod->prepareForIteration();
	assembleStaticSystem();
	return stepper->captureEvaluation();
}

StaticEvaluation world::evaluatePhysicalStaticSystem()
{
	rod->prepareForIteration();
	assembleStaticSystem(false);
	return stepper->captureEvaluation();
}

ContactKktLinearization world::evaluatePlanarKktLinearization(
	std::vector<ContactCandidate> &contacts,
	const VectorXd &multipliers)
{
	if (m_contactDomain == nullptr)
	{
		throw logic_error("Contact KKT linearization requires a confined domain");
	}
	if (multipliers.size() != static_cast<int>(contacts.size()))
	{
		throw invalid_argument("Active contacts and KKT multipliers do not match");
	}
	ContactKktLinearization linearization;
	const StaticEvaluation physical = evaluatePhysicalStaticSystem();
	linearization.hessian = physical.denseJacobian();
	linearization.residual = physical.residual;
	const int count = static_cast<int>(contacts.size());
	linearization.constraintJacobian = MatrixXd::Zero(count, rod->uncons);
	linearization.gaps = VectorXd(count);
	linearization.multipliers = multipliers;
	for (int contact = 0; contact < count; ++contact)
	{
		ContactCandidate &candidate = contacts[contact];
		const Vector3d center = rod->getVertex(candidate.rodVertex);
		const DomainQuery query = m_contactDomain->queryBoundary(
			center, candidate.boundaryId);
		candidate.rodCenter = center;
		candidate.centerlineClearance = query.clearance;
		candidate.gap = query.clearance - rodRadius;
		candidate.closestBoundaryPoint = query.closestPoint;
		candidate.normal = query.clearanceNormal;
		candidate.gapHessian = query.clearanceHessian;
		linearization.gaps[contact] = candidate.gap;
		for (int component = 0; component < 3; ++component)
		{
			const int fullDof = 4 * candidate.rodVertex + component;
			if (rod->getIfConstrained(fullDof) == 0)
			{
				linearization.constraintJacobian(
					contact, rod->fullToUnconsMap[fullDof]) =
					candidate.normal[component];
			}
		}
		for (int row = 0; row < 3; ++row)
		{
			const int rowDof = 4 * candidate.rodVertex + row;
			if (rod->getIfConstrained(rowDof) != 0)
			{
				continue;
			}
			const int freeRow = rod->fullToUnconsMap[rowDof];
			for (int column = 0; column < 3; ++column)
			{
				const int columnDof = 4 * candidate.rodVertex + column;
				if (rod->getIfConstrained(columnDof) == 0)
				{
					linearization.hessian(
						freeRow, rod->fullToUnconsMap[columnDof]) -=
						multipliers[contact] * query.clearanceHessian(row, column);
				}
			}
		}
	}
	return linearization;
}

EquilibriumResult world::solveStaticEquilibrium()
{
	const RodState inputState = rod->captureState();
	EquilibriumResult result;
	result.state = inputState;

	for (int iteration = 0; iteration <= maxIter; ++iteration)
	{
		StaticEvaluation evaluation;
		try
		{
			evaluation = evaluateStaticSystem();
		}
		catch (const InfeasibleContactError &)
		{
			result.lineSearchFailed = true;
			break;
		}
		const double residualNorm = evaluation.residualNorm();
		if (iteration == 0)
		{
			result.initialResidualNorm = residualNorm;
		}
		result.finalResidualNorm = residualNorm;
		result.newtonIterations = iteration;

		if (!std::isfinite(residualNorm))
		{
			break;
		}

		if (residualNorm <= forceTol)
		{
			rod->commitStaticState();
			result.finalResidualNorm = evaluateStaticSystem().residualNorm();
			result.success = std::isfinite(result.finalResidualNorm) &&
				result.finalResidualNorm <= forceTol;
			if (result.success)
			{
				result.state = rod->captureState();
				return result;
			}
			break;
		}

		if (iteration == maxIter)
		{
			break;
		}

		VectorXd correction;
		result.linearSolverInfo = stepper->solveBandedSystem(evaluation, correction);
		if (result.linearSolverInfo != 0 || !correction.allFinite())
		{
			break;
		}
		if (m_planarBarrierContactForce == nullptr)
		{
			// Preserve the validated legacy/contact-free Newton trajectory. The
			// feasibility globalization below is specific to the new barrier mode.
			rod->applyFreeDofIncrement(-correction);
			continue;
		}

		const VectorXd direction = -correction;
		const bool currentContactActive =
			!m_planarBarrierContactForce->lastDetection().bodyCandidates.empty();
		const double merit = 0.5 * evaluation.residual.squaredNorm();
		const double meritSlope = evaluation.residual.dot(
			evaluation.multiplyJacobian(direction));
		if (!std::isfinite(meritSlope) || meritSlope >= 0.0)
		{
			result.lineSearchFailed = true;
			break;
		}

		const RodState iterationState = rod->captureState();
		double stepLength = 1.0;
		bool accepted = false;
		for (int trial = 0; trial <= maxLineSearchIter; ++trial)
		{
			rod->restoreState(iterationState);
			rod->applyFreeDofIncrement(stepLength * direction);

			try
			{
				const StaticEvaluation trialEvaluation = evaluateStaticSystem();
				const double trialMerit =
					0.5 * trialEvaluation.residual.squaredNorm();
				const double sufficientDecrease =
					merit + lineSearchArmijo * stepLength * meritSlope;
				const bool trialContactActive =
					!m_planarBarrierContactForce->lastDetection().bodyCandidates.empty();
				const bool contactFreeTransition =
					!currentContactActive && !trialContactActive;
				if (std::isfinite(trialMerit) &&
					(contactFreeTransition || trialMerit <= sufficientDecrease))
				{
					accepted = true;
					result.finalResidualNorm = trialEvaluation.residualNorm();
					result.minimumAcceptedStepLength = std::min(
						result.minimumAcceptedStepLength, stepLength);
					break;
				}
			}
			catch (const InfeasibleContactError &)
			{
				// A barrier trial outside the open feasible set is rejected.
				++result.infeasibleTrialRejections;
			}

			if (trial < maxLineSearchIter)
			{
				stepLength *= lineSearchReduction;
				++result.lineSearchBacktracks;
			}
		}

		if (!accepted)
		{
			rod->restoreState(iterationState);
			result.lineSearchFailed = true;
			break;
		}
	}

	rod->restoreState(inputState);
	result.state = inputState;
	return result;
}

const ContactDetectionResult &world::getLastContactDetection() const
{
	if (m_planarBarrierContactForce == nullptr)
	{
		throw logic_error(
			"Contact detection diagnostics require contactModel=planar_barrier");
	}
	return m_planarBarrierContactForce->lastDetection();
}

PlanarContactKktSeed world::buildPlanarContactKktSeed()
{
	if (m_planarBarrierContactForce == nullptr)
	{
		throw logic_error(
			"Planar KKT seeding requires contactModel=planar_barrier");
	}

	// Refresh candidates and barrier multipliers at the current configuration.
	(void)evaluateStaticSystem();
	const std::vector<ContactCandidate> detected =
		m_planarBarrierContactForce->lastDetection().bodyCandidates;

	std::vector<ContactCandidate> activeContacts;
	std::vector<VectorXd> constraintRows;
	std::vector<double> gaps;
	std::vector<double> multipliers;
	for (const ContactCandidate &candidate : detected)
	{
		const BarrierPotentialEvaluation barrier = evaluateBarrierPotential(
			candidate.gap, dBar, stiffness);
		if (!barrier.valid)
		{
			throw InfeasibleContactError(
				"Cannot seed a KKT system from a nonpositive contact gap");
		}
		if (!barrier.active)
		{
			continue;
		}

		VectorXd row = VectorXd::Zero(rod->uncons);
		for (int component = 0; component < 3; ++component)
		{
			const int fullDof = 4 * candidate.rodVertex + component;
			if (rod->getIfConstrained(fullDof) == 0)
			{
				row[rod->fullToUnconsMap[fullDof]] = candidate.normal[component];
			}
		}
		if (row.norm() == 0.0)
		{
			continue;
		}
		activeContacts.push_back(candidate);
		constraintRows.push_back(row);
		gaps.push_back(candidate.gap);
		multipliers.push_back(-barrier.firstDerivative);
	}
	if (activeContacts.empty())
	{
		throw runtime_error("No active planar barrier contacts are available for KKT seeding");
	}

	// Remove the barrier residual and Hessian. They provide the warm start and
	// initial compression estimate, but are not part of exact KKT stationarity.
	const StaticEvaluation physicalEvaluation = evaluatePhysicalStaticSystem();

	PlanarContactKktSeed seed;
	seed.contacts = activeContacts;
	seed.linearization.hessian = physicalEvaluation.denseJacobian();
	seed.linearization.residual = physicalEvaluation.residual;
	const int contactCount = static_cast<int>(activeContacts.size());
	seed.linearization.constraintJacobian =
		MatrixXd::Zero(contactCount, rod->uncons);
	seed.linearization.gaps = VectorXd(contactCount);
	seed.linearization.multipliers = VectorXd(contactCount);
	for (int contact = 0; contact < contactCount; ++contact)
	{
		seed.linearization.constraintJacobian.row(contact) =
			constraintRows[contact].transpose();
		seed.linearization.gaps[contact] = gaps[contact];
		seed.linearization.multipliers[contact] = multipliers[contact];
		const ContactCandidate &candidate = activeContacts[contact];
		for (int row = 0; row < 3; ++row)
		{
			const int rowDof = 4 * candidate.rodVertex + row;
			if (rod->getIfConstrained(rowDof) != 0)
			{
				continue;
			}
			const int freeRow = rod->fullToUnconsMap[rowDof];
			for (int column = 0; column < 3; ++column)
			{
				const int columnDof = 4 * candidate.rodVertex + column;
				if (rod->getIfConstrained(columnDof) == 0)
				{
					seed.linearization.hessian(
						freeRow, rod->fullToUnconsMap[columnDof]) -=
						multipliers[contact] *
						candidate.gapHessian(row, column);
				}
			}
		}
	}
	return seed;
}

ContactKktEquilibriumResult world::solvePlanarContactKktEquilibrium()
{
	if (m_planarBarrierContactForce == nullptr || m_contactDomain == nullptr)
	{
		throw logic_error(
			"Planar KKT equilibrium requires contactModel=planar_barrier");
	}

	const RodState inputState = rod->captureState();
	ContactKktEquilibriumResult result;
	result.state = inputState;

	std::vector<ContactCandidate> activeContacts;
	VectorXd multipliers(0);
	if (m_kktWarmStart != nullptr)
	{
		activeContacts = m_kktWarmStart->activeContacts;
		multipliers = m_kktWarmStart->multipliers;
	}
	else
	{
		const EquilibriumResult barrierResult = solveStaticEquilibrium();
		if (!barrierResult.success)
		{
			result.rolledBack = true;
			return result;
		}

		const ContactDetectionResult barrierDetection =
			m_planarBarrierContactForce->detectCurrent();
		bool hasActiveBarrierContact = false;
		for (const ContactCandidate &candidate : barrierDetection.bodyCandidates)
		{
			const BarrierPotentialEvaluation barrier = evaluateBarrierPotential(
				candidate.gap, dBar, stiffness);
			hasActiveBarrierContact = hasActiveBarrierContact || barrier.active;
		}
		if (hasActiveBarrierContact)
		{
			const PlanarContactKktSeed seed = buildPlanarContactKktSeed();
			activeContacts = seed.contacts;
			multipliers = seed.linearization.multipliers;
		}
	}

	// Convert active-gap error to force units using the contact activation
	// length. Scaling by the full rod length underweights feasibility and makes
	// the filter reject useful KKT steps that rapidly close the gap.
	const double constraintForceScale = characteristicForce /
		std::max(dBar, std::numeric_limits<double>::min());
	for (int iteration = 0; iteration <= maxIter; ++iteration)
	{
		result.nonlinearIterations = iteration;
		ContactKktLinearization linearization =
			evaluatePlanarKktLinearization(activeContacts, multipliers);
		const ContactDetectionResult detection =
			m_planarBarrierContactForce->detectCurrent();

		const ContactActiveSetUpdate activeSet = updateContactActiveSet(
			activeContacts, multipliers, detection.bodyCandidates,
			kktGapTolerance, kktMultiplierTolerance);
		if (activeSet.changed())
		{
			activeContacts = activeSet.contacts;
			multipliers = activeSet.multipliers;
			++result.activeSetUpdates;
			if (result.activeSetUpdates > kktMaxActiveSetUpdates)
			{
				break;
			}
			continue;
		}

		const VectorXd stationarity = linearization.residual -
			linearization.constraintJacobian.transpose() * multipliers;
		result.stationarityNorm = stationarity.norm();
		result.activeConstraintNorm = linearization.gaps.norm();
		result.complementarityNorm = activeContacts.empty() ? 0.0 :
			(linearization.gaps.array() * multipliers.array()).matrix().norm();
		result.minimumBodyGap = detection.minimumBodyGap;
		result.tipClearance = detection.tip.clearance;
		result.minimumMultiplier = activeContacts.empty() ?
			std::numeric_limits<double>::infinity() : multipliers.minCoeff();
		result.tipSafe = detection.tip.safe;

		const bool inactiveFeasible =
			std::isinf(detection.minimumBodyGap) ||
			detection.minimumBodyGap >= -kktGapTolerance;
		const bool multipliersFeasible = activeContacts.empty() ||
			result.minimumMultiplier >= -kktMultiplierTolerance;
		if (result.stationarityNorm <= forceTol &&
			result.activeConstraintNorm <= kktGapTolerance &&
			result.complementarityNorm <= kktComplementarityTolerance &&
			inactiveFeasible && multipliersFeasible && result.tipSafe)
		{
			rod->commitStaticState();
			result.success = true;
			result.state = rod->captureState();
			result.activeContacts = activeContacts;
			result.multipliers = multipliers;
			return result;
		}
		if (iteration == maxIter || !result.tipSafe)
		{
			break;
		}

		const ContactKktSolution correction =
			solveContactKktSystem(linearization);
		if (!correction.success)
		{
			break;
		}

		const RodState iterationState = rod->captureState();
		const double currentMerit = 0.5 * (
			stationarity.squaredNorm() +
			std::pow(constraintForceScale, 2) *
				linearization.gaps.squaredNorm());
		double stepLength = 1.0;
		bool accepted = false;
		for (int trial = 0; trial <= maxLineSearchIter; ++trial)
		{
			rod->restoreState(iterationState);
			rod->applyFreeDofIncrement(
				stepLength * correction.configurationIncrement);
			const VectorXd trialMultipliers =
				multipliers + stepLength * correction.multiplierIncrement;
			std::vector<ContactCandidate> trialContacts = activeContacts;
			ContactKktLinearization trialLinearization =
				evaluatePlanarKktLinearization(trialContacts, trialMultipliers);
			const ContactDetectionResult trialDetection =
				m_planarBarrierContactForce->detectCurrent();
			const VectorXd trialStationarity = trialLinearization.residual -
				trialLinearization.constraintJacobian.transpose() * trialMultipliers;
			const double trialMerit = 0.5 * (
				trialStationarity.squaredNorm() +
				std::pow(constraintForceScale, 2) *
					trialLinearization.gaps.squaredNorm());

			const bool activeFeasibleTrial = trialLinearization.gaps.size() == 0 ||
				trialLinearization.gaps.minCoeff() >= -kktGapTolerance;
			if (trialDetection.tip.safe && activeFeasibleTrial &&
				std::isfinite(trialMerit) &&
				trialMerit <= (1.0 - lineSearchArmijo * stepLength) * currentMerit)
			{
				activeContacts = trialContacts;
				multipliers = trialMultipliers;
				accepted = true;
				break;
			}

			if (trial < maxLineSearchIter)
			{
				stepLength *= lineSearchReduction;
				++result.lineSearchBacktracks;
			}
		}
		if (!accepted)
		{
			rod->restoreState(iterationState);
			break;
		}
	}

	rod->restoreState(inputState);
	result.rolledBack = true;
	result.state = inputState;
	result.activeContacts = activeContacts;
	result.multipliers = multipliers;
	return result;
}

ContactKktEquilibriumResult world::correctPlanarContactKktEquilibrium(
	const ContactKktEquilibriumResult &warmStart)
{
	if (!warmStart.success)
	{
		throw invalid_argument("Planar KKT correction requires a converged warm start");
	}
	if (m_kktWarmStart != nullptr)
	{
		throw logic_error("Nested planar KKT warm-start correction is not supported");
	}
	rod->restoreState(warmStart.state);
	m_kktWarmStart = &warmStart;
	try
	{
		ContactKktEquilibriumResult result =
			solvePlanarContactKktEquilibrium();
		m_kktWarmStart = nullptr;
		return result;
	}
	catch (...)
	{
		m_kktWarmStart = nullptr;
		throw;
	}
}

ContactStabilityAnalysis world::analyzePlanarContactStability(
	const ContactKktEquilibriumResult &equilibrium)
{
	if (!equilibrium.success)
	{
		throw invalid_argument("Stability analysis requires a converged KKT equilibrium");
	}
	const RodState savedState = rod->captureState();
	try
	{
		rod->restoreState(equilibrium.state);
		std::vector<ContactCandidate> contacts = equilibrium.activeContacts;
		const ContactKktLinearization linearization =
			evaluatePlanarKktLinearization(contacts, equilibrium.multipliers);
		const ContactStabilityAnalysis result = analyzeContactStability(
			linearization.hessian, linearization.constraintJacobian);
		rod->restoreState(savedState);
		return result;
	}
	catch (...)
	{
		rod->restoreState(savedState);
		throw;
	}
}

ContactEquilibriumSensitivity world::computePlanarContactFieldSensitivity(
	const ContactKktEquilibriumResult &equilibrium)
{
	if (!equilibrium.success)
	{
		throw invalid_argument(
			"Contact field sensitivity requires a converged KKT equilibrium");
	}
	if (m_tipMagneticForce == nullptr)
	{
		throw logic_error(
			"Contact field sensitivity requires magneticModel=axial_tip");
	}

	const RodState savedState = rod->captureState();
	try
	{
		rod->restoreState(equilibrium.state);
		std::vector<ContactCandidate> contacts = equilibrium.activeContacts;
		const ContactKktLinearization linearization =
			evaluatePlanarKktLinearization(contacts, equilibrium.multipliers);
		const MatrixXd residualDerivative =
			m_tipMagneticForce->computeFreeDofResidualDerivativeField();
		const ConstrainedSensitivityResult freeSensitivity =
			solveConstrainedSensitivity(
				linearization.hessian, linearization.constraintJacobian,
				residualDerivative);

		ContactEquilibriumSensitivity result;
		result.success = freeSensitivity.success;
		result.activeContacts = contacts;
		result.multiplierDerivative = freeSensitivity.multiplierDerivative;
		result.linearResidualNorm = freeSensitivity.linearResidualNorm;
		result.configurationDerivative = MatrixXd::Zero(rod->ndof, 3);
		if (freeSensitivity.success)
		{
			for (int freeDof = 0; freeDof < rod->uncons; ++freeDof)
			{
				result.configurationDerivative.row(
					rod->unconstrainedMap[freeDof]) =
					freeSensitivity.configurationDerivative.row(freeDof);
			}
		}
		rod->restoreState(savedState);
		return result;
	}
	catch (...)
	{
		rod->restoreState(savedState);
		throw;
	}
}

ActuationEquilibriumSensitivity world::computePlanarContactActuationSensitivity(
	const ContactKktEquilibriumResult &equilibrium)
{
	if (!equilibrium.success)
	{
		throw invalid_argument(
			"Actuation sensitivity requires a converged KKT equilibrium");
	}
	if (m_tipMagneticForce == nullptr || m_insertionModel == nullptr)
	{
		throw logic_error(
			"Actuation sensitivity requires axial_tip magnetics and an insertion model");
	}

	const RodState savedState = rod->captureState();
	try
	{
		rod->restoreState(equilibrium.state);
		std::vector<ContactCandidate> contacts = equilibrium.activeContacts;
		const ContactKktLinearization linearization =
			evaluatePlanarKktLinearization(contacts, equilibrium.multipliers);
		MatrixXd residualDerivative = MatrixXd::Zero(rod->uncons, 4);
		residualDerivative.col(0) =
			m_insertionModel->freeResidualDerivativeCoordinate();
		residualDerivative.rightCols(3) =
			m_tipMagneticForce->computeFreeDofResidualDerivativeField();
		const ConstrainedSensitivityResult freeSensitivity =
			solveConstrainedSensitivity(
				linearization.hessian, linearization.constraintJacobian,
				residualDerivative);

		ActuationEquilibriumSensitivity result;
		result.success = freeSensitivity.success;
		result.activeContacts = contacts;
		result.multiplierDerivative = freeSensitivity.multiplierDerivative;
		result.linearResidualNorm = freeSensitivity.linearResidualNorm;
		result.configurationDerivative = MatrixXd::Zero(rod->ndof, 4);
		if (freeSensitivity.success)
		{
			for (int freeDof = 0; freeDof < rod->uncons; ++freeDof)
			{
				result.configurationDerivative.row(
					rod->unconstrainedMap[freeDof]) =
					freeSensitivity.configurationDerivative.row(freeDof);
			}
		}
		rod->restoreState(savedState);
		return result;
	}
	catch (...)
	{
		rod->restoreState(savedState);
		throw;
	}
}

Actuation world::getActuation() const
{
	Actuation actuation;
	actuation.xi = m_insertionModel == nullptr ? 0.0 :
		m_insertionModel->coordinate();
	actuation.field = baVector;
	return actuation;
}

void world::setActuation(const Actuation &actuation)
{
	if (!std::isfinite(actuation.xi) || actuation.xi < 0.0 ||
		!actuation.field.allFinite())
	{
		throw invalid_argument("Actuation must contain finite, nonnegative insertion");
	}
	if (m_insertionModel == nullptr)
	{
		if (actuation.xi != 0.0)
		{
			throw logic_error("Nonzero insertion requires an insertion model");
		}
	}
	else
	{
		m_insertionModel->setCoordinate(actuation.xi);
		insertionCoordinate = actuation.xi;
	}
	setAppliedField(actuation.field);
}

PlannerState world::capturePlannerState(
	const ContactKktEquilibriumResult &equilibrium,
	const ContactStabilityAnalysis &stability) const
{
	if (!equilibrium.success || !stability.valid)
	{
		throw invalid_argument(
			"PlannerState requires a converged equilibrium and valid stability analysis");
	}
	PlannerState state;
	state.rodState = equilibrium.state;
	state.actuation = getActuation();
	state.activeContacts = equilibrium.activeContacts;
	state.multipliers = equilibrium.multipliers;
	state.stabilityMargin = stability.minimumEigenvalue;
	state.stationarityNorm = equilibrium.stationarityNorm;
	state.activeConstraintNorm = equilibrium.activeConstraintNorm;
	state.complementarityNorm = equilibrium.complementarityNorm;
	state.minimumBodyGap = equilibrium.minimumBodyGap;
	state.tipClearance = equilibrium.tipClearance;
	state.tipSafe = equilibrium.tipSafe;
	return state;
}

void world::restorePlannerState(const PlannerState &state)
{
	setActuation(state.actuation);
	rod->restoreState(state.rodState);
}

FieldContinuationResult world::continuePlanarContactField(
	const Vector3d &targetField,
	const FieldContinuationOptions &options)
{
	if (m_tipMagneticForce == nullptr)
	{
		throw logic_error(
			"Field continuation requires magneticModel=axial_tip");
	}
	if (m_planarBarrierContactForce == nullptr)
	{
		throw logic_error(
			"Planar contact continuation requires contactModel=planar_barrier");
	}
	if (!targetField.allFinite() ||
		!std::isfinite(options.initialStepFraction) ||
		!std::isfinite(options.minimumStepFraction) ||
		!std::isfinite(options.maximumStepFraction) ||
		!std::isfinite(options.stepReduction) ||
		!std::isfinite(options.stepGrowth) ||
		!std::isfinite(options.stabilityTolerance) ||
		options.initialStepFraction <= 0.0 ||
		options.minimumStepFraction <= 0.0 ||
		options.maximumStepFraction < options.minimumStepFraction ||
		options.initialStepFraction > options.maximumStepFraction ||
		options.stepReduction <= 0.0 || options.stepReduction >= 1.0 ||
		options.stepGrowth <= 1.0 || options.easyCorrectorIterations < 0 ||
		options.maximumAttempts <= 0)
	{
		throw invalid_argument("Invalid field-continuation options");
	}

	const RodState inputState = rod->captureState();
	const Vector3d startField = baVector;
	FieldContinuationResult result;
	result.startField = startField;
	result.targetField = targetField;

	auto rollBack = [&]() {
		rod->restoreState(inputState);
		setAppliedField(startField);
		result.rolledBack = true;
	};
	auto sameId = [](const ContactCandidate &left,
					 const ContactCandidate &right) {
		return left.rodVertex == right.rodVertex &&
			left.boundaryId == right.boundaryId;
	};
	auto countMissing = [&sameId](
		const std::vector<ContactCandidate> &source,
		const std::vector<ContactCandidate> &destination) {
		int missing = 0;
		for (const ContactCandidate &candidate : source)
		{
			bool found = false;
			for (const ContactCandidate &other : destination)
			{
				found = found || sameId(candidate, other);
			}
			missing += found ? 0 : 1;
		}
		return missing;
	};

	ContactKktEquilibriumResult current =
		solvePlanarContactKktEquilibrium();
	if (!current.success)
	{
		rollBack();
		return result;
	}
	ContactStabilityAnalysis currentStability =
		analyzePlanarContactStability(current);
	if (!currentStability.valid ||
		currentStability.minimumEigenvalue <= options.stabilityTolerance)
	{
		rollBack();
		return result;
	}

	FieldContinuationPoint initialPoint;
	initialPoint.field = startField;
	initialPoint.equilibrium = current;
	initialPoint.stability = currentStability;
	result.points.push_back(initialPoint);

	const Vector3d totalFieldChange = targetField - startField;
	if (totalFieldChange.norm() == 0.0)
	{
		result.success = true;
		return result;
	}

	double pathFraction = 0.0;
	double stepFraction = options.initialStepFraction;
	while (pathFraction < 1.0 && result.attemptedSteps < options.maximumAttempts)
	{
		stepFraction = std::min(stepFraction, 1.0 - pathFraction);
		result.minimumAttemptedStepFraction = std::min(
			result.minimumAttemptedStepFraction, stepFraction);
		++result.attemptedSteps;
		const Vector3d fieldIncrement = stepFraction * totalFieldChange;
		const Vector3d trialField =
			startField + (pathFraction + stepFraction) * totalFieldChange;

		const ContactEquilibriumSensitivity sensitivity =
			computePlanarContactFieldSensitivity(current);
		if (!sensitivity.success)
		{
			rollBack();
			return result;
		}
		ContactKktEquilibriumResult predicted = current;
		predicted.state.configuration +=
			sensitivity.configurationDerivative * fieldIncrement;
		predicted.multipliers +=
			sensitivity.multiplierDerivative * fieldIncrement;

		setAppliedField(trialField);
		const ContactKktEquilibriumResult corrected =
			correctPlanarContactKktEquilibrium(predicted);
		bool accepted = corrected.success;
		ContactStabilityAnalysis correctedStability;
		if (accepted)
		{
			correctedStability = analyzePlanarContactStability(corrected);
			accepted = correctedStability.valid &&
				correctedStability.minimumEigenvalue > options.stabilityTolerance;
		}

		if (!accepted)
		{
			++result.rejectedSteps;
			rod->restoreState(current.state);
			setAppliedField(startField + pathFraction * totalFieldChange);
			stepFraction *= options.stepReduction;
			if (stepFraction < options.minimumStepFraction)
			{
				rollBack();
				return result;
			}
			continue;
		}

		FieldContinuationPoint point;
		point.field = trialField;
		point.equilibrium = corrected;
		point.stability = correctedStability;
		point.pathFraction = pathFraction + stepFraction;
		point.acceptedStepFraction = stepFraction;
		point.configurationPredictorError =
			(corrected.state.configuration - predicted.state.configuration).norm();
		point.contactsAdded = countMissing(
			corrected.activeContacts, current.activeContacts);
		point.contactsReleased = countMissing(
			current.activeContacts, corrected.activeContacts);
		point.contactSetChanged =
			point.contactsAdded > 0 || point.contactsReleased > 0;
		double multiplierErrorSquared = 0.0;
		for (int correctedContact = 0;
			 correctedContact < static_cast<int>(corrected.activeContacts.size());
			 ++correctedContact)
		{
			bool matched = false;
			for (int predictedContact = 0;
				 predictedContact < static_cast<int>(predicted.activeContacts.size());
				 ++predictedContact)
			{
				if (sameId(corrected.activeContacts[correctedContact],
						predicted.activeContacts[predictedContact]))
				{
					const double difference = corrected.multipliers[correctedContact] -
						predicted.multipliers[predictedContact];
					multiplierErrorSquared += difference * difference;
					matched = true;
					break;
				}
			}
			if (!matched)
			{
				multiplierErrorSquared +=
					std::pow(corrected.multipliers[correctedContact], 2);
			}
		}
		for (int predictedContact = 0;
			 predictedContact < static_cast<int>(predicted.activeContacts.size());
			 ++predictedContact)
		{
			bool matched = false;
			for (const ContactCandidate &correctedContact : corrected.activeContacts)
			{
				matched = matched || sameId(
					predicted.activeContacts[predictedContact], correctedContact);
			}
			if (!matched)
			{
				multiplierErrorSquared +=
					std::pow(predicted.multipliers[predictedContact], 2);
			}
		}
		point.multiplierPredictorError = std::sqrt(multiplierErrorSquared);
		result.points.push_back(point);

		pathFraction += stepFraction;
		current = corrected;
		currentStability = correctedStability;
		if (!point.contactSetChanged &&
			corrected.nonlinearIterations <= options.easyCorrectorIterations)
		{
			stepFraction = std::min(
				options.maximumStepFraction,
				stepFraction * options.stepGrowth);
		}
	}

	if (pathFraction >= 1.0)
	{
		setAppliedField(targetField);
		result.success = true;
		return result;
	}
	rollBack();
	return result;
}

ActuationContinuationResult world::continuePlanarContactActuation(
	const Actuation &targetActuation,
	const FieldContinuationOptions &options)
{
	return continuePlanarContactActuationImpl(
		targetActuation, options, nullptr);
}

ActuationContinuationResult world::continuePlanarContactActuation(
	const PlannerState &startState,
	const Actuation &targetActuation,
	const FieldContinuationOptions &options)
{
	if (startState.multipliers.size() !=
		static_cast<int>(startState.activeContacts.size()))
	{
		throw invalid_argument(
			"PlannerState contact and multiplier counts must agree");
	}
	restorePlannerState(startState);
	ContactKktEquilibriumResult warmStart;
	warmStart.success = true;
	warmStart.state = startState.rodState;
	warmStart.activeContacts = startState.activeContacts;
	warmStart.multipliers = startState.multipliers;
	warmStart.stationarityNorm = startState.stationarityNorm;
	warmStart.activeConstraintNorm = startState.activeConstraintNorm;
	warmStart.complementarityNorm = startState.complementarityNorm;
	warmStart.minimumBodyGap = startState.minimumBodyGap;
	warmStart.tipClearance = startState.tipClearance;
	warmStart.tipSafe = startState.tipSafe;
	return continuePlanarContactActuationImpl(
		targetActuation, options, &warmStart);
}

ActuationContinuationResult world::continuePlanarContactActuationImpl(
	const Actuation &targetActuation,
	const FieldContinuationOptions &options,
	const ContactKktEquilibriumResult *warmStart)
{
	const Actuation startActuation = getActuation();
	if (m_insertionModel == nullptr || m_tipMagneticForce == nullptr ||
		m_planarBarrierContactForce == nullptr)
	{
		throw logic_error(
			"Actuation continuation requires proximal insertion, axial_tip, and planar_barrier");
	}
	if (!std::isfinite(targetActuation.xi) ||
		targetActuation.xi < startActuation.xi ||
		!targetActuation.field.allFinite() ||
		!std::isfinite(options.initialStepFraction) ||
		!std::isfinite(options.minimumStepFraction) ||
		!std::isfinite(options.maximumStepFraction) ||
		!std::isfinite(options.stepReduction) ||
		!std::isfinite(options.stepGrowth) ||
		!std::isfinite(options.stabilityTolerance) ||
		options.initialStepFraction <= 0.0 ||
		options.minimumStepFraction <= 0.0 ||
		options.maximumStepFraction < options.minimumStepFraction ||
		options.initialStepFraction > options.maximumStepFraction ||
		options.stepReduction <= 0.0 || options.stepReduction >= 1.0 ||
		options.stepGrowth <= 1.0 || options.easyCorrectorIterations < 0 ||
		options.maximumAttempts <= 0)
	{
		throw invalid_argument(
			"Invalid or retracting actuation-continuation request");
	}

	const RodState inputState = rod->captureState();
	ActuationContinuationResult result;
	result.startActuation = startActuation;
	result.targetActuation = targetActuation;
	auto rollBack = [&]() {
		setActuation(startActuation);
		rod->restoreState(inputState);
		result.rolledBack = true;
	};
	auto sameId = [](const ContactCandidate &left, const ContactCandidate &right) {
		return left.rodVertex == right.rodVertex &&
			left.boundaryId == right.boundaryId;
	};
	auto countMissing = [&sameId](
		const std::vector<ContactCandidate> &source,
		const std::vector<ContactCandidate> &destination) {
		int missing = 0;
		for (const ContactCandidate &candidate : source)
		{
			bool found = false;
			for (const ContactCandidate &other : destination)
			{
				found = found || sameId(candidate, other);
			}
			missing += found ? 0 : 1;
		}
		return missing;
	};
	auto multiplierPredictionError = [&sameId](
		const ContactKktEquilibriumResult &corrected,
		const ContactKktEquilibriumResult &predicted) {
		double squaredError = 0.0;
		for (int i = 0; i < static_cast<int>(corrected.activeContacts.size()); ++i)
		{
			bool matched = false;
			for (int j = 0; j < static_cast<int>(predicted.activeContacts.size()); ++j)
			{
				if (sameId(corrected.activeContacts[i], predicted.activeContacts[j]))
				{
					const double difference =
						corrected.multipliers[i] - predicted.multipliers[j];
					squaredError += difference * difference;
					matched = true;
					break;
				}
			}
			if (!matched)
			{
				squaredError += corrected.multipliers[i] * corrected.multipliers[i];
			}
		}
		for (int j = 0; j < static_cast<int>(predicted.activeContacts.size()); ++j)
		{
			bool matched = false;
			for (const ContactCandidate &candidate : corrected.activeContacts)
			{
				matched = matched || sameId(candidate, predicted.activeContacts[j]);
			}
			if (!matched)
			{
				squaredError += predicted.multipliers[j] * predicted.multipliers[j];
			}
		}
		return std::sqrt(squaredError);
	};

	ContactKktEquilibriumResult current = warmStart == nullptr ?
		solvePlanarContactKktEquilibrium() :
		correctPlanarContactKktEquilibrium(*warmStart);
	if (!current.success)
	{
		result.failureReason =
			std::isfinite(current.tipClearance) && !current.tipSafe ?
			"tip_safety" : "initial_equilibrium";
		rollBack();
		return result;
	}
	ContactStabilityAnalysis stability = analyzePlanarContactStability(current);
	if (!stability.valid || stability.minimumEigenvalue <= options.stabilityTolerance)
	{
		result.failureReason = "stability";
		rollBack();
		return result;
	}
	ActuationContinuationPoint initial;
	initial.actuation = startActuation;
	initial.equilibrium = current;
	initial.stability = stability;
	result.points.push_back(initial);

	Vector4d totalChange;
	totalChange << targetActuation.xi - startActuation.xi,
		(targetActuation.field - startActuation.field);
	if (totalChange.norm() == 0.0)
	{
		result.success = true;
		return result;
	}

	double pathFraction = 0.0;
	double stepFraction = options.initialStepFraction;
	string lastRejectionReason;
	while (pathFraction < 1.0 && result.attemptedSteps < options.maximumAttempts)
	{
		stepFraction = std::min(stepFraction, 1.0 - pathFraction);
		result.minimumAttemptedStepFraction = std::min(
			result.minimumAttemptedStepFraction, stepFraction);
		++result.attemptedSteps;
		const Vector4d increment = stepFraction * totalChange;
		Actuation trialActuation;
		trialActuation.xi = startActuation.xi +
			(pathFraction + stepFraction) * totalChange[0];
		trialActuation.field = startActuation.field +
			(pathFraction + stepFraction) * totalChange.tail<3>();

		const ActuationEquilibriumSensitivity sensitivity =
			computePlanarContactActuationSensitivity(current);
		if (!sensitivity.success)
		{
			result.failureReason = "sensitivity";
			rollBack();
			return result;
		}
		ContactKktEquilibriumResult predicted = current;
		predicted.state.configuration +=
			sensitivity.configurationDerivative * increment;
		predicted.multipliers += sensitivity.multiplierDerivative * increment;
		setActuation(trialActuation);
		const ContactKktEquilibriumResult corrected =
			correctPlanarContactKktEquilibrium(predicted);
		bool accepted = corrected.success;
		ContactStabilityAnalysis correctedStability;
		if (accepted)
		{
			correctedStability = analyzePlanarContactStability(corrected);
			accepted = correctedStability.valid &&
				correctedStability.minimumEigenvalue > options.stabilityTolerance;
			if (!accepted)
			{
				lastRejectionReason = "stability";
			}
		}
		else
		{
			lastRejectionReason =
				std::isfinite(corrected.tipClearance) && !corrected.tipSafe ?
				"tip_safety" : "corrector";
		}
		if (!accepted)
		{
			++result.rejectedSteps;
			Actuation acceptedActuation;
			acceptedActuation.xi = startActuation.xi +
				pathFraction * totalChange[0];
			acceptedActuation.field = startActuation.field +
				pathFraction * totalChange.tail<3>();
			setActuation(acceptedActuation);
			rod->restoreState(current.state);
			stepFraction *= options.stepReduction;
			if (stepFraction < options.minimumStepFraction)
			{
				result.failureReason = lastRejectionReason;
				rollBack();
				return result;
			}
			continue;
		}

		ActuationContinuationPoint point;
		point.actuation = trialActuation;
		point.equilibrium = corrected;
		point.stability = correctedStability;
		point.pathFraction = pathFraction + stepFraction;
		point.acceptedStepFraction = stepFraction;
		point.configurationPredictorError =
			(corrected.state.configuration - predicted.state.configuration).norm();
		point.multiplierPredictorError =
			multiplierPredictionError(corrected, predicted);
		point.contactsAdded = countMissing(
			corrected.activeContacts, current.activeContacts);
		point.contactsReleased = countMissing(
			current.activeContacts, corrected.activeContacts);
		point.contactSetChanged =
			point.contactsAdded > 0 || point.contactsReleased > 0;
		result.points.push_back(point);
		pathFraction += stepFraction;
		current = corrected;
		if (!point.contactSetChanged &&
			corrected.nonlinearIterations <= options.easyCorrectorIterations)
		{
			stepFraction = std::min(
				options.maximumStepFraction,
				stepFraction * options.stepGrowth);
		}
	}

	if (pathFraction >= 1.0)
	{
		setActuation(targetActuation);
		result.success = true;
		return result;
	}
	result.failureReason = lastRejectionReason.empty() ?
		"attempt_limit" : lastRejectionReason;
	rollBack();
	return result;
}

void world::setAppliedField(const Vector3d &field)
{
	if (m_tipMagneticForce == nullptr)
	{
		throw logic_error(
			"Runtime field updates require magneticModel=axial_tip");
	}
	baVector = field;
	m_tipMagneticForce->setField(field);
}

Vector3d world::getAppliedField() const
{
	return baVector;
}

MatrixXd world::computeConfigurationFieldSensitivity()
{
	if (m_tipMagneticForce == nullptr)
	{
		throw logic_error(
			"Field sensitivity requires magneticModel=axial_tip");
	}

	const StaticEvaluation evaluation = evaluateStaticSystem();
	const MatrixXd residualDerivative =
		m_tipMagneticForce->computeFreeDofResidualDerivativeField();
	MatrixXd fullSensitivity = MatrixXd::Zero(rod->ndof, 3);

	for (int component = 0; component < 3; ++component)
	{
		StaticEvaluation rightHandSide = evaluation;
		rightHandSide.residual = -residualDerivative.col(component);
		VectorXd freeSensitivity;
		const int solveInfo =
			stepper->solveBandedSystem(rightHandSide, freeSensitivity);
		if (solveInfo != 0 || !freeSensitivity.allFinite())
		{
			throw runtime_error(
				"Failed to solve equilibrium field-sensitivity system");
		}
		for (int freeIndex = 0; freeIndex < rod->uncons; ++freeIndex)
		{
			fullSensitivity(rod->unconstrainedMap[freeIndex], component) =
				freeSensitivity[freeIndex];
		}
	}
	return fullSensitivity;
}

RodState world::captureRodState() const
{
	return rod->captureState();
}

void world::restoreRodState(const RodState &state)
{
	rod->restoreState(state);
}

void world::applyFreeDofIncrement(const VectorXd &increment)
{
	rod->applyFreeDofIncrement(increment);
}

int world::numFreeDofs() const
{
	return rod->uncons;
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

double world::getVelocityNorm() const
{
	return rod->u.norm();
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
