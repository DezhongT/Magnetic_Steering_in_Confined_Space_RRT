#include "timeStepper.h"

timeStepper::timeStepper(elasticRod &m_rod)
{
	rod = &m_rod;
	kl = 10; // lower diagonals
	ku = 10; // upper diagonals
	freeDOF = rod->uncons;
	ldb = freeDOF;
	NUMROWS = 2 * kl + ku + 1;
	totalForce = new double[freeDOF];
	jacobianLen = (2 * kl + ku + 1) * freeDOF;
	jacobian = new double [jacobianLen];
	nrhs = 1;
    ipiv = new int[freeDOF];
    info = 0;
}

timeStepper::~timeStepper()
{
	delete[] totalForce;
	delete[] jacobian;
	delete[] ipiv;
}

double* timeStepper::getForce()
{
	return totalForce;
}

double* timeStepper::getJacobian()
{
	return jacobian;
}

void timeStepper::addForce(int ind, double p)
{
	if (rod->getIfConstrained(ind) == 0) // free dof
	{
		mappedInd = rod->fullToUnconsMap[ind];
		totalForce[mappedInd] = totalForce[mappedInd] + p; // subtracting elastic force
	}
}

void timeStepper::addJacobian(int ind1, int ind2, double p)
{
	mappedInd1 = rod->fullToUnconsMap[ind1];
	mappedInd2 = rod->fullToUnconsMap[ind2];
	if (rod->getIfConstrained(ind1) == 0 && rod->getIfConstrained(ind2) == 0) // both are free
	{
		row = kl + ku + mappedInd2 - mappedInd1;
        col = mappedInd1;
        offset = row + col * NUMROWS;
        jacobian[offset] = jacobian[offset] + p;
	}
}

void timeStepper::setZero()
{
	for (int i=0; i < freeDOF; i++)
		totalForce[i] = 0;
	for (int i=0; i < jacobianLen; i++)
		jacobian[i] = 0;
}

StaticEvaluation timeStepper::captureEvaluation() const
{
	StaticEvaluation evaluation;
	evaluation.residual = Map<const VectorXd>(totalForce, freeDOF);
	evaluation.bandedJacobian =
		Map<const Matrix<double, Dynamic, Dynamic, ColMajor>>(jacobian, NUMROWS, freeDOF);
	evaluation.lowerBandwidth = kl;
	evaluation.upperBandwidth = ku;
	return evaluation;
}

int timeStepper::solveBandedSystem(
	const StaticEvaluation &evaluation,
	VectorXd &solution) const
{
	const int systemSize = static_cast<int>(evaluation.residual.size());
	const int expectedRows =
		2 * evaluation.lowerBandwidth + evaluation.upperBandwidth + 1;
	if (systemSize != freeDOF ||
		evaluation.bandedJacobian.rows() != expectedRows ||
		evaluation.bandedJacobian.cols() != systemSize)
	{
		return -100;
	}

	MatrixXd jacobianWorkspace = evaluation.bandedJacobian;
	solution = evaluation.residual;
	VectorXi pivotWorkspace(systemSize);

	int n = systemSize;
	int lower = evaluation.lowerBandwidth;
	int upper = evaluation.upperBandwidth;
	int rightHandSides = 1;
	int leadingDimension = expectedRows;
	int rightHandSideLeadingDimension = systemSize;
	int solveInfo = 0;

	dgbsv_(&n, &lower, &upper, &rightHandSides,
		jacobianWorkspace.data(), &leadingDimension, pivotWorkspace.data(),
		solution.data(), &rightHandSideLeadingDimension, &solveInfo);
	return solveInfo;
}

int timeStepper::integrator()
{
	dgbsv_(&freeDOF, &kl, &ku, &nrhs, jacobian, &NUMROWS, ipiv, totalForce, &ldb, &info);
	return info;
}
