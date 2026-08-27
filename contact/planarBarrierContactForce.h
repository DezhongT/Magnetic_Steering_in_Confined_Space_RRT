#ifndef PLANARBARRIERCONTACTFORCE_H
#define PLANARBARRIERCONTACTFORCE_H

#include "contact/barrierPotential.h"
#include "contact/contactCandidate.h"
#include "contact/contactDetector.h"
#include "elasticRod.h"
#include "geometry/confinedDomain.h"
#include "timeStepper.h"

struct PlanarBarrierContactEvaluation
{
	bool valid = false;
	bool active = false;
	double energy = 0.0;
	Vector3d residual = Vector3d::Zero();
	Matrix3d jacobian = Matrix3d::Zero();
};

PlanarBarrierContactEvaluation evaluatePlanarBarrierContact(
	const ContactCandidate &candidate,
	double activationDistance,
	double stiffness);

class PlanarBarrierContactForce
{
public:
	PlanarBarrierContactForce(
		elasticRod &rod,
		timeStepper &stepper,
		const ConfinedDomain &domain,
		double rodRadius,
		double activationDistance,
		double stiffness,
		double tipSafeDistance);

	void computeFc();
	ContactDetectionResult detectCurrent() const;
	const ContactDetectionResult &lastDetection() const;

private:
	elasticRod *rod;
	timeStepper *stepper;
	ContactDetector detector;
	double activationDistance;
	double stiffness;
	ContactDetectionResult detection;
};

#endif
