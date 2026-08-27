#ifndef CONTACTDETECTOR_H
#define CONTACTDETECTOR_H

#include "contact/contactCandidate.h"
#include "geometry/confinedDomain.h"

#include <vector>

class ContactDetector
{
public:
	ContactDetector(
		const ConfinedDomain &domain,
		double rodRadius,
		double activationDistance,
		double tipSafeDistance);

	ContactDetectionResult detect(
		const std::vector<Vector3d> &rodVertices) const;

private:
	const ConfinedDomain *domain;
	double rodRadius;
	double activationDistance;
	double tipSafeDistance;
};

#endif
