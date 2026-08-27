#ifndef CONFINEDDOMAIN_H
#define CONFINEDDOMAIN_H

#include "geometry/domainQuery.h"

class ConfinedDomain
{
public:
	virtual ~ConfinedDomain() = default;
	virtual DomainQuery query(const Vector3d &point) const = 0;
	virtual DomainQuery queryBoundary(
		const Vector3d &point, int boundaryId) const = 0;
};

#endif
