#ifndef ACTUATION_H
#define ACTUATION_H

#include "eigenIncludes.h"

struct Actuation
{
	double xi = 0.0;
	Vector3d field = Vector3d::Zero();
};

#endif
