#ifndef INSERTIONMODEL_H
#define INSERTIONMODEL_H

#include "eigenIncludes.h"

#include <string>

class InsertionModel
{
public:
	virtual ~InsertionModel() = default;
	virtual void setCoordinate(double xi) = 0;
	virtual double coordinate() const = 0;
	virtual void assemble() = 0;
	virtual VectorXd freeResidualDerivativeCoordinate() const = 0;
	virtual Vector3d guideTarget() const = 0;
	virtual std::string modelName() const = 0;
};

#endif
