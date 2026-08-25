#include "dampingForce.h"
#include <iostream>

dampingForce::dampingForce(elasticRod &m_rod, timeStepper &m_stepper, double m_viscosity)
{
	rod = &m_rod;
	stepper = &m_stepper;
	viscosity = m_viscosity;
	dt = rod->dt;
         
}

dampingForce::~dampingForce()
{
	;
}

void dampingForce::computeFd()
{
	for (int i = 0; i < rod->ndof; i++)
	{
		double localForce = - viscosity * (rod->x(i) - rod->x0(i)) / rod->dt;

		stepper->addForce(i, - localForce); // subtracting external force
	}
}

void dampingForce::computeJd()
{
	for (int i = 0; i < rod->ndof; i++)
	{
		double localJacobian = - viscosity / rod->dt;

		stepper->addJacobian(i, i, - localJacobian);
	}
}
