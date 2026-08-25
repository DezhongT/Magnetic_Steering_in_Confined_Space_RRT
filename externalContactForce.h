#ifndef EXTERNALCONTACTFORCE_H
#define EXTERNALCONTACTFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class externalContactForce
{
public:
	externalContactForce(elasticRod &m_rod, timeStepper &m_stepper, 
        double m_dBar, double m_stiffness, double m_thickness);
	~externalContactForce();

	void computeFc();
	void computeJc();
    
private:
	elasticRod *rod;
	timeStepper *stepper;

    double thick;
    double db;
    double stiff;

   
    Vector3d p1;
    Vector3d p2;
    Vector3d p3;

    Vector3d point;


    Vector3d computeContactForce(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double p3x, double p3y, double p3z, 
  double x, double y, double z, double dbar, double Kc, double thickness);
    Matrix3d computeContactJacobian(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double p3x, double p3y, double p3z, 
  double x, double y, double z, double dbar, double Kc, double thickness);
    Vector3d ListVec(double a1, double a2, double a3);
    Matrix3d ListMat(Vector3d a1, Vector3d a2, Vector3d a3);
};

#endif
