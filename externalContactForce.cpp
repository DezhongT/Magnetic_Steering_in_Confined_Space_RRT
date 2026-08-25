#include "externalContactForce.h"

externalContactForce::externalContactForce(elasticRod &m_rod, timeStepper &m_stepper, 
        double m_dBar, double m_stiffness, double m_thickness)
{
	rod = &m_rod;
	stepper = &m_stepper;

	db = m_dBar;
	stiff = m_stiffness; 
	thick = m_thickness;
}

externalContactForce::~externalContactForce()
{
	;
}

void externalContactForce::computeFc()
{
	for (int i = 0; i < rod->nv; i++)
	{
		Vector3d xCurrent = rod->getVertex(i);

		double minDistance = 100000.00;
		int contactIndex = 0;

		for (int kk = 0; kk < rod->plateTri; kk++)
		{
      Vector3i plateIndex = rod->v_triangular[kk];

			Vector3d plateNode1 = rod->v_nodes[plateIndex(0)];
			Vector3d plateNode2 = rod->v_nodes[plateIndex(1)];
      Vector3d plateNode3 = rod->v_nodes[plateIndex(2)];

			double d1 = (xCurrent - plateNode1).norm();
			double d2 = (xCurrent - plateNode2).norm();
      double d3 = (xCurrent - plateNode3).norm();

			double xCurrentDis = (d1 + d2 + d3) / 3;

			if (xCurrentDis < minDistance)
			{
        minDistance = xCurrentDis;
        contactIndex = kk;
			}
		}

		point = xCurrent;

    Vector3i plateContactIndex = rod->v_triangular[contactIndex];

		p1 = rod->v_nodes[plateContactIndex(0)];
		p2 = rod->v_nodes[plateContactIndex(1)];
    p3 = rod->v_nodes[plateContactIndex(2)];

		Vector3d e1 = p2 - p1;
		Vector3d e2 = p3 - p1;
    Vector3d pp = point - p1;

    Vector3d surfN = e1.cross(e2);

    surfN = surfN / surfN.norm();

		double dis = pp.dot(surfN);

    if (dis < 0)
    {
      p1 = rod->v_nodes[plateContactIndex(0)];
      p2 = rod->v_nodes[plateContactIndex(2)];
      p3 = rod->v_nodes[plateContactIndex(1)];

      e1 = p2 - p1;
      e2 = p3 - p1;
      pp = point - p1;

      surfN = e1.cross(e2);
      surfN = surfN / surfN.norm();

      dis = pp.dot(surfN);
    }

    double delta = thick - dis;

		//cout << dis << endl;

		if (delta < db)
		{
			Vector3d force;

      force = - computeContactForce(p1(0), p1(1), p1(2), p2(0), p2(1), p2(2), p3(0), p3(1), p3(2),
          point(0), point(1), point(2), db, stiff, thick);

			for (int ii = 0; ii < 3; ii++)
			{
				int ind = 4 * i + ii;
				stepper->addForce(ind, - force(ii));
			}

  


			Matrix3d jaco;

      jaco = computeContactJacobian(p1(0), p1(1), p1(2), p2(0), p2(1), p2(2), p3(0), p3(1), p3(2),
          point(0), point(1), point(2), db, stiff, thick);

			for (int ii = 0; ii < 3; ii++)
			{
				for (int jj = 0; jj < 3; jj++)
				{
					int ind1 = 4 * i + ii;
					int ind2 = 4 * i + jj;

					stepper->addJacobian(ind1, ind2, jaco(ii, jj));
				}
			}
		}

	}
}


Vector3d externalContactForce::computeContactForce(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double p3x, double p3y, double p3z, 
  double x, double y, double z, double dbar, double Kc, double thickness)
{
  Vector3d vecResult;

  vecResult = ListVec((Kc*(-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
       pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
               p1y*(-p2z + p3z))*(p1x - x) + 
            (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
            (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
          sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
            pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
            pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
     (sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
         pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
         pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
       (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
             (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
             (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
             (p1z - z))/
          sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
            pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
            pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
    (2*Kc*(-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
       (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
             (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
             (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
             (p1z - z))/
          sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
            pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
            pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))*
       log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
               (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
               (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
               (p1z - z))/
            sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
              pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
              pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar))
      /sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
       pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
       pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),
   (Kc*(p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
       pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
               p1y*(-p2z + p3z))*(p1x - x) + 
            (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
            (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
          sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
            pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
            pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
     (sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
         pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
         pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
       (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
             (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
             (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
             (p1z - z))/
          sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
            pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
            pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
    (2*Kc*(p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
       (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
             (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
             (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
             (p1z - z))/
          sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
            pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
            pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))*
       log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
               (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
               (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
               (p1z - z))/
            sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
              pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
              pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar))
      /sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
       pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
       pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),
   (Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
       pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
               p1y*(-p2z + p3z))*(p1x - x) + 
            (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
            (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
          sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
            pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
            pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
     (sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
         pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
         pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
       (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
             (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
             (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
             (p1z - z))/
          sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
            pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
            pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
    (2*Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
       (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
             (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
             (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
             (p1z - z))/
          sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
            pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
            pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))*
       log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
               (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
               (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
               (p1z - z))/
            sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
              pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
              pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar))
      /sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
       pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
       pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)));

  return vecResult;
}

Matrix3d externalContactForce::computeContactJacobian(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double p3x, double p3y, double p3z, 
  double x, double y, double z, double dbar, double Kc, double thickness)
{
  Matrix3d matResult;

  matResult = ListMat(ListVec((-4*Kc*pow(-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z),2)*
        (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
     (Kc*pow(-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z),2)*
        pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
                p1y*(-p2z + p3z))*(p1x - x) + 
             (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
             (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        pow(thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2)) - 
     (2*Kc*pow(-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z),2)*
        log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
                (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
                (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
                (p1z - z))/
             sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
               pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
               pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar)
        )/(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
        pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
        pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),
    (-4*Kc*(p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
     (Kc*(p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
                p1y*(-p2z + p3z))*(p1x - x) + 
             (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
             (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        pow(thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2)) - 
     (2*Kc*(p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
                (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
                (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
                (p1z - z))/
             sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
               pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
               pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar)
        )/(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
        pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
        pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),
    (-4*Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
     (Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
                p1y*(-p2z + p3z))*(p1x - x) + 
             (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
             (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        pow(thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2)) - 
     (2*Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
                (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
                (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
                (p1z - z))/
             sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
               pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
               pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar)
        )/(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
        pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
        pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))),
   ListVec((-4*Kc*(p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
     (Kc*(p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
                p1y*(-p2z + p3z))*(p1x - x) + 
             (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
             (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        pow(thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2)) - 
     (2*Kc*(p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
                (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
                (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
                (p1z - z))/
             sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
               pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
               pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar)
        )/(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
        pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
        pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),
    (-4*Kc*pow(p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z),2)*
        (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
     (Kc*pow(p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z),2)*
        pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
                p1y*(-p2z + p3z))*(p1x - x) + 
             (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
             (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        pow(thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2)) - 
     (2*Kc*pow(p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z),2)*
        log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
                (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
                (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
                (p1z - z))/
             sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
               pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
               pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar)
        )/(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
        pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
        pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),
    (-4*Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
     (Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
                p1y*(-p2z + p3z))*(p1x - x) + 
             (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
             (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        pow(thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2)) - 
     (2*Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
                (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
                (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
                (p1z - z))/
             sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
               pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
               pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar)
        )/(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
        pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
        pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))),
   ListVec((-4*Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
     (Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
                p1y*(-p2z + p3z))*(p1x - x) + 
             (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
             (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        pow(thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2)) - 
     (2*Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (-(p1z*(p2y - p3y)) - p2z*p3y + p2y*p3z - p1y*(-p2z + p3z))*
        log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
                (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
                (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
                (p1z - z))/
             sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
               pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
               pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar)
        )/(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
        pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
        pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),
    (-4*Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
     (Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
                p1y*(-p2z + p3z))*(p1x - x) + 
             (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
             (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        pow(thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2)) - 
     (2*Kc*(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y))*
        (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
        log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
                (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
                (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
                (p1z - z))/
             sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
               pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
               pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar)
        )/(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
        pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
        pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),
    (-4*Kc*pow(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y),2)*
        (-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        (thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))) + 
     (Kc*pow(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y),2)*
        pow(-dbar + thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + 
                p1y*(-p2z + p3z))*(p1x - x) + 
             (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*(-p1y + y) + 
             (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*(p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2))/
      ((pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
          pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
          pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))*
        pow(thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
              (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
              (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
              (p1z - z))/
           sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
             pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
             pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)),2)) - 
     (2*Kc*pow(-(p1y*(p2x - p3x)) - p2y*p3x + p2x*p3y - p1x*(-p2y + p3y),2)*
        log((thickness - ((p1z*(p2y - p3y) + p2z*p3y - p2y*p3z + p1y*(-p2z + p3z))*
                (p1x - x) + (p1z*(p2x - p3x) + p2z*p3x - p2x*p3z + p1x*(-p2z + p3z))*
                (-p1y + y) + (p1y*(p2x - p3x) + p2y*p3x - p2x*p3y + p1x*(-p2y + p3y))*
                (p1z - z))/
             sqrt(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
               pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
               pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2)))/dbar)
        )/(pow(p1y*p2x - p1x*p2y - p1y*p3x + p2y*p3x + p1x*p3y - p2x*p3y,2) + 
        pow(p1z*p2x - p1x*p2z - p1z*p3x + p2z*p3x + p1x*p3z - p2x*p3z,2) + 
        pow(p1z*p2y - p1y*p2z - p1z*p3y + p2z*p3y + p1y*p3z - p2y*p3z,2))));

  return matResult;
}

Vector3d externalContactForce::ListVec(double a1, double a2, double a3)
{
  Vector3d vecResult;

  vecResult(0) = a1;
  vecResult(1) = a2;
  vecResult(2) = a3;
  
  return vecResult;
}

Matrix3d externalContactForce::ListMat(Vector3d a1, Vector3d a2, Vector3d a3)
{
  Matrix3d matResult;

  matResult.col(0) = a1;
  matResult.col(1) = a2;
  matResult.col(2) = a3;
  
  return matResult;
}
