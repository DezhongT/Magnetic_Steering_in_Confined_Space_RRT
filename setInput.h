#ifndef SETINPUT_H
#define SETINPUT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "Option.h"
#include "eigenIncludes.h"

class setInput
{
public:

  typedef std::map<std::string, Option> OptionMap;
  OptionMap m_options;

  setInput();
  ~setInput();

  template <typename T>
	int AddOption(const std::string& name, const std::string& desc, const T& def);

  Option* GetOption(const std::string& name);
  bool& GetBoolOpt(const std::string& name);
  int& GetIntOpt(const std::string& name);
  double& GetScalarOpt(const std::string& name);
  Vector3d& GetVecOpt(const std::string& name);
  string& GetStringOpt(const std::string& name);

  int LoadOptions(const char* filename);
  int LoadOptions(const std::string& filename)
  {
    return LoadOptions(filename.c_str());
  }
  int LoadOptions(int argc, char** argv);

private:
	double RodLength = 1.0;
	double rodRadius = 1.0e-2;
	int numVertices = 100;
	double youngM = 1.0e7;
	double Poisson = 0.5;
	double deltaTime = 5.0e-4;
	double totalTime = 1.0;
	double tol = 1.0e-6;
	double stol = 1.0e-6;
	int maxIter = 100; // maximum number of iterations
	double density = 1000.0;
	Vector3d gVector = Vector3d::Zero();
	double viscosity = 1.0e-2;
	bool render = false;
	bool saveData = false;

	double thickness = 1.0e-1;
	Vector3d shellCenter = Vector3d::Zero();
	double shellRadius = 1.0;
	double shellMinusThickness = 0.2;
	double shellPlusThickness = 0.2;
	Vector3d cavityCenter = Vector3d::Zero();
	double cavityRadius = 2.0;
	Vector3d obstacleCenter = Vector3d(0.0, 0.5, 0.0);
	double obstacleRadius = 0.1;
	Vector3d secondObstacleCenter = Vector3d(0.0, -0.5, 0.0);
	double secondObstacleRadius = 0.1;
    double dBar = 5.0e-3;
    double stiffness = 1.0e4;
	string contactModel = "legacy";
	double tipSafeDistance = 0.0;
	int maxLineSearchIter = 20;
	double lineSearchReduction = 0.5;
	double lineSearchArmijo = 1.0e-4;
	double kktGapTolerance = 1.0e-8;
	double kktMultiplierTolerance = 1.0e-8;
	double kktComplementarityTolerance = 1.0e-10;
	int kktMaxActiveSetUpdates = 20;
	string insertionModel = "none";
	double insertionCoordinate = 0.0;
	double insertionStiffness = 1.0e3;
	Vector3d insertionAxis = Vector3d::UnitX();

	Vector3d baVector = Vector3d::Zero();
    Vector3d brVector = Vector3d::Zero();
    double muZero = 1.0;
	string magneticModel = "legacy";
	double tipDipoleMoment = 0.0;

    double scaleRendering = 1.0;
};

#include "setInput.tcc"

#endif // PROBLEMBASE_H
