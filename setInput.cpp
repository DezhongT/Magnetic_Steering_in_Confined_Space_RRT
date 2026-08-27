#include "setInput.h"

#include <stdexcept>

using namespace std;

setInput::setInput()
{
	AddOption("render", "visualization", render);
	AddOption("saveData", "should results be saved", saveData);
	AddOption("rodRadius", "Radius of Rod", rodRadius);
  AddOption("RodLength", "RodLength", RodLength);
  AddOption("numVertices", "numVertices", numVertices);
	AddOption("youngM", "Young's Modulus", youngM);
	AddOption("Poisson", "Poisson Ratio", Poisson);
	AddOption("deltaTime", "Time Step Length", deltaTime);
	AddOption("totalTime", "Total Running Time", totalTime);
	AddOption("tol", "Tolerance of Newton Method", tol);
	AddOption("stol", "Ratio between initial and final error", stol);
	AddOption("maxIter", "Maximum Running Times of Each Stepper", maxIter);
	AddOption("density", "Density of the Rod", density);
	AddOption("viscosity", "Viscous Froce", viscosity);
	AddOption("gVector", "Gravity", gVector);
  AddOption("baVector", "baVector", baVector);
  AddOption("brVector", "brVector", brVector);
  AddOption("muZero", "muZero", muZero);
  AddOption("magneticModel", "legacy or axial_tip", magneticModel);
  AddOption("tipDipoleMoment", "axial tip dipole magnitude", tipDipoleMoment);
  AddOption("scaleRendering", "scaleRendering", scaleRendering);
  AddOption("thickness", "thickness", thickness);
  AddOption("shellCenter", "spherical shell center", shellCenter);
  AddOption("shellRadius", "spherical shell reference radius", shellRadius);
  AddOption("shellMinusThickness", "spherical shell inward thickness", shellMinusThickness);
  AddOption("shellPlusThickness", "spherical shell outward thickness", shellPlusThickness);
  AddOption("cavityCenter", "outer spherical cavity center", cavityCenter);
  AddOption("cavityRadius", "outer spherical cavity radius", cavityRadius);
  AddOption("obstacleCenter", "excluded spherical obstacle center", obstacleCenter);
  AddOption("obstacleRadius", "excluded spherical obstacle radius", obstacleRadius);
  AddOption("secondObstacleCenter", "second excluded spherical obstacle center", secondObstacleCenter);
  AddOption("secondObstacleRadius", "second excluded spherical obstacle radius", secondObstacleRadius);
  AddOption("dBar", "dBar", dBar);
  AddOption("stiffness", "stiffness", stiffness);
  AddOption("contactModel", "legacy, planar_barrier, spherical_shell_barrier, spherical_obstacle_barrier, double_spherical_obstacle_barrier, or none", contactModel);
  AddOption("tipSafeDistance", "minimum permitted tip clearance", tipSafeDistance);
  AddOption("maxLineSearchIter", "maximum static-Newton backtracks", maxLineSearchIter);
  AddOption("lineSearchReduction", "static-Newton step reduction", lineSearchReduction);
  AddOption("lineSearchArmijo", "static-Newton Armijo coefficient", lineSearchArmijo);
  AddOption("kktGapTolerance", "planar KKT gap tolerance", kktGapTolerance);
  AddOption("kktMultiplierTolerance", "planar KKT multiplier tolerance", kktMultiplierTolerance);
  AddOption("kktComplementarityTolerance", "planar KKT complementarity tolerance", kktComplementarityTolerance);
  AddOption("kktMaxActiveSetUpdates", "maximum planar KKT active-set changes", kktMaxActiveSetUpdates);
  AddOption("insertionModel", "none or proximal_guide", insertionModel);
  AddOption("insertionCoordinate", "initial nonnegative insertion coordinate", insertionCoordinate);
  AddOption("insertionStiffness", "proximal guide axial stiffness", insertionStiffness);
  AddOption("insertionAxis", "proximal guide insertion axis", insertionAxis);
}

setInput::~setInput()
{
	;
}

Option* setInput::GetOption(const string& name)
{
  const OptionMap::iterator option = m_options.find(name);
  if (option == m_options.end())
  {
    throw invalid_argument("Option " + name + " does not exist");
  }
  return &(option->second);
}

bool& setInput::GetBoolOpt(const string& name)
{
  return GetOption(name)->b;
}

int& setInput::GetIntOpt(const string& name)
{
  return GetOption(name)->i;
}

double& setInput::GetScalarOpt(const string& name)
{
  return GetOption(name)->r;
}

Vector3d& setInput::GetVecOpt(const string& name)
{
  return GetOption(name)->v;
}

string& setInput::GetStringOpt(const string& name)
{
  return GetOption(name)->s;
}

int setInput::LoadOptions(const char* filename)
{
  ifstream input(filename);
  if (!input.is_open()) 
  {
    cerr << "ERROR: File " << filename << " not found" << endl;
    return -1;
  }

  string line, option;
  istringstream sIn;
  string tmp;
  for (getline(input, line); !input.eof(); getline(input, line)) 
  {
    sIn.clear();
    option.clear();
    sIn.str(line);
    sIn >> option;
    if (option.size() == 0 || option.c_str()[0] == '#') continue;
    OptionMap::iterator itr;
    itr = m_options.find(option);
    if (itr == m_options.end()) 
    {
      cerr << "Invalid option: " << option << endl;
      continue;
    }
    if (itr->second.type == Option::BOOL) 
    {
      sIn >> tmp;
      if (tmp == "true" || tmp == "1") itr->second.b = true;
      else if (tmp == "false" || tmp == "0") itr->second.b = false;
    } 
    else if (itr->second.type == Option::INT) 
    {
      sIn >> itr->second.i;
    } 
    else if (itr->second.type == Option::DOUBLE) 
    {
      sIn >> itr->second.r;
    } 
    else if (itr->second.type == Option::VEC) 
    {
      Vector3d& v = itr->second.v;
      sIn >> v[0];
      sIn >> v[1];
      sIn >> v[2];
    } 
    else if (itr->second.type == Option::STRING) 
    {
      sIn >> itr->second.s;
    } else 
    {
      cerr << "Invalid option type" << endl;
    }
  }
  input.close();

  return 0;
}

int setInput::LoadOptions(int argc, char** argv)
{
  string option, tmp;
  int start = 0;
  while (start < argc && string(argv[start]) != "--") ++start;
  for (int i = start + 1; i < argc; ++i) 
  {
    option = argv[i];
    OptionMap::iterator itr;
    itr = m_options.find(option);
    if (itr == m_options.end()) 
    {
      cerr << "Invalid option on command line: " << option << endl;
      continue;
    }
    if (i == argc - 1) 
    {
      cerr << "Too few arguments on command line" << endl;
      break;
    }
    if (itr->second.type == Option::BOOL) 
    {
      tmp = argv[i+1]; ++i;
      if (tmp == "true" || tmp == "1") itr->second.b = true;
      if (tmp == "false" || tmp == "0") itr->second.b = false;
    } 
    else if (itr->second.type == Option::INT) 
    {
      itr->second.i = atoi(argv[i+1]); ++i;
    } 
    else if (itr->second.type == Option::DOUBLE) 
    {
      itr->second.r = atof(argv[i+1]); ++i;
    } 
    else if (itr->second.type == Option::VEC) 
    {
      if (i >= argc - 3) 
      {
        cerr << "Too few arguments on command line" << endl;
        break;
      }
      Vector3d& v = itr->second.v;
      v[0] = atof(argv[i+1]); ++i;
      v[1] = atof(argv[i+1]); ++i;
      v[2] = atof(argv[i+1]); ++i;
    } 
    else if (itr->second.type == Option::STRING) 
    {
      itr->second.s = argv[i+1]; ++i;
    } 
    else 
    {
      //cerr << "Invalid option type" << endl;
    }
  }
  return 0;
}
