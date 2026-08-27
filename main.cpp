/**
 * simDER
 * simDER stands for "[sim]plified [D]iscrete [E]lastic [R]ods"
 * Dec 2017
 * This code is based on previous iterations. 
 * */

//This line is for mac
//#include <GLUT/glut.h>

//This is for linux
#include <GL/glut.h>

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include "eigenIncludes.h"

// Rod and stepper are included in the world
#include "world.h"
#include "setInput.h"

std::unique_ptr<world> myWorld;
int NPTS;
ofstream outfile;

static void Key(unsigned char key, int, int)
{
  switch (key) // ESCAPE to quit
  {
	case 27:
		exit(0);
  }
}

/* Initialize OpenGL Graphics */
void initGL() 
{
	glClearColor(0.7f, 0.7f, 0.7f, 0.0f); // Set background color to black and opaque
	glClearDepth(10.0f);                   // Set background depth to farthest
	//glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	//glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	glShadeModel(GL_SMOOTH);   // Enable smooth shading
	//glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections

	glLoadIdentity();
	//gluLookAt(0.00, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	gluLookAt(0.05, 0.05, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
	glPushMatrix();

	//glMatrixMode(GL_MODELVIEW);
}

void display(void)
{
	while ( myWorld->simulationRunning() > 0)
	{
		//  Clear screen and Z-buffer
		glClear(GL_COLOR_BUFFER_BIT);


		/*

		// draw axis
		double axisLen = 1;
		glLineWidth(0.5);
		
		glBegin(GL_LINES);
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f(-axisLen, 0.0, 0.0);
			glVertex3f(axisLen, 0.0, 0.0);

			glColor3f(0.0, 1.0, 0.0);
			glVertex3f(0.0, -axisLen, 0.0);
			glVertex3f(0.0, axisLen, 0.0);

			glColor3f(0.0, 0.0, 1.0);
			glVertex3f(0.0, 0.0, -axisLen);
			glVertex3f(0.0, 0.0, axisLen);
		glEnd();

		*/
		
		//draw a line
		glColor3f(0.1, 0.1, 0.1);
		glLineWidth(3.0);
		
		glBegin(GL_LINES);
		for (int i=0; i < NPTS-1; i++)
		{
			glVertex3f( myWorld->getScaledCoordinate(4*i), myWorld->getScaledCoordinate(4*i+1), myWorld->getScaledCoordinate(4*i+2));
			glVertex3f( myWorld->getScaledCoordinate(4*(i+1)), myWorld->getScaledCoordinate(4*(i+1)+1), myWorld->getScaledCoordinate(4*(i+1)+2));
		}
		glEnd();

		glColor3f(1.0, 0.0, 0.0);
		glLineWidth(2.0);
		glBegin(GL_LINES);
		for (int i=0; i < myWorld->numTriangle(); i++)
		{
			Vector3d xCurrent1 = myWorld->getScaledCoordinateSurface(i, 0);
			Vector3d xCurrent2 = myWorld->getScaledCoordinateSurface(i, 1);
			Vector3d xCurrent3 = myWorld->getScaledCoordinateSurface(i, 2);
			
			glVertex3f( xCurrent1(0), xCurrent1(1), xCurrent1(2));
			glVertex3f( xCurrent2(0), xCurrent2(1), xCurrent2(2));

			glVertex3f( xCurrent1(0), xCurrent1(1), xCurrent1(2));
			glVertex3f( xCurrent3(0), xCurrent3(1), xCurrent3(2));

			glVertex3f( xCurrent2(0), xCurrent2(1), xCurrent2(2));
			glVertex3f( xCurrent3(0), xCurrent3(1), xCurrent3(2));
		}
		glEnd();
		
		glFlush();

		// Update step
		myWorld->updateTimeStep();
		myWorld->CoutData(outfile);
	}
	
	myWorld->CloseFile(outfile);
	exit(EXIT_SUCCESS);
}

int main(int argc,char *argv[])
{
	if (argc < 2)
	{
		cerr << "Usage: " << argv[0] << " <option-file> [-- option value ...]\n";
		return EXIT_FAILURE;
	}

	setInput inputData;
	if (inputData.LoadOptions(argv[1]) != 0)
	{
		return EXIT_FAILURE;
	}
	inputData.LoadOptions(argc,argv);
	//read input parameters from txt file and cmd

	myWorld = std::make_unique<world>(inputData);
	myWorld->setRodStepper();

	myWorld->OpenFile(outfile);

	bool render = myWorld->isRender();
	if (render) // if OpenGL visualization is on
	{
		NPTS = myWorld->numPoints();
	
		glutInit(&argc,argv);
		glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize (1000, 1000);
		glutInitWindowPosition (100, 100);
		glutCreateWindow ("simDER");
		initGL();
		glutKeyboardFunc(Key);
		glutDisplayFunc(display);
		glutMainLoop();
	}	
	else
	{
		while ( myWorld->simulationRunning() > 0)
		{
			myWorld->updateTimeStep(); // update time step
			myWorld->CoutData(outfile); // write data to file
		}
	}

	// Close (if necessary) the data file
	myWorld->CloseFile(outfile);
	
	return 0;
}
