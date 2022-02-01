
#define _USE_MATH_DEFINES
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <random>
#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Window.h"

//this struct is used to hold the state of interactive portion of this program
struct State {

	//sceneNum attribute holds the scene id that should be represented on the screen
	//set to 1 to render Sierpinski Triangles scene
	//set to 2 to render Squares and Diamonds scene
	//set to 3 to render Koch Snowflake scene
	//Set to 4 to render Dragon Curve scene
	//Set to 5 to render Hilbert Curve scene
	int sceneNum = 1;
	int subdivisionLvl = 0;//holds the current subdivision level that should be represented on the screen

	//overwrite the == comparator operator for the State struct
	bool operator==(State const& other) const {
		//compare the subdivision_lvl attribute to another State object's
		return (subdivisionLvl == other.subdivisionLvl) || (sceneNum == other.sceneNum);
	}

	//overwrite the != comparator operator for the State struct 
	bool operator!=(State const& other) const {
		//compare the subdivision_lvl attribute to another State object's
		return (subdivisionLvl != other.subdivisionLvl) || (sceneNum != other.sceneNum);
	}
};

//Callbacks Definition
class MyCallbacks : public CallbackInterface {

public:
	MyCallbacks(ShaderProgram& shader) : shader(shader) {}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (key == GLFW_KEY_UP && action == GLFW_PRESS) {//increase the subdivision number
			if (state.subdivisionLvl < 20) {
				state.subdivisionLvl++;
			}
		}
		if (key == GLFW_KEY_DOWN && action == GLFW_PRESS) {//decrease the subdivision number
			if (state.subdivisionLvl > 0) {
				state.subdivisionLvl--;
			}
		}
		if (key == GLFW_KEY_RIGHT && action == GLFW_PRESS) {//increase the scene number
			if (state.sceneNum < 5) {
				state.sceneNum++;
			}
			else if (state.sceneNum == 5) {
				state.sceneNum = 1;
			}
		}
		if (key == GLFW_KEY_LEFT && action == GLFW_PRESS) {//increase the scene number
			if (state.sceneNum > 1) {
				state.sceneNum--;
			}
			else if (state.sceneNum == 1) {
				state.sceneNum = 5;
			}
		}
	}

	State getState() {//return the State object that is stored in this MyCallbacks class
		return state;
	}

	void setState(int newSubdivisionLvl, int newSceneNum) {//update the State object that is stores in this MyCallbacksClass
		state.subdivisionLvl = newSubdivisionLvl;
		state.sceneNum = newSceneNum;
	}

private:
	State state;
	ShaderProgram& shader;
};

//recursive fuction to generate Sierpinski Triangle for Part I of the assigment
void recursiveSierpinskiTriangle(CPU_Geometry& cpuGeom, glm::vec3 ptA, glm::vec3 ptB, glm::vec3 ptC, int subdivision_lvl) {

	//declare three new points for an iteration of subdivision
	glm::vec3 ptD, ptE, ptF;

	//perform an iteration of the subdivision
	if (subdivision_lvl > 0) {
		ptD = (0.5f * ptA) + (0.5f * ptB);//D = 0.5*A+0.5*B
		ptE = (0.5f * ptB) + (0.5f * ptC);//E = 0.5*B+0.5*C
		ptF = (0.5f * ptC) + (0.5f * ptA);//F = 0.5*C+0.5*A
		recursiveSierpinskiTriangle(cpuGeom, ptD, ptB, ptE, subdivision_lvl - 1);
		recursiveSierpinskiTriangle(cpuGeom, ptA, ptD, ptF, subdivision_lvl - 1);
		recursiveSierpinskiTriangle(cpuGeom, ptF, ptE, ptC, subdivision_lvl - 1);
	}

	else {//insert vertices and colors into the cpuGeom object
		//insert the coordinates of triangle verticies into the cpuGeom object
		cpuGeom.verts.push_back(ptA);
		cpuGeom.verts.push_back(ptB);
		cpuGeom.verts.push_back(ptC);

		//generate a random integer to assign random colors in the appropriate shade for each tip of the triangle
		//the following method to generate a random number was inspired by the following source:
		//https://stackoverflow.com/questions/5008804/generating-random-integer-from-a-range
		int min = 0;
		int max = 2;
		std::random_device rd;// only used once to initialise (seed) engine
		std::mt19937 rng(rd());// random-number engine used (Mersenne-Twister in this case)
		std::uniform_int_distribution<int> uni(min, max); // guaranteed unbiased
		int random_integer = uni(rng);

		//insert the colors dependednt on triangle coordinates
		if (ptB.x > 0 && ptC.y < 0) {//if bottom right quadrant
			//use different shades of green for triangle color
			//make an array of vec3 for 3 distinct shades of green
			glm::vec3 arr[3];
			arr[0] = glm::vec3(0.0f,  0.1f, 0.0f);
			arr[1] = glm::vec3(0.0f,  0.9f, 0.0f);
			arr[2] = glm::vec3(0.01f,  0.4f, 0.1f);
			glm::vec3 randomColor = arr[random_integer];
			cpuGeom.cols.push_back(randomColor);
			cpuGeom.cols.push_back(randomColor);
			cpuGeom.cols.push_back(randomColor);
		}
		else if (ptB.y > 0) {//if top half of the window
			//use different shades of blue for triangle color
			//make an array of vec3 for 3 distinct shades of blue
			glm::vec3 arr[3];
			arr[0] = glm::vec3(0.0f, 0.0f, 0.1f);
			arr[1] = glm::vec3(0.0f, 0.0f, 0.9f);
			arr[2] = glm::vec3(0.01f, 0.01f, 0.8f);
			glm::vec3 randomColor = arr[random_integer];
			cpuGeom.cols.push_back(randomColor);
			cpuGeom.cols.push_back(randomColor);
			cpuGeom.cols.push_back(randomColor);
		}
		else {//if bottom left quadrant
			//use different shades of red for triangle colors
			//make an array of vec3 for 3 distinct shades of red
			glm::vec3 arr[3];
			arr[0] = glm::vec3(0.1f, 0.0f, 0.0f);
			arr[1] = glm::vec3(0.9f, 0.0f, 0.0f);
			arr[2] = glm::vec3(0.5f, 0.05f, 0.01f);
			glm::vec3 randomColor = arr[random_integer];
			cpuGeom.cols.push_back(randomColor);
			cpuGeom.cols.push_back(randomColor);
			cpuGeom.cols.push_back(randomColor);
		}
	}
}

//create five points that represent a Koch snowflake side for a given input pair of points
void kochLine(CPU_Geometry& cpuGeom, glm::vec3 pt1, glm::vec3 pt5, int subdivision_lvl, glm::vec3 color, State state) {

	glm::vec3 pt2, pt3, pt4;//these will be assigned the endpoints of the newly created segments

	//create an array of 6 different colors to use for each different iteration of the snowflake
	glm::vec3 arr[6];
	arr[0] = glm::vec3(0.0f, 1.0f, 0.0f);
	arr[1] = glm::vec3(0.0f, 0.0f, 1.0f);
	arr[2] = glm::vec3(1.0f, 1.0f, 1.0f);
	arr[3] = glm::vec3(1.0f, 0.2f, 1.0f);
	arr[4] = glm::vec3(1.0f, 1.0f, 0.00f);
	arr[5] = glm::vec3(0.3f, 0.0f, 0.6f);

	if (subdivision_lvl > 0) {//split the provided line into 4 segments
		//compute endpoint of the first segment
		pt2 = pt1 + (1.0f / 3.0f) * (pt5 - pt1);

		//compute endpoint of the third segment
		pt4 = pt1 + (2.0f / 3.0f) * (pt5 - pt1);

		//compute endpoint of the second segment, use a rotation matrix to rotate a vector and also shift the origin to pt2
		pt3.x = pt2.x + (pt4.x - pt2.x) * (float)cos(M_PI / 3.0f) - (pt4.y - pt2.y) * (float)sin(M_PI / 3.0f);
		pt3.y = pt2.y + (pt4.x - pt2.x) * (float)sin(M_PI / 3.0f) + (pt4.y - pt2.y) * (float)cos(M_PI / 3.0f);
		pt3.z = 0.0f;

		//recursive calls to create Koch Lines on newly generated segments
		kochLine(cpuGeom, pt1, pt2, subdivision_lvl - 1, color, state);
		kochLine(cpuGeom, pt2, pt3, subdivision_lvl - 1, arr[state.subdivisionLvl-subdivision_lvl], state);
		kochLine(cpuGeom, pt3, pt4, subdivision_lvl - 1, arr[state.subdivisionLvl-subdivision_lvl], state);
		kochLine(cpuGeom, pt4, pt5, subdivision_lvl - 1, color, state);
	}

	else {
		//base case
		//create new segment vertices and color those vertices
		cpuGeom.verts.push_back(pt1); cpuGeom.verts.push_back(pt5);	

		cpuGeom.cols.push_back(color);
        cpuGeom.cols.push_back(color);
	}
}

//kick off the recursive generation of new segments by using the initial triangle segments
void recursiveKochSnowlake(CPU_Geometry& cpuGeom, glm::vec3 ptA, glm::vec3 ptB, glm::vec3 ptC, int subdivision_lvl, State state) {

		kochLine(cpuGeom, ptA, ptB, subdivision_lvl, glm::vec3(1.0f, 0.0f, 0.0f), state);
		kochLine(cpuGeom, ptB, ptC, subdivision_lvl, glm::vec3(1.0f, 0.0f, 0.0f), state);
		kochLine(cpuGeom, ptC, ptA, subdivision_lvl, glm::vec3(1.0f, 0.0f, 0.0f), state);
}

//recursive fuction to generate Squares and Diamonds for Part II of the assigment
void recursiveSquaresAndDiamonds(CPU_Geometry& cpuGeom, glm::vec3 ptA, glm::vec3 ptB, glm::vec3 ptC, glm::vec3 ptD, int subdivision_lvl) {
		
	//create coordinates of vertices for the diamond
	glm::vec3 ptE, ptF, ptG, ptH;//declare points of the diamond		
	ptE = ptA + (ptB - ptA) * (1.0f / 2.0f);
	ptF = ptB + (ptC - ptB) * (1.0f / 2.0f);
	ptG = ptC + (ptD - ptC) * (1.0f / 2.0f);
	ptH = ptD + (ptA - ptD) * (1.0f / 2.0f);

	//perform an iteration of the subdivision
	if (subdivision_lvl > 0) {

		//compute new points to represent another square
		glm::vec3 ptA2, ptB2, ptC2, ptD2;
		ptA2 = ptH + (ptE - ptH) * (1.0f / 2.0f);
		ptB2 = ptE + (ptF - ptE) * (1.0f / 2.0f);
		ptC2 = ptF + (ptG - ptF) * (1.0f / 2.0f);
		ptD2 = ptG + (ptH - ptG) * (1.0f / 2.0f);

		//call the next iteration of square and diamond
		recursiveSquaresAndDiamonds(cpuGeom, ptA2, ptB2, ptC2, ptD2, subdivision_lvl-1);
	}

	//draw the square and diamond pair
	//place the vertices of line segments into cpuGeom
	cpuGeom.verts.push_back(ptA); cpuGeom.verts.push_back(ptB);
	cpuGeom.verts.push_back(ptB); cpuGeom.verts.push_back(ptC);
	cpuGeom.verts.push_back(ptC); cpuGeom.verts.push_back(ptD);
	cpuGeom.verts.push_back(ptD); cpuGeom.verts.push_back(ptA);
	cpuGeom.verts.push_back(ptE); cpuGeom.verts.push_back(ptF);
	cpuGeom.verts.push_back(ptF); cpuGeom.verts.push_back(ptG);
	cpuGeom.verts.push_back(ptG); cpuGeom.verts.push_back(ptH);
	cpuGeom.verts.push_back(ptH); cpuGeom.verts.push_back(ptE);
	//place the color of the vertices into cpuGeom
	cpuGeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));cpuGeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
	cpuGeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));cpuGeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
	cpuGeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));cpuGeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
	cpuGeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));cpuGeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
	cpuGeom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));cpuGeom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
	cpuGeom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));cpuGeom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
	cpuGeom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));cpuGeom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
	cpuGeom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));cpuGeom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
}

//this function recursively generates a dragon curve 
void recursiveDragonCurve(CPU_Geometry& cpuGeom, glm::vec3 ptA, glm::vec3 ptB, int subdivisionLvl) {

	glm::vec3 ptC, ptD;
	if (subdivisionLvl > 0) {
		//find the mid-point between ptA and ptB and call it ptC
		ptC = ptA + 0.5f * (ptB - ptA);

		//multiply the vector from ptC to ptA by a rotation matrix, with theta = 90 degrees
		//i.e rotate the C to B vector by 90 degrees counterclockwise.
		ptD.x = (ptA.x - ptC.x) * (float)cos(M_PI / 2.0f) - (ptA.y - ptC.y) * (float)sin(M_PI / 2.0f);
		ptD.y = (ptA.x - ptC.x) * (float)sin(M_PI / 2.0f) + (ptA.y - ptC.y) * (float)cos(M_PI / 2.0f);
		ptD.z = 0.0f;//initialize the z component to zero

		//transform the rotated vector to have origin of ptC
		ptD.x = ptD.x + ptC.x;
		ptD.y = ptD.y + ptC.y;

		//make recursive calls to subdivide the generated dragon curve
		recursiveDragonCurve(cpuGeom, ptA, ptD, subdivisionLvl - 1);
		recursiveDragonCurve(cpuGeom, ptB, ptD, subdivisionLvl - 1);
	}

	else {
		//base case
		//add vertices
		cpuGeom.verts.push_back(ptA);
		cpuGeom.verts.push_back(ptB);
		//add color of the verticies
		cpuGeom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
		cpuGeom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
	}
}

//dragon curve method from Assignment 1 walkthrough tutorial
void createDragonCurve(CPU_Geometry& cgeom, glm::vec3 startPt, glm::vec3 endPt, glm::vec3 colour, int levels) {
	//Approach in tutorial
	if (levels == 0) {
		cgeom.verts.push_back(startPt);
		cgeom.verts.push_back(endPt);

		cgeom.cols.push_back(colour);
		cgeom.cols.push_back(colour);
	}
	else {
		glm::mat4 rotMat(glm::rotate(glm::mat4(1.0f), glm::radians(-45.0f), glm::vec3(0.f, 0.f, 1.f)));
		glm::vec3 shortenedVec = (float)(1.0f / sqrt(2.0f))*(endPt - startPt);

		//do the rotation
		glm::vec4 rotatedVec4 = rotMat * glm::vec4(shortenedVec, 0.0f);
		//multi platform safe use of abs
		float x = -1.5f;
		std::cout << "My float: " << std::abs(x) << std::endl;

		//truncate to vec 3 so we can add with out start point
		glm::vec3 rotated = glm::vec3(rotatedVec4);//constructor in glm lets us switch between vec3 and vec4 easily

		glm::vec3 midPt = startPt + rotated;

		//next do a recursive call
		createDragonCurve(cgeom, startPt, midPt, colour, levels - 1);
		createDragonCurve(cgeom, endPt, midPt, colour, levels - 1);
	}
}

//this function recursively generates a hilbert curve
void recursiveHilbertCurve(CPU_Geometry& cpuGeom, float x0, float y0, float xi, float xj, float yi, float yj, int subdivisionLvl) {
	
	if (subdivisionLvl > 0) {
		//make recursive calls to add points for the growing hilbert curve
		recursiveHilbertCurve(cpuGeom, x0, y0, yi/2.0f, yj/2.0f, xi/2.0f, xj/2.0f, subdivisionLvl-1);
		recursiveHilbertCurve(cpuGeom, x0+(xi/2.0f), y0+(xj/2.0f), xi/2, xj/2, yi/2, yj/2, subdivisionLvl-1);
		recursiveHilbertCurve(cpuGeom, x0+(xi/2.0f)+(yi/2.0f), y0+(xj/2.0f)+(yj/2.0f), xi/2.0f, xj/2.0f, yi/2.0f, yj/2.0f, subdivisionLvl-1);
		recursiveHilbertCurve(cpuGeom, x0+(xi/2.0f)+yi,  y0+(xj/2.0f)+yj, -yi/2.0f, -yj/2.0f, -xi/2.0f, -xj/2, subdivisionLvl-1);
	}
	else {
		//base case
		float x = x0 + (xi + yi) / 2.0f;
		float y = y0 + (xj + yj) / 2.0f;

		glm::vec3 newPt = glm::vec3(x, y, 0.0f);//package as a vec3 type

		cpuGeom.verts.push_back(newPt);
		cpuGeom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
	}
}

int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit();
	Window window(900, 800, "CPSC 453"); // can set callbacks at construction if desired

	GLDebug::enable();

	// SHADERS
	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	// CALLBACKS
	std::shared_ptr<MyCallbacks> callbacks = std::make_shared<MyCallbacks>(shader);
	window.setCallbacks(callbacks); // can also update callbacks to new ones

	//create Geometry objects
	CPU_Geometry cpuGeom;
	GPU_Geometry gpuGeom;

	//Create and initialize our state object
	State state;
	state.sceneNum = 1;
	state.subdivisionLvl = 0;

	//Set up the initial scene
	cpuGeom.verts.clear();
	cpuGeom.cols.clear();
	recursiveSierpinskiTriangle(cpuGeom,	glm::vec3(-1.0f, -1.0f, 0.0f),
											glm::vec3(0.0f, 1.0f, 0.0f),
											glm::vec3(1.0f, -1.0f, 0.0f),
											state.subdivisionLvl);
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);

	// RENDER LOOP
	while (!window.shouldClose()) {
		glfwPollEvents();

		//check if user changed any parameters and update the cpuGeom objects accordingle
		if (callbacks->getState() != state) {
			//update the state object inside main to match the state private variable of the callback
			state = callbacks->getState();

			//clear the cpuGeom object vertices and their colors
			cpuGeom.verts.clear();
			cpuGeom.cols.clear();

			//Generate cpuGeom object that reflects the user selected scene
			if (state.sceneNum == 1) {
				//pass the coordinates of the initial triangle into the function call
				recursiveSierpinskiTriangle(cpuGeom,	glm::vec3(-1.0f, -1.0f, 0.0f),
														glm::vec3(0.0f, 1.0f, 0.0f),
														glm::vec3(1.0f, -1.0f, 0.0f),
														state.subdivisionLvl);
			}
			else if (state.sceneNum == 2) {
				//pass the coordinates of the initial square into the function call
				recursiveSquaresAndDiamonds(cpuGeom,	glm::vec3(-0.9f, -0.9f, 0.0f),
														glm::vec3(-0.9f, 0.9f, 0.0f),
														glm::vec3(0.9f, 0.9f, 0.0f),
														glm::vec3(0.9f, -0.9f, 0.0f),
														state.subdivisionLvl);
			}
			else if (state.sceneNum == 3) {
				//pass the coordinates of the initial triangle into the function call
				recursiveKochSnowlake(cpuGeom,	glm::vec3(-0.5f, -0.5f, 0.0f),
												glm::vec3( 0.0f,  0.5f, 0.0f),
												glm::vec3( 0.5f, -0.5f, 0.0f),
												state.subdivisionLvl, state);
			}
			else if (state.sceneNum == 4) {
				//pass a set of starting points into the recursive function
				recursiveDragonCurve(cpuGeom,	glm::vec3(-0.5f, 0.0f, 0.0f),
												glm::vec3(0.5f, 0.0f, 0.0f),
												state.subdivisionLvl);
			}
			else if (state.sceneNum == 5) {
				//pass a set of starting points into the recursive function
				recursiveHilbertCurve(cpuGeom, -1.0f, -1.0f, 0.0f, 2.0f, 2.0f, 0.0f, state.subdivisionLvl+1);
			}

			//bind the updated gpuGeom object
			gpuGeom.setVerts(cpuGeom.verts);
			gpuGeom.setCols(cpuGeom.cols);
		}	

		shader.use();

		glEnable(GL_FRAMEBUFFER_SRGB);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//bind and draw on the screen
		gpuGeom.bind();

		//Select the draw method to reflect the current scene
		if(state.sceneNum == 1)
			glDrawArrays(GL_TRIANGLES, 0, (GLsizei)cpuGeom.verts.size());
		else if (state.sceneNum == 2)
			glDrawArrays(GL_LINES, 0, (GLsizei)cpuGeom.verts.size());
		else if (state.sceneNum == 3)
			glDrawArrays(GL_LINES, 0, (GLsizei)cpuGeom.verts.size());
		else if (state.sceneNum == 4)
			glDrawArrays(GL_LINES, 0, (GLsizei)cpuGeom.verts.size());
		else if (state.sceneNum == 5)
			glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)cpuGeom.verts.size());


		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}

