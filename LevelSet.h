#ifndef LEVELSET_H
#define LEVELSET_H

#include <math.h>

#include "Grid.h"

#include "util.h"

class LevelSet : public Grid<float> {
public:
	LevelSet(int size_x, int size_y, float length_x, float length_y);
	~LevelSet();

	float distance(int from_i, int from_j, int to_i, int to_j);

	float computeUpwindGradientX(int i, int j, float vel_x);
	float computeUpwindGradientY(int i, int j, float vel_y);
};

#endif