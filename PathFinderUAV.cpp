// PathFinderUAV.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stack>
#include <math.h>
#include <cmath>
#include <cstdlib>

using namespace std;

#define WORLD_SIZE 64

double LiftOff(double x_pos, double y_pos, double z_pos, double x_vel, double y_vel, double z_vel, double x_acc, double y_acc, double z_acc, double heading, double ang_vel, ofstream& outfile)
{
	double z_pos_prev = z_pos;
	double z_vel_prev = z_vel;
	double z_acc_prev = 0.0;
	int z_count = 0;

	double z_pos_next = 0.0;
	double z_vel_next = 0.0;
	double z_acc_next = z_acc;

	bool phase_1 = true;
	bool phase_2 = true;
	bool phase_3 = true;

	outfile << x_pos << " " << y_pos << " " << z_pos_prev << " " << x_vel << " " << y_vel << " " << z_vel_prev << " " << x_acc << " " << y_acc << " " << z_acc_prev << " " << heading << " " << ang_vel << endl;

	while (phase_1)
	{
		z_vel_next = 0.5*(z_acc_prev + z_acc_next)*0.05 + z_vel_prev;
		z_pos_next = 0.5*(z_vel_prev + z_vel_next)*0.05 + z_pos_prev;

		outfile << x_pos << " " << y_pos << " " << z_pos_next << " " << x_vel << " " << y_vel << " " << z_vel_next << " " << x_acc << " " << y_acc << " " << z_acc_next << " " << heading << " " << ang_vel << endl;

		z_pos_prev = z_pos_next;
		z_vel_prev = z_vel_next;
		z_acc_prev = z_acc_next;
		z_count++;

		if (z_vel_next > 0.35)
			z_acc_next = 0.0;

		if (z_acc_prev < 0.01)
			phase_1 = false;
	}

	z_acc_next = 0.0;
	double z_pos_covered = z_pos_prev;
	while (phase_2)
	{
		z_pos_next = 0.5*(z_vel_prev + z_vel_next)*0.05 + z_pos_prev;

		outfile << x_pos << " " << y_pos << " " << z_pos_next << " " << x_vel << " " << y_vel << " " << z_vel_next << " " << x_acc << " " << y_acc << " " << z_acc_next << " " << heading << " " << ang_vel << endl;
		z_pos_prev = z_pos_next;

		if (z_pos_next > 1.0 - z_pos_covered)
			phase_2 = false;
	}

	z_acc_next = -0.5;
	for (int i = 0; i < z_count; i++)
	{
		z_vel_next = 0.5*(z_acc_prev + z_acc_next)*0.05 + z_vel_prev;
		z_pos_next = 0.5*(z_vel_prev + z_vel_next)*0.05 + z_pos_prev;

		outfile << x_pos << " " << y_pos << " " << z_pos_next << " " << x_vel << " " << y_vel << " " << z_vel_next << " " << x_acc << " " << y_acc << " " << z_acc_next << " " << heading << " " << ang_vel << endl;

		z_pos_prev = z_pos_next;
		z_vel_prev = z_vel_next;
		z_acc_prev = z_acc_next;

		if (i == z_count - 2)
			z_acc_next = 0.0;

	}

	z_acc_next = 0.0;
	//outfile << x_pos << " " << y_pos << " " << z_p << " " << x_vel << " " << y_vel << " " << z_v << " " << x_acc << " " << y_acc << " " << z_a << " " << heading << " " << ang_vel << endl;

	return z_pos_next;
}

void Landing(double x_pos, double y_pos, double z_pos, double x_vel, double y_vel, double z_vel, double x_acc, double y_acc, double z_acc, double heading, double ang_vel, ofstream& outfile)
{
	double z_pos_prev = z_pos;
	double z_vel_prev = z_vel;
	double z_acc_prev = 0.0;
	int z_count = 0;

	double z_pos_next = 0.0;
	double z_vel_next = 0.0;
	double z_acc_next = z_acc;

	bool phase_1 = true;
	bool phase_2 = true;
	bool phase_3 = true;

	outfile << x_pos << " " << y_pos << " " << z_pos_prev << " " << x_vel << " " << y_vel << " " << z_vel_prev << " " << x_acc << " " << y_acc << " " << z_acc_prev << " " << heading << " " << ang_vel << endl;

	while (phase_1)
	{
		z_vel_next = 0.5*(z_acc_prev + z_acc_next)*0.05 + z_vel_prev;
		z_pos_next = 0.5*(z_vel_prev + z_vel_next)*0.05 + z_pos_prev;

		outfile << x_pos << " " << y_pos << " " << z_pos_next << " " << x_vel << " " << y_vel << " " << z_vel_next << " " << x_acc << " " << y_acc << " " << z_acc_next << " " << heading << " " << ang_vel << endl;

		z_pos_prev = z_pos_next;
		z_vel_prev = z_vel_next;
		z_acc_prev = z_acc_next;
		z_count++;

		if (z_vel_next < -0.35)
			z_acc_next = 0.0;

		if (z_acc_prev > -0.01)
			phase_1 = false;
	}

	z_acc_next = 0.0;
	double z_pos_covered = 1.0 - z_pos_prev;
	while (phase_2)
	{
		z_pos_next = 0.5*(z_vel_prev + z_vel_next)*0.05 + z_pos_prev;

		outfile << x_pos << " " << y_pos << " " << z_pos_next << " " << x_vel << " " << y_vel << " " << z_vel_next << " " << x_acc << " " << y_acc << " " << z_acc_next << " " << heading << " " << ang_vel << endl;
		z_pos_prev = z_pos_next;

		if (z_pos_next < z_pos_covered)
			phase_2 = false;
	}

	z_acc_next = 0.5;
/*
	while (phase_3)
	{
		z_vel_next = 0.5*(z_acc_prev + z_acc_next)*0.05 + z_vel_prev;
		z_pos_next = 0.5*(z_vel_prev + z_vel_next)*0.05 + z_pos_prev;

		outfile << x_pos << " " << y_pos << " " << z_pos_next << " " << x_vel << " " << y_vel << " " << z_vel_next << " " << x_acc << " " << y_acc << " " << z_acc_next << " " << heading << " " << ang_vel << endl;

		z_pos_prev = z_pos_next;
		z_vel_prev = z_vel_next;
		z_acc_prev = z_acc_next;
		
	}*/
	for (int i = 0; i < z_count; i++)
	{
		z_vel_next = 0.5*(z_acc_prev + z_acc_next)*0.05 + z_vel_prev;
		z_pos_next = 0.5*(z_vel_prev + z_vel_next)*0.05 + z_pos_prev;

		outfile << x_pos << " " << y_pos << " " << z_pos_next << " " << x_vel << " " << y_vel << " " << z_vel_next << " " << x_acc << " " << y_acc << " " << z_acc_next << " " << heading << " " << ang_vel << endl;

		z_pos_prev = z_pos_next;
		z_vel_prev = z_vel_next;
		z_acc_prev = z_acc_next;

		if (i == z_count - 2)
			z_acc_next = 0.0;

	}

	z_acc_next = 0.0;
	//outfile << x_pos << " " << y_pos << " " << z_p << " " << x_vel << " " << y_vel << " " << z_v << " " << x_acc << " " << y_acc << " " << z_a << " " << heading << " " << ang_vel << endl;

}

void Move(double& x_pos, double y_pos, double z_pos, double x_vel, double y_vel, double z_vel, double x_acc, double y_acc, double z_acc, double heading, double ang_vel, ofstream& outfile)
{

	double x_pos_prev = x_pos;
	double x_vel_prev = 0.0;
	double x_acc_prev = 0.0;
	int x_count = 0;

	double x_pos_next = 0.0;
	double x_vel_next = 0.0;
	double x_acc_next = x_acc;

	bool phase_1 = true;
	bool phase_2 = true;
	bool phase_3 = true;

	outfile << x_pos_prev << " " << y_pos << " " << z_pos << " " << x_vel_prev << " " << y_vel << " " << z_vel << " " << x_acc_prev << " " << y_acc << " " << z_acc << " " << heading << " " << ang_vel << endl;

	while (phase_1)
	{
		x_vel_next = 0.5*(x_acc_prev + x_acc_next)*0.05 + x_vel_prev;
		x_pos_next = 0.5*(x_vel_prev + x_vel_next)*0.05 + x_pos_prev;

		outfile << x_pos_next << " " << y_pos << " " << z_pos << " " << x_vel_next << " " << y_vel << " " << z_vel << " " << x_acc_next << " " << y_acc << " " << z_acc << " " << heading << " " << ang_vel << endl;

		x_pos_prev = x_pos_next;
		x_vel_prev = x_vel_next;
		x_acc_prev = x_acc_next;
		x_count++;

		if (x_vel_next > 0.35)
			x_acc_next = 0.0;

		if (x_acc_prev < 0.01)
			phase_1 = false;
	}

	x_acc_next = 0.0;
	double x_pos_covered = x_pos_prev-x_pos;
	while (phase_2)
	{
		x_pos_next = 0.5*(x_vel_prev + x_vel_next)*0.05 + x_pos_prev;

		outfile << x_pos_next << " " << y_pos << " " << z_pos << " " << x_vel_next << " " << y_vel << " " << z_vel << " " << x_acc_next << " " << y_acc << " " << z_acc << " " << heading << " " << ang_vel << endl;
		
		x_pos_prev = x_pos_next;

		if (x_pos_next > x_pos + 1.0 - fabs(x_pos_covered))
			phase_2 = false;
	}

	x_acc_next = -0.5;
	for (int i = 0; i < x_count; i++)
	{
		x_vel_next = 0.5*(x_acc_prev + x_acc_next)*0.05 + x_vel_prev;
		x_pos_next = 0.5*(x_vel_prev + x_vel_next)*0.05 + x_pos_prev;

		outfile << x_pos_next << " " << y_pos << " " << z_pos << " " << x_vel_next << " " << y_vel << " " << z_vel << " " << x_acc_next << " " << y_acc << " " << z_acc << " " << heading << " " << ang_vel << endl;

		x_pos_prev = x_pos_next;
		x_vel_prev = x_vel_next;
		x_acc_prev = x_acc_next;

		if (i == x_count - 2)
			x_acc_next = 0.0;

	}

	x_acc_next = 0.0;
	x_pos = x_pos_next;
	//y_pos = y_pos_next;
	//outfile << x_pos << " " << y_pos << " " << z_p << " " << x_vel << " " << y_vel << " " << z_v << " " << x_acc << " " << y_acc << " " << z_a << " " << heading << " " << ang_vel << endl;

}

void MoveDiagonal(double& x_pos, double& y_pos, double z_pos, double x_vel, double y_vel, double z_vel, double x_acc, double y_acc, double z_acc, double heading, double ang_vel, ofstream& outfile)
{

	double x_pos_prev = x_pos;
	double x_vel_prev = 0.0;
	double x_acc_prev = 0.0;
	int x_count = 0;

	double y_pos_prev = y_pos;
	double y_vel_prev = 0.0;
	double y_acc_prev = 0.0;
	int y_count = 0;

	double x_pos_next = 0.0;
	double x_vel_next = 0.0;
	double x_acc_next = x_acc;

	double y_pos_next = 0.0;
	double y_vel_next = 0.0;
	double y_acc_next = y_acc;

	bool phase_1 = true;
	bool phase_2 = true;
	bool phase_3 = true;

	outfile << x_pos_prev << " " << y_pos_prev << " " << z_pos << " " << x_vel_prev << " " << y_vel_prev << " " << z_vel << " " << x_acc_prev << " " << y_acc_prev << " " << z_acc << " " << heading << " " << ang_vel << endl;

	while (phase_1)
	{
		x_vel_next = 0.5*(x_acc_prev + x_acc_next)*0.05 + x_vel_prev;
		x_pos_next = 0.5*(x_vel_prev + x_vel_next)*0.05 + x_pos_prev;

		y_vel_next = 0.5*(y_acc_prev + y_acc_next)*0.05 + y_vel_prev;
		y_pos_next = 0.5*(y_vel_prev + y_vel_next)*0.05 + y_pos_prev;

		outfile << x_pos_next << " " << y_pos_next << " " << z_pos << " " << x_vel_next << " " << y_vel_next << " " << z_vel << " " << x_acc_next << " " << y_acc_next << " " << z_acc << " " << heading << " " << ang_vel << endl;

		x_pos_prev = x_pos_next;
		x_vel_prev = x_vel_next;
		x_acc_prev = x_acc_next;
		x_count++;

		y_pos_prev = y_pos_next;
		y_vel_prev = y_vel_next;
		y_acc_prev = y_acc_next;

		if (x_vel_next > 0.35)
		{
			x_acc_next = 0.0;
			y_acc_next = 0.0;
		}

		if (x_acc_prev < 0.01)
			phase_1 = false;
	}

	x_acc_next = 0.0;
	y_acc_next = 0.0;
	double x_pos_covered = x_pos_prev-x_pos;
	while (phase_2)
	{
		x_pos_next = 0.5*(x_vel_prev + x_vel_next)*0.05 + x_pos_prev;
		y_pos_next = 0.5*(y_vel_prev + y_vel_next)*0.05 + y_pos_prev;

		outfile << x_pos_next << " " << y_pos_next << " " << z_pos << " " << x_vel_next << " " << y_vel_next << " " << z_vel << " " << x_acc_next << " " << y_acc_next << " " << z_acc << " " << heading << " " << ang_vel << endl;

		x_pos_prev = x_pos_next;
		y_pos_prev = y_pos_next;

		if (x_pos_next > x_pos + 1.0 - fabs(x_pos_covered))
			phase_2 = false;
	}

	x_acc_next = -0.5;
	y_acc_next = -0.5;
	for (int i = 0; i < x_count; i++)
	{
		x_vel_next = 0.5*(x_acc_prev + x_acc_next)*0.05 + x_vel_prev;
		x_pos_next = 0.5*(x_vel_prev + x_vel_next)*0.05 + x_pos_prev;

		y_vel_next = 0.5*(y_acc_prev + y_acc_next)*0.05 + y_vel_prev;
		y_pos_next = 0.5*(y_vel_prev + y_vel_next)*0.05 + y_pos_prev;

		outfile << x_pos_next << " " << y_pos_next << " " << z_pos << " " << x_vel_next << " " << y_vel_next << " " << z_vel << " " << x_acc_next << " " << y_acc_next << " " << z_acc << " " << heading << " " << ang_vel << endl;

		y_pos_prev = y_pos_next;
		y_vel_prev = y_vel_next;
		y_acc_prev = y_acc_next;

		if (i == x_count - 2)
		{
			x_acc_next = 0.0;
			y_acc_next = 0.0;
		}

	}

	x_acc_next = 0.0;
	x_pos = x_pos_next;
	y_acc_next = 0.0;
	y_pos = y_pos_next;
	//y_pos = y_pos_next;
	//outfile << x_pos << " " << y_pos << " " << z_p << " " << x_vel << " " << y_vel << " " << z_v << " " << x_acc << " " << y_acc << " " << z_a << " " << heading << " " << ang_vel << endl;

}

const int x_Size = 12;
const int y_Size = 12;


struct SearchCell
{
	double x_pos, y_pos;
	int m_id;
	SearchCell *parent;
	int G;
	int H;

	SearchCell() { parent = 0; }
	SearchCell(double x, double y, SearchCell *_parent = 0) { x_pos = x, y_pos = y, parent = _parent, m_id = y*WORLD_SIZE + y, G = 0, H = 0; }

	double F() { return G + H; }
	double ManHattanDistance(SearchCell *nodeEnd)
	{
		double x = (fabs(this->x_pos - nodeEnd->x_pos));
		double y = (fabs(this->y_pos - nodeEnd->y_pos));

		return x + y;
	}
};

void FindPath(SearchCell start, SearchCell goal)
{
	bool foundGoal = false;

	stack<SearchCell> ClosedList;
	stack<SearchCell> OpenList;

	ClosedList.push(start);
	SearchCell currentCell = start;

	/*while (!foundGoal)
	{
	    SearchCell NewCell_1  = SearchCell(currentCell.x_pos-0.5, currentCell.y_pos, 0);
	    SearchCell NewCell_2  = SearchCell(currentCell.x_pos-0.5, currentCell.y_pos-0.5, 0);
	    SearchCell NewCell_3  = SearchCell(currentCell.x_pos, currentCell.y_pos-0.5, 0);
	    SearchCell NewCell_4  = SearchCell(currentCell.x_pos+0.5, currentCell.y_pos-0.5, 0);
	    SearchCell NewCell_5  = SearchCell(currentCell.x_pos+0.5, currentCell.y_pos, 0);
	    SearchCell NewCell_6  = SearchCell(currentCell.x_pos+0.5, currentCell.y_pos+0.5, 0);
	    SearchCell NewCell_7  = SearchCell(currentCell.x_pos, currentCell.y_pos+0.5, 0);
	    SearchCell NewCell_8  = SearchCell(currentCell.x_pos-0.5, currentCell.y_pos+0.5, 0);

	}
*/
}

int main()
{
		cout << "starting" << endl;

	SearchCell start = SearchCell(-1.5, 1.5, 0);
	SearchCell goal = SearchCell(2, 0, 0);
cout << "init start and end" << endl;
	FindPath(start, goal);
cout << "find path" << endl;
	vector<double> coordinates;

	ifstream infile;
	ofstream outfile;
	double z_pos = 0.0;
	double x_pos = -1.5;
	double y_pos = 1.5;

	infile.open("obstacles.txt");
	outfile.open("path.txt");
	outfile.setf(ios::fixed);
	outfile.setf(ios::showpoint);
	outfile.precision(3);
		
	cout << "opening files" << endl;
	
	if (infile.fail())
	{
		cerr << "Error opening file" << endl;
		exit(1);
	}
	
	cout << "opening files" << endl;
	
	double x, y;

	while (!infile.eof())
	{
		infile >> x >> y;

		cout.setf(ios::fixed);
		cout.setf(ios::showpoint);
		cout.precision(3);
		cout << x << " " << y << endl;

		coordinates.push_back(x);
		coordinates.push_back(y);
		cout << "pushed coordinates" << endl;

		//outfile << x << " " << y << endl;
	}
	cout << "Done opening obs.txt" << endl;
	z_pos = LiftOff(-1.5, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, outfile);
	cout << "Done running LiftOff" << endl;
	Move(x_pos, y_pos, z_pos, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, outfile);
	MoveDiagonal(x_pos, y_pos, z_pos, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, outfile);
	Landing(x_pos, y_pos, z_pos, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.0, 0.0, outfile);

	//cout << coordinates.size() << endl;

	for (int i = 0; i < coordinates.size()-1; i += 2)
	{
		cout << coordinates[i] << " " << coordinates[i + 1] << endl;
	}

	infile.close();
	outfile.close();

	system("PAUSE");

    return 0;
}


