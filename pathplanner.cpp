#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stack>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <queue>
#include <stack>

using namespace std;

double LiftOff(double x_pos, double y_pos, double z_pos, double x_vel, double y_vel, double z_vel, double x_acc, double y_acc, double z_acc, double heading, double ang_vel, ofstream& outfile)	{
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

	for (int i = 0; i < 200; i++)
	{
		outfile << x_pos << " " << y_pos << " " << z_pos_next << " " << x_vel << " " << y_vel << " " << z_vel_next << " " << x_acc << " " << y_acc << " " << z_acc_next << " " << heading << " " << ang_vel << endl;
	}

	z_acc_next = 0.0;
	//outfile << x_pos << " " << y_pos << " " << z_p << " " << x_vel << " " << y_vel << " " << z_v << " " << x_acc << " " << y_acc << " " << z_a << " " << heading << " " << ang_vel << endl;

	return z_pos_next;
}

void Landing(double x_pos, double y_pos, double z_pos, double x_vel, double y_vel, double z_vel, double x_acc, double y_acc, double z_acc, double heading, double ang_vel, ofstream& outfile)	{
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



void Move(double x_pos_1, double y_pos_1, double x_pos_2, double y_pos_2, double z_pos, ofstream& outfile)
{
	double diff_x = x_pos_2 - x_pos_1;
	double diff_y = y_pos_2 - y_pos_1;

	double x_target = fabs(diff_x);
	double y_target = fabs(diff_y);

	double x_pos_prev = x_pos_1;
	double x_vel_prev = 0.0;
	double x_acc_prev = 0.0;
	//int x_count = 0;

	double x_pos_next = 0.0;
	double x_vel_next = 0.0;
	double x_acc_next = 0.0;
	double x_acc_initial = 0.0;

	double y_pos_prev = y_pos_1;
	double y_vel_prev = 0.0;
	double y_acc_prev = 0.0;
	//int y_count = 0;

	int count = 0;
	double y_pos_next = 0.0;
	double y_vel_next = 0.0;
	double y_acc_next = 0.0;
	double y_acc_initial = 0.0;

	double z_vel = 0.0;
	double z_acc = 0.0;
	double heading = 0.0;
	double ang_vel = 0.0;
	bool move_in_x = false;
	bool move_in_y = false;

	if (diff_x > 0.45 && (diff_y < 0.05 && diff_y > -0.05))
	{
		x_acc_initial = 0.5;
		y_acc_initial = 0.0;
		move_in_x = true;
		move_in_y = false;
	}
	else if (diff_x > 0.45 && diff_y > 0.45)
	{
		x_acc_initial = 0.5;
		y_acc_initial = 0.5;
		move_in_x = true;
		move_in_y = true;
	}
	else if ((diff_x < 0.05 && diff_x > -0.05) && diff_y > 0.45)
	{
		x_acc_initial = 0.0;
		y_acc_initial = 0.5;
		move_in_x = false;
		move_in_y = true;
	}
	else if (diff_x < -0.45 && diff_y > 0.45)
	{
		x_acc_initial = -0.5;
		y_acc_initial = 0.5;
		move_in_x = true;
		move_in_y = true;
	}
	else if (diff_x < -0.45 && (diff_y < 0.05 && diff_y > -0.05))
	{
		x_acc_initial = -0.5;
		y_acc_initial = 0.0;
		move_in_x = true;
		move_in_y = false;
	}
	else if (diff_x < -0.45 && diff_y < -0.45)
	{
		x_acc_initial = -0.5;
		y_acc_initial = -0.5;
		move_in_x = true;
		move_in_y = true;
	}
	else if ((diff_x < 0.05 && diff_x > -0.05) && diff_y < -0.45)
	{
		x_acc_initial = 0.0;
		y_acc_initial = -0.5;
		move_in_x = false;
		move_in_y = true;
	}
	else if (diff_x > 0.45 && diff_y < -0.45)
	{
		x_acc_initial = 0.5;
		y_acc_initial = -0.5;
		move_in_x = true;
		move_in_y = true;
	}

	x_acc_next = x_acc_initial;
	y_acc_next = y_acc_initial;

	bool phase_1 = true;
	bool phase_2 = true;
	bool phase_3 = true;
	
	while (phase_1)
	{
		//x-direction
		x_vel_next = 0.5*(x_acc_prev + x_acc_next)*0.05 + x_vel_prev;
		x_pos_next = 0.5*(x_vel_prev + x_vel_next)*0.05 + x_pos_prev;
		
		//y-direction
		y_vel_next = 0.5*(y_acc_prev + y_acc_next)*0.05 + y_vel_prev;
		y_pos_next = 0.5*(y_vel_prev + y_vel_next)*0.05 + y_pos_prev;
		
		outfile << x_pos_next << " " << y_pos_next << " " << z_pos << " " << x_vel_next << " " << y_vel_next << " " << z_vel << " " << x_acc_next << " " << y_acc_next << " " << z_acc << " " << heading << " " << ang_vel << endl;

		x_pos_prev = x_pos_next;
		x_vel_prev = x_vel_next;
		x_acc_prev = x_acc_next;
		//x_count++;
		count++;

		y_pos_prev = y_pos_next;
		y_vel_prev = y_vel_next;
		y_acc_prev = y_acc_next;
		//y_count++;
		if (move_in_x)
		{
			if (fabs(x_vel_next) > 0.35)
			{
				x_acc_next = 0.0;
				y_acc_next = 0.0;
			}
		}
		else if (move_in_y)
		{
			if (fabs(y_vel_next) > 0.35)
			{
				x_acc_next = 0.0;
				y_acc_next = 0.0;
			}
		}

		if (move_in_x)
		{
			if (fabs(x_acc_prev) < 0.01)
				phase_1 = false;
		}
		else if (move_in_y)
		{
			if (fabs(y_acc_prev) < 0.01)
				phase_1 = false;
		}
	}

	x_acc_next = 0.0;
	y_acc_next = 0.0;
	double x_pos_covered = x_pos_prev - x_pos_1;
	double y_pos_covered = y_pos_prev - y_pos_1;
	double x_pos_curr = 0;
	double y_pos_curr = 0;
	
	while (phase_2)
	{
		
		x_pos_next = 0.5*(x_vel_prev + x_vel_next)*0.05 + x_pos_prev;
		y_pos_next = 0.5*(y_vel_prev + y_vel_next)*0.05 + y_pos_prev;

		outfile << x_pos_next << " " << y_pos_next << " " << z_pos << " " << x_vel_next << " " << y_vel_next << " " << z_vel << " " << x_acc_next << " " << y_acc_next << " " << z_acc << " " << heading << " " << ang_vel << endl;
		
		x_pos_prev = x_pos_next;
		y_pos_prev = y_pos_next;
		x_pos_curr = fabs(x_pos_next - x_pos_1);
		y_pos_curr = fabs(y_pos_next - y_pos_1);
		
		if (x_target > 0.1)
		{
			if (x_pos_curr > x_target - fabs(x_pos_covered))
				phase_2 = false;
		}
		else if (y_target > 0.1)
		{
			if (y_pos_curr > y_target - fabs(y_pos_covered))
				phase_2 = false;
		}
	}

	x_acc_next = -x_acc_initial;
	y_acc_next = -y_acc_initial;

	for (int i = 0; i < count; i++)
	{
		x_vel_next = 0.5*(x_acc_prev + x_acc_next)*0.05 + x_vel_prev;
		x_pos_next = 0.5*(x_vel_prev + x_vel_next)*0.05 + x_pos_prev;
		y_vel_next = 0.5*(y_acc_prev + y_acc_next)*0.05 + y_vel_prev;
		y_pos_next = 0.5*(y_vel_prev + y_vel_next)*0.05 + y_pos_prev;

		outfile << x_pos_next << " " << y_pos_next << " " << z_pos << " " << x_vel_next << " " << y_vel_next << " " << z_vel << " " << x_acc_next << " " << y_acc_next << " " << z_acc << " " << heading << " " << ang_vel << endl;

		x_pos_prev = x_pos_next;
		x_vel_prev = x_vel_next;
		x_acc_prev = x_acc_next;

		y_pos_prev = y_pos_next;
		y_vel_prev = y_vel_next;
		y_acc_prev = y_acc_next;

		if (i == count - 2)	{
			x_acc_next = 0.0;
			y_acc_next = 0.0;
		}

	}
}

struct cell	{
	//Coordinates
	double x;
	double y;
	double id;
	
	//Cost 
	int g;	
	int h;
	int f; 
	
	//Obstacle?
	bool obs = false; 

	//Parent celll
	cell* parentPtr = NULL;

	//Label each cell as XXYY
	int getID()	{
		id = (int) (x*1000) + (int) (y*10);
		return id;
	}
	
	//Manhattan Distance
	int getH(cell end)	{
		int x_h = (int)(fabs(this->x - end.x)*20);
		int y_h = (int)(fabs(this->y - end.y)*20);
		h = x_h + y_h;
		return h;
	}
	
	int getF()	{
		f = g+h;
		return f;
	}
};

bool checkGoal(cell curr, cell end)	{
	if(curr.id == end.id)	{
		return true;
	}
	else	{
		return false;
	}
	
}

stack<cell> path(cell end, cell start)	{
	cell current = end;
	cell *tmp;
	tmp=&current;
	stack<cell> visitedPath;

	while (current.parentPtr != NULL)	{
		//cout << "ParentCell is not NULL " << current.parentPtr  << " " << &current << endl;
		visitedPath.push(current);
		tmp = current.parentPtr;
		current = *tmp;
	}
	
	visitedPath.push(start);
	return visitedPath;
}


stack<cell> pathSearch(cell start, cell end, vector<int> obs)	{
	//Initialize start and end cells 
	start.g = 0;
	start.id = start.getID();
	end.id = end.getID();
	start.h = start.getH(end);
	stack<cell> empty_stack;
	
	//Create an open and closed list
	vector<cell> open_list;
	open_list.clear();
	vector<cell> closed_list;
	closed_list.clear();
	queue<cell> path_q;
	bool inOpen = false;
	bool inClosed = false;
	int o_idx = 0;
	int tmp = 0;
	
	//Add start to open list
	open_list.push_back(start);
	cell tempo=start;
	
	
	/*
		For checking adjacent cells	
		8 Directions are:
		(-1, 1)   (0, 1)  (1, 1)
		(-1, 0)   (0 ,0)  (1, 0)
		(-1, -1) (0, -1) (1,-1)
	*/
	int dir_x[8] = {1,1,0,-1,-1,-1,0,1};
	int dir_y[8] = {0,-1,-1,-1,0,1,1,1};
	
	while(!open_list.empty())	{
		cell* curr= new cell;
		//cout << "------------------------------" << endl;
		
		//Search for cell with smallest F cost
		cell temp = open_list[0];
		tmp = 0;
		for(int idx = 0; idx < open_list.size(); idx++)	{
			if(open_list[idx].getF() < temp.getF())	{
				temp = open_list[idx];
				tmp = idx;
			}
		}
		
		(*curr) = temp;
		//cout << (*curr).x << " " << (*curr).y << " " << (*curr).id << " " << &curr  << " " << (*curr).parentPtr << endl;
		if(checkGoal((*curr),end))	{
			cout << "Goal reached" << endl;
			return path((*curr), start);
		} 
		
		//Remove chosen cell from open_list and add it to closed_list
		closed_list.push_back((*curr));
		open_list.erase(open_list.begin()+tmp);
		//cout << open_list.size() << " " << closed_list.size() << endl;
		//Check each neighbor of the current cell
		for(int k = 0; k < 8; k++)	{
			
			cell* new_cell= new cell;
			(*new_cell).x = (*curr).x + dir_x[k]*0.5;
			(*new_cell).y = (*curr).y + dir_y[k]*0.5;
			(*new_cell).getID();
			(*new_cell).getH(end);
			(*new_cell).getF();
			//cout << (*new_cell).x << " " << (*new_cell).y << " " << (*new_cell).id << endl;
			
			//If the cells is a wall, mark them as obstacles
			if ((*new_cell).x <-2.45 || (*new_cell).x > 2.45 || (*new_cell).y <-2.45 || (*new_cell).y > 2.45)	{
				(*new_cell).obs = true;
				closed_list.push_back((*new_cell));
				//cout << "Wall" << endl;
			}		
			
			//Check if cell is in the obs vector
			for(int idx = 0; idx < obs.size(); idx++)	{
				if((*new_cell).id == obs[idx])	{
					(*new_cell).obs = true;
					closed_list.push_back((*new_cell));
					//cout << "cell is obstacle" << endl;
				}
			} 
			
			//Check if cell is in closed_list
			for(int idx = 0; idx < closed_list.size(); idx++)	{
				if((*new_cell).id == closed_list[idx].id)	{
					inClosed = true;
					//cout << "cell is in closed list" << endl;
					break;
				}
			}
			
			//Skip cell if in closed_list
			if(inClosed)	{
				inClosed = false;
				continue;
			}
			
			//Get tentative g score
			int tmp_g;
			if(k%2 == 0)	{
				tmp_g = (*curr).g + 10;
			}
			else	{
				tmp_g = (*curr).g + 14;
			}
			
			//Check if cell is in open_list
			for(int idx = 0; idx < open_list.size(); idx++)	{
				if((*new_cell).id == open_list[idx].id)	{
					inOpen = true;
					o_idx = idx;
					//cout << "cell is in open list" << endl;
					break;
				}
			}
			//If not in open_list, add to list
			if(!inOpen)	{
				(*new_cell).parentPtr = curr;
				(*new_cell).g = tmp_g;
				(*new_cell).getH(end);
				(*new_cell).getF();
				open_list.push_back((*new_cell));
				//cout << "cell is not in open list" << endl;
			}
			//If cell is in open list, check if F score is lower than current path
			else{
				inOpen = false;
				//If F score is larger, ignore cell
				if(tmp_g >= open_list[o_idx].g)	{
					//cout << "cost is higher" << endl;
					continue;
				}
				//Else, update score of cell in the open list
				else{
					//cout << "updating scores" << endl;
					open_list[o_idx].g = tmp_g;
					int f = open_list[o_idx].getF();
					open_list[o_idx].parentPtr = curr;
					//curr = open_list[o_idx];
					//path_q.push(curr);
				}
			}
		}
	}
	
	return empty_stack;
	
}

float getGradient(cell A, cell B)
{
	float grad;
if((A.y-B.y)==0.0)
	{
		grad=0.0;	
	}
else if ((A.x-B.x)==0.0)
		{
			grad =100.0;	
		}else
			{
				grad=(A.y-B.y)/(A.x-B.x);
			}
return grad;
}

bool Combine(cell A, cell B, cell C)
{
	bool combine=false;
	float gradient;
	float gradient2;

	gradient=getGradient(A,B);
	gradient2=getGradient(B,C);
	if(gradient2==gradient)
	{
		combine=true; 
	}
 return combine;
}

int main()	{
	
	vector<int> obs;
	//For marking obstacles
	int odir_x[4] = {1,0,-1,0};
	int odir_y[4] = {0,-1,0,1};
	double obs_x;
	double obs_y;
	
	//Open obs.txt
	ifstream infile;
	infile.open("obs.txt");
	
	//Go through obs.txt and add obstacles to the obstacles vector
	double x, y;
	while(!infile.eof())	{
		infile >> x >> y;
		obs.push_back((int) (x*1000) + (int) (y*10));
		for(int idx = 0; idx < 4; idx++)	{
			obs_x = x + odir_x[idx]*0.5;
			obs_y = y + odir_y[idx]*0.5;
			obs.push_back((int) (obs_x*1000) + (int) (obs_y*10));
		}
		
	}
	
	for(int idx = 0; idx < obs.size(); idx++)	{
		cout << obs[idx] << endl;
	}
	
	//A* Algorithm
	cell start, end;
	start.x = -1.5;
	start.y = 1.5;
	end.x = 2.0;
	end.y = 0;

	stack<cell> visited_path;
	cout << "starting path search" << endl;
	visited_path = pathSearch(start,end,obs);
	cout << "done searching" << endl;
	
	//Generate drone path
	ofstream outfile;
	outfile.open("path.txt");
	cout.setf(ios::showpoint);
	cout.precision(3);
	outfile.setf(ios::showpoint);
	outfile.setf(ios::fixed);
	outfile.precision(3);
	int z_pos;
	int ind=0;
	//Take Off
	z_pos = LiftOff(start.x, start.y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, outfile);
	cout<<"printing coordinates"<<endl;
	vector<cell> PathSimplified;
	while(!visited_path.empty()){
		
		PathSimplified.push_back(visited_path.top()); 
		//cout << visited_path.top().x << " " << visited_path.top().y << endl;
		visited_path.pop(); 
		}

	for (int idx = 0; idx < PathSimplified.size()-1; idx++)
	{
		int idx_a = idx;
		int idx_b = idx + 1;
		for(int i = idx_a; i < PathSimplified.size()-2; i++)	
		{
			if(Combine(PathSimplified[i],PathSimplified[i+1],PathSimplified[i+2]))
			{
				idx_b = i+2;
			}
			else
				break;
		}
		
		Move(PathSimplified[idx_a].x, PathSimplified[idx_a].y, PathSimplified[idx_b].x, PathSimplified[idx_b].y, z_pos, outfile);
		cout << PathSimplified[idx_a].x << " " << PathSimplified[idx_a].y << endl;
		idx=idx_b-1;
	}	


	cell start_2, end_2;
	start_2.x = 2.0;
	start_2.y = 0.0;
	end_2.x = -1.5;
	end_2.y = -1.5;

	stack<cell>visited_path_2;
	cout << "starting path search" << endl;
	visited_path_2 = pathSearch(start_2,end_2,obs);
	cout << "done searching" << endl;

	vector<cell> PathSimplified_2;
	while(!visited_path_2.empty()){
		
		PathSimplified_2.push_back(visited_path_2.top()); 
		//cout << visited_path_2.top().x << " " << visited_path_2.top().y << endl;
		visited_path_2.pop(); 

	}

	for (int idx = 0; idx < PathSimplified_2.size()-1; idx++)
	{

		int idx_a = idx;
		int idx_b = idx + 1;
		for(int i = idx_a; i < PathSimplified_2.size()-2; i++)	
		{
			if(Combine(PathSimplified_2[i],PathSimplified_2[i+1],PathSimplified_2[i+2]))
			{
				idx_b = i+2;
			}
			else
				break;
		}
		
		Move(PathSimplified_2[idx_a].x, PathSimplified_2[idx_a].y, PathSimplified_2[idx_b].x, PathSimplified_2[idx_b].y, z_pos, outfile);
		cout << PathSimplified_2[idx_a].x << " " << PathSimplified_2[idx_a].y << endl;
		idx=idx_b-1;
	}	
	
	Landing(end_2.x, end_2.y, z_pos, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.0, 0.0, outfile);
}




