#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stack>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <queue>

using namespace std;

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
	/*	
		Label each cell as XXYY
	*/
	int getID()	{
		/*
		int sum = fabs(x*1000) + fabs(y*10);
		if(x < 0)	{
			sum = sum*10 + 1;
		}
		else	{
			sum = sum*10;
		}
		if(y < 0)	{
			sum = sum*10 + 1;
		}
		else	{
			sum = sum*10;
		}
		return sum;
		*/
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


queue<cell> pathSearch(cell start, cell end, vector<int> obs)	{
	//Initialize start and end cells 
	start.g = 0;
	start.id = start.getID();
	end.id = end.getID();
	start.h = start.getH(end);
	
	//Create an open and closed list
	vector<cell> open_list;
	vector<cell> closed_list;
	queue<cell> path_q;
	bool inOpen = false;
	bool inClosed = false;
	int o_idx = 0;
	int tmp = 0;
	
	//Add start to open list
	cell curr;
	open_list.push_back(start);
	curr = start;
	
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

		//Check if current cell is the end
		/* 		if(checkGoal(curr,end))	{
			cout << "Goal reached" << endl;
			return path_q;
		} */
		
		//Search for cell with smallest F cost
		cell temp = open_list[0];
		for(int idx = 0; idx < open_list.size(); idx++)	{
			if(open_list[idx].getF() < temp.getF())	{
				temp = open_list[idx];
				tmp = idx;
			}
		}
		
		curr = temp;
		cout << curr.x << " " << curr.y << " " << curr.id << endl;
		if(checkGoal(curr,end))	{
			cout << "Goal reached" << endl;
			return path_q;
		} 
		
		//Remove chosen cell from open_list and add it to closed_list
		closed_list.push_back(curr);
		open_list.erase(open_list.begin()+tmp);
		
		//Check each neighbor of the current cell
		for(int k = 0; k < 8; k++)	{
			
			cell new_cell;
			new_cell.x = curr.x + dir_x[k]*0.5;
			new_cell.y = curr.y + dir_y[k]*0.5;
			new_cell.getID();
			new_cell.getH(end);
			new_cell.getF();
			//cout << new_cell.x << " " << new_cell.y << " " << new_cell.id << endl;
			
			//If the cells is a wall, mark them as obstacles
			if (new_cell.x <-2.45 || new_cell.x > 2.45 || new_cell.y <-2.45 || new_cell.y > 2.45)	{
				new_cell.obs = true;
				closed_list.push_back(new_cell);
				//cout << "Wall" << endl;
			}		
			
			//Check if cell is in the obs vector
			for(int idx = 0; idx < obs.size(); idx++)	{
				if(new_cell.id == obs[idx])	{
					new_cell.obs = true;
					closed_list.push_back(new_cell);
				}
			} 
			
			//Check if cell is in closed_list
			for(int idx = 0; idx < closed_list.size(); idx++)	{
				if(new_cell.id == closed_list[idx].id)	{
					inClosed = true;
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
				tmp_g = curr.g + 10;
			}
			else	{
				tmp_g = curr.g + 14;
			}
			
			//Check if cell is in open_list
			for(int idx = 0; idx < open_list.size(); idx++)	{
				if(new_cell.id == open_list[idx].id)	{
					inOpen = true;
					o_idx = idx;
					break;
				}
			}
			//If not in open_list, add to list
			if(!inOpen)	{
				new_cell.parentPtr = &curr;
				new_cell.g = tmp_g;
				new_cell.getH(end);
				new_cell.getF();
				open_list.push_back(new_cell);
			}
			//If cell is in open list, check if F score is lower than current path
			else{
				inOpen = false;
				//If F score is larger, ignore cell
				if(tmp_g >= open_list[o_idx].g)	{
					continue;
				}
				//Else, update score of cell in the open list
				else{
					open_list[o_idx].g = tmp_g;
					int f = open_list[o_idx].getF();
					open_list[o_idx].parentPtr = &curr;
					//curr = open_list[o_idx];
					path_q.push(curr);
				}
			}
		}
	}
	
	return path_q;
	
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
	end.y = 0.0;
	queue<cell> visited_path;
	cout << "starting path search" << endl;
	visited_path = pathSearch(start,end,obs);
	cout << "done searching" << endl;
	
}
