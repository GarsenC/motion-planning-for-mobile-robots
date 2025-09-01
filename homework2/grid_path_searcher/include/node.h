#ifndef _NODE_H_
#define _NODE_H_

#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include "backward.hpp"

#define inf 1>>20
struct GridNode;
typedef GridNode* GridNodePtr;

struct GridNode
{    
    // 数据成员
    int id;        // 1--> open set, -1 --> closed set, 0 --> unvisited
    Eigen::Vector3d coord; // 三维实数向量 world 3D position
    Eigen::Vector3i dir;   // 整数向量 direction of expanding
    Eigen::Vector3i index; // 整数索引 grid 3D position
	
    double gScore, fScore;
    GridNodePtr cameFrom; // mark the father node
    std::multimap<double, GridNodePtr>::iterator nodeMapIt;

    // 构造函数（带参数）
    GridNode(Eigen::Vector3i _index, Eigen::Vector3d _coord){  
		// 初始化节点各个属性
    id = 0;
		index = _index;
		coord = _coord;
		dir   = Eigen::Vector3i::Zero();// 初始化为 (0,0,0)

		gScore = inf;
		fScore = inf;// 初始化为无穷大（inf，表示初始未走过）
		cameFrom = NULL;
    }

    GridNode(){};
    ~GridNode(){};
};


#endif
