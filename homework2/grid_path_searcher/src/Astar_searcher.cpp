#include "Astar_searcher.h"

using namespace std;
using namespace Eigen;

void AstarPathFinder::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
    
    GridNodeMap = new GridNodePtr ** [GLX_SIZE];
    for(int i = 0; i < GLX_SIZE; i++){
        GridNodeMap[i] = new GridNodePtr * [GLY_SIZE];
        for(int j = 0; j < GLY_SIZE; j++){
            GridNodeMap[i][j] = new GridNodePtr [GLZ_SIZE];
            for( int k = 0; k < GLZ_SIZE;k++){
                Vector3i tmpIdx(i,j,k);
                Vector3d pos = gridIndex2coord(tmpIdx);
                GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
            }
        }
    }
}

void AstarPathFinder::resetGrid(GridNodePtr ptr)
{
    ptr->id = 0;
    ptr->cameFrom = NULL;
    ptr->gScore = inf;
    ptr->fScore = inf;
}

void AstarPathFinder::resetUsedGrids()
{   
    for(int i=0; i < GLX_SIZE ; i++)
        for(int j=0; j < GLY_SIZE ; j++)
            for(int k=0; k < GLZ_SIZE ; k++)
                resetGrid(GridNodeMap[i][j][k]);
}

void AstarPathFinder::setObs(const double coord_x, const double coord_y, const double coord_z)
{
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      

    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

vector<Vector3d> AstarPathFinder::getVisitedNodes()
{   
    vector<Vector3d> visited_nodes;
    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++){   
                //if(GridNodeMap[i][j][k]->id != 0) // visualize all nodes in open and close list
                if(GridNodeMap[i][j][k]->id == -1)  // visualize nodes in close list only
                    visited_nodes.push_back(GridNodeMap[i][j][k]->coord);
            }

    ROS_WARN("visited_nodes size : %d", visited_nodes.size());
    return visited_nodes;
}

Vector3d AstarPathFinder::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i AstarPathFinder::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d AstarPathFinder::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

inline bool AstarPathFinder::isOccupied(const Eigen::Vector3i & index) const
{
    return isOccupied(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isFree(const Eigen::Vector3i & index) const
{
    return isFree(index(0), index(1), index(2));
}

// 越界检查+障碍物检查
inline bool AstarPathFinder::isOccupied(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return  (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
            (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == 1));// 三维转一维公式
}

inline bool AstarPathFinder::isFree(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

inline void AstarPathFinder::AstarGetSucc(GridNodePtr currentPtr, vector<GridNodePtr> & neighborPtrSets, vector<double> & edgeCostSets)
{   
    neighborPtrSets.clear();
    edgeCostSets.clear();
    /*
    *
    ✅STEP 4: finish AstarPathFinder::AstarGetSucc yourself 
    please write your code below
    26邻域搜索
    *
    *
    */
    if (currentPtr == nullptr)
        std::cout << "Error: Current pointer is null" << endl;

   for(int i=-1;i<=1;i++){
       for(int j=-1;j<=1;j++){
           for(int k=-1;k<=1;k++){
            if(i!=0 || j!=0 || k!=0){
                // int tem_x = currentPtr->index(0) + i;
                // int tem_y = currentPtr->index(1) + j;
                // int tem_z = currentPtr->index(2) + k;
                Vector3i neighborIdx = currentPtr->index + Vector3i(i,j,k);

                if(isFree(neighborIdx)){
                    if(GridNodeMap[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)]->id != -1){
                    
                    GridNodePtr neighborPtr = GridNodeMap[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];
                    //错误写法:new新建了另一个“假的节点”，跟地图无关，没人知道它是谁，也没办法更新它的状态
                    //GridNodePtr neighborPtr = new GridNode(neighborIdx, gridIndex2coord(neighborIdx));
                    neighborPtrSets.push_back(neighborPtr);
                    edgeCostSets.push_back( (neighborPtr->coord - currentPtr->coord).norm() );// 基于物理坐标（coord）距离
                    // edgeCostSets.push_back(pow((i*i + j*j + k*k), 0.5));栅格层面的坐标距离
                    }

                }
            }
           }
        }
    }
}

double AstarPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2)
{
    /* 
    choose possible heuristic function you want
    Manhattan, Euclidean, Diagonal, or 0 (Dijkstra)
    Remember tie_breaker learned in lecture, add it here ?
    *
    *
    *
    ✅STEP 1: finish the AstarPathFinder::getHeu , which is the heuristic function
    please write your code below
    *
    *
    */

    double dx = abs(node1->coord(0) - node2->coord(0));
    double dy = abs(node1->coord(1) - node2->coord(1));
    double dz = abs(node1->coord(2) - node2->coord(2));
    double manha=dx+dy+dz;
    double euc=(node1->coord - node2->coord).norm();

    double dmin = std::min({dx, dy, dz});
    double dmax = std::max({dx, dy, dz});
    double dmid = dx + dy + dz - dmin - dmax;
    double diag=dx+dy+dz-(3-sqrt(3))*dmin-(2-sqrt(2))*(dmid-dmin);

    double tie_breaker = 1.0+1.0/1000.0;
    return tie_breaker * diag;
}

void AstarPathFinder::AstarGraphSearch(Vector3d start_pt, Vector3d end_pt)
{   
    ros::Time time_1 = ros::Time::now();    

    //index of start_point and end_point
    Vector3i start_idx = coord2gridIndex(start_pt);
    Vector3i end_idx   = coord2gridIndex(end_pt);
    goalIdx = end_idx;

    //position of start_point and end_point
    start_pt = gridIndex2coord(start_idx);
    end_pt   = gridIndex2coord(end_idx);

    //Initialize the pointers of struct GridNode which represent start node and goal node
    /*
    *
    GridNodePtr startPtr = new GridNode(start_idx, start_pt);
    GridNodePtr endPtr   = new GridNode(end_idx,   end_pt);
    上面是错误写法，new出来的会导致
    1.GridNodeMap[start_idx] 的状态（如 id, gScore, cameFrom 等）没有被更新；
    2.getVisitedNodes() 里会漏掉起点；
    3.后续需要清空地图状态，清不了 new 出来的节点；
    4.以后调用 resetUsedGrids()，这个 startPtr 不会被 reset；
    5.多次调用 A* 会残留错误状态，造成隐式 Bug 或内存泄露；
    [TRACE] terminatePtr = 0x56400157bb60, GridNodeMap = 0x564000c1c570, equal? ❌ NO
    节点（起点）不来自GridNodeMap！
    *
    */
    
    GridNodePtr startPtr = GridNodeMap[start_idx.x()][start_idx.y()][start_idx.z()];
    GridNodePtr endPtr = GridNodeMap[end_idx.x()][end_idx.y()][end_idx.z()];

    //openSet is the open_list implemented through multimap in STL library
    openSet.clear();
    // currentPtr represents the node with lowest f(n) in the open_list
    GridNodePtr currentPtr  = nullptr;
    GridNodePtr neighborPtr = nullptr;

    //put start node in open set
    startPtr -> gScore = 0;
    startPtr -> fScore = getHeu(startPtr,endPtr);   
    //STEP 1: finish the AstarPathFinder::getHeu , which is the heuristic function
    startPtr -> id = 1; 
    startPtr -> coord = start_pt;
    openSet.insert( make_pair(startPtr -> fScore, startPtr) );

    // 测试：②输出两个指针地址对比
    std::cout << "[DEBUG] startPtr from : " << startPtr << std::endl;
    std::cout << "[DEBUG] GridNodeMap pointer: " << GridNodeMap[start_idx(0)][start_idx(1)][start_idx(2)] << std::endl;

    /*
    *
    STEP 2 :  some else preparatory works which should be done before while loop
    please write your code below
    *
    *
    */
    startPtr->cameFrom = nullptr;
    terminatePtr = nullptr;

    vector<GridNodePtr> neighborPtrSets;
    vector<double> edgeCostSets;

    // this is the main loop
    while ( !openSet.empty() ){
        /*
        *
        *
        ✅step 3: Remove the node with lowest cost function from open set to closed set
        please write your code below
        
        IMPORTANT NOTE!!!
        This part you should use the C++ STL: multimap, more details can be find in Homework description
        *
        *
        */

        // get the node with min f
        currentPtr=openSet.begin()->second;
        // Remove the node into closeSet
        currentPtr->id=-1;
        // Remove the node with min f
        openSet.erase(openSet.begin());

        // if the current node is the goal 
        if( currentPtr->index == goalIdx ){
            ros::Time time_2 = ros::Time::now();
            terminatePtr = currentPtr;
            ROS_WARN("[A*]{sucess}  Time in A*  is %f ms, path cost if %f m", (time_2 - time_1).toSec() * 1000.0, currentPtr->gScore * resolution );            
            return;
        }
        //get the succetion
        AstarGetSucc(currentPtr, neighborPtrSets, edgeCostSets);  //STEP 4: finish AstarPathFinder::AstarGetSucc yourself     
        
        /*
        *
        *
        ✅STEP 5:  For all unexpanded neigbors "m" of node "n", please finish this for loop
        please write your code below
        *        
        */         
        for(int i = 0; i < (int)neighborPtrSets.size(); i++){
            /*
            *
            *
            Judge if the neigbors have been expanded
            please write your code below
            
            IMPORTANT NOTE!!!
            neighborPtrSets[i]->id = -1 : expanded, equal to this node is in close set
            neighborPtrSets[i]->id = 1 : unexpanded, equal to this node is in open set
            *        
            */
            neighborPtr = neighborPtrSets[i]; // 获取邻居节点
            double gCost = currentPtr->gScore + edgeCostSets[i];
            double fCost = gCost + getHeu(neighborPtr, endPtr);

            if(neighborPtr-> id == 0){ //discover a new node, which is not in the closed set and open set
                /*
                *
                *
                ✅STEP 6:  As for a new node, do what you need do ,and then put neighbor in open set and record it
                please write your code below
                *        
                */
                neighborPtr->gScore = gCost;        // 计算当前节点到邻居的代价
                neighborPtr->fScore = fCost;        // fScore = gScore + h(n)
                neighborPtr->cameFrom = currentPtr; // 设置父节点
                openSet.insert(make_pair(neighborPtr->fScore, neighborPtr));
                neighborPtr->id = 1;
                continue;
            }
            else if(neighborPtr-> id == 1){ //this node is in open set and need to judge if it needs to update, the "0" should be deleted when you are coding
                if (gCost < neighborPtr->gScore) {
                    neighborPtr->gScore = gCost;
                    neighborPtr->fScore = fCost;
                    neighborPtr->cameFrom = currentPtr;
                    openSet.insert(make_pair(neighborPtr->fScore, neighborPtr));
                }
                /*
                *
                *
                ✅STEP 7:  As for a node in open set, update it , maintain the openset ,and then put neighbor in open set and record it
                please write your code below
                *        
                */
                continue;
            }
            else{//this node is in closed set
                /*
                *
                please write your code below
                *        
                */
                continue;
            }
        }      
    }
    //if search fails
    ros::Time time_2 = ros::Time::now();
    if((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("Time consume in Astar path finding is %f", (time_2 - time_1).toSec() );
}


vector<Vector3d> AstarPathFinder::getPath() 
{   
    vector<Vector3d> path;
    // vector<GridNodePtr> gridPath;
    /*
    *
    *
    ✅STEP 8:  trace back from the curretnt nodePtr to get all nodes along the path
    please write your code below
    *      
    */
    while(terminatePtr != nullptr){
        int x =terminatePtr->index(0);
        int y =terminatePtr->index(1);
        int z =terminatePtr->index(2);

        // 测试：③输出比较 terminatePtr 和 GridNodeMap[x][y][z] 中的节点是否是同一个
        std::cout << "[TRACE] terminatePtr = " << terminatePtr
                  << ", GridNodeMap = " << GridNodeMap[x][y][z]
                  << ", equal? " << (terminatePtr == GridNodeMap[x][y][z] ? "✅ YES" : "❌ NO") << std::endl;

        path.push_back(GridNodeMap[x][y][z]->coord);
        terminatePtr = terminatePtr->cameFrom;
    }

    // for (auto ptr: gridPath)
    //     path.push_back(ptr->coord);
        
    reverse(path.begin(),path.end());

    return path;
}