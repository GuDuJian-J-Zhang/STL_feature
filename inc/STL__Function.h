#ifndef STL_FUNCTION_H
#define STL_FUNCTION_H

#include <vector>
#include "STL__Kernel.h"
#define  PI 3.1415926535898
using namespace std;

ostream & operator<<(ostream &os,const Point_3 &P);
ostream & operator<<(ostream &os,const Face_Tri &Tri);
ostream & operator<<(ostream &os,const Edge &e);
bool operator<(Point_3 P1, Point_3 P2);

//STL文件读取函数
int Read_STL(vector<Point_3> &Point_List, vector<HalfEdge> &Half_Edge_List, vector<Face_Tri> &Tri_List);//读取STL文件,并建立三角面链表Tri_List和模型节点树_rb_tree_node


/*
//查找函数
int Search_Edge_Of_Node(Point_3& P, Edge& e);//判断边e的反向边是否在点p的入射半边表中
int Search_Edge(vector<Edge> &Edge_List, Edge e);//查找是否有相同边*/

/*
//红黑树算法相关函数
int Compare_Node(const Point_3 p1, const Point_3 p2);//比较两点p1、p2的大小;  p1 > p2 返回1, p1 = p2 返回0, p1 < p2 返回-1*/


//拓扑重建函数
Model Build_ToPo_Of_Model(vector<Point_3> &Point_List, vector<HalfEdge> &H_E_L, vector<Face_Tri> &Tri_List);//创建模型的拓扑结构
void Build_Edge_List(vector<Point_3> &Point_List, vector<HalfEdge> &H_E_L, vector<Edge> &Edge_List);//遍历各点,构建模型边表
void Build_TriandEdge_List_Of_Face(vector<Face_Tri>& Tri_List, vector<Edge> Edge_List);//遍历三角形，建立三角形的邻接面表
/*
void Del_FaceListOfNode(int tri_num,Point& p);//将面tri从点p的邻接面表中删除
void Add_FaceListOfNode(Triangle tri,Point* p);//将面tri加入点p的邻接面表
void Set_Triangle_List_OfNode(ListOfTriangle* Tri_List);//设置各节点的面表
void ReNew_EdgeListOfNode(Edge tri,Point* p);//更新节点的对边表
void Set_Edge_List_OfNode(ListOfTriangle* Tri_List);//设置各节点的入射半边表*/

//网格重剖分
vector<int> Find_Tri_List_OfEdge(vector<Edge>::iterator e, vector<HalfEdge>& H_L);//查找并返回边e的邻接面编号
void Get_Feature_Line(Model& MyModel);//提取模型的特征边和特征点
vector<Sub_Set> Surface_div(Model& MyModel);//模型子域分解
int Search_Face(vector<Face_Tri> F_L);//查找有无未分组的面片，如果有 返回其编号  否则返回0
/*void Mesh(Model& MyModel, double MySize);//网格剖分
int Greater_Length(Edge e1, Edge e2);//边长比较函数  e1 > e2: 返回1  否则 返回0

//模型修复函数
void Repair_Overlay(Model* model1);//修复重叠错误
void Merge_Node(Model* model1, Point P1, Point p2);//合并顶点,将两个距离 < 参数值的顶点合并
void ReSet_Edge_List_OfNode(Point& P, Point New_Point, Point P1, Point P2);//更新点P的入射半边表*/

//输出函数
void Show_Node_List(vector<Point_3>& Point_List);//输出模型节点链表
void Show_Edge_List(Model& MyModel);//输出模型边链表
void Show_Triangle_List(Model& MyModel);//输出模型三角面链表
void Show_Edge_List_OfNode(Model& MyModel);//输出节点的入射半边链表
void Show_Tri_List_OfNode(Model& MyModel);//输出各节点的邻接面表
void Show_Tri_List_OfEdge(Model& MyModel);//输出边的邻接三角形
void Show_Tri_List_OfFace(Model& MyModel);//输出三角形的邻接面片
void Show_Half_Edge_List(vector<HalfEdge> H_E_L);//输出模型的半边表
void Write_Data(Model& MyModel);//将无重复顶点和三角面片写入中间文件
void Write_Data_F(Model& MyModel);

#endif