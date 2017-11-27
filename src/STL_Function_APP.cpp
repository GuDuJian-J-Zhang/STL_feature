#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
#include <stack>
#include <set>
#include <queue>
#include <iterator>
#include "STL__Function.h"
using namespace std;

ostream & operator<<(ostream &os,const Point_3 &P)
{
	//os << P.V_ID() << setw(15) << P.X() << setw(15) << P.Y() << setw(15) << P.Z();
	os << P.X() << setw(15) << P.Y() << setw(15) << P.Z();
	return os;
}
ostream & operator<<(ostream &os,const Face_Tri &Tri)
{
	//os << Tri.F_ID() << setw(15) << Tri.P_ID(0) << setw(15) << Tri.P_ID(1) << setw(15) << Tri.P_ID(2);
	//os << setw(15) << Tri.P_ID(0) << setw(15) << Tri.P_ID(1) << setw(15) << Tri.P_ID(2);
	os << "3" << setw(15) << Tri.P_ID(0) - 1 << setw(15) << Tri.P_ID(1) - 1 << setw(15) << Tri.P_ID(2) - 1;
	return os;
}
ostream & operator<<(ostream &os,const Edge &e)
{
	os << e.P0() << ", " << e.P1();
	return os;
}
bool operator<(Point_3 P1, Point_3 P2)
{
	if (P1.X() < P2.X()-Tolerance || abs(P1.X() - P2.X()) < Tolerance && P1.Y() < P2.Y()-Tolerance || abs(P1.X() - P2.X()) < Tolerance && abs(P1.Y() - P2.Y()) < Tolerance && P1.Z() < P2.Z()-Tolerance)
		return 1;
	else return 0;
}
bool operator==(Point_3 P1, Point_3 P2)
{
	if (abs(P1.X() - P2.X()) < Tolerance && abs(P1.Y() - P2.Y()) && abs(P1.Z() - P2.Z()) < Tolerance)
		return 1;
	else return 0;
}

//STL文件读取
int Read_STL(vector<Point_3> &Point_List, vector<HalfEdge> &Half_Edge_List, vector<Face_Tri> &Tri_List)
{
	//读取STL文件,并建立三角面链表Tri_List和模型节点链表Point_List
	set<Point_3> Point_Set;
	ifstream read;
	read.open("123.stl",ios::in);
	string line;
	getline(read, line);
	string end = "end" + line;
	int n = 0;//number of Triangle
	int m = 0;//number of Node
	int k = 0;//number of Half_Edge
	set<Point_3>::const_iterator point_itr;
	Point_List.reserve(10000);
	Half_Edge_List.reserve(10000);
	Tri_List.reserve(10000);
	std::cout << "文件读取中...\n";
	while (1)
	{
		getline(read, line);
		if (line == end) break;
		getline(read, line);
		Face_Tri Tri;
		Point_3 p[3],q[3];//存储三角形三个点
		n++;
		Tri.F_ID() = n; // 设置面片编号
		if (n > 1)
		{
			for (int j = 0; j < 3; j++)
				q[j] = p[j];
		}
		for (int i = 0; i < 3; i++)
		{
			read >> line;
			Point_3 pt;
			read >> pt.X() >> pt.Y() >> pt.Z();
			int flag = 0;
			for (int j = 0; j < 3; j++)
			{
				if (pt == q[j])
			    {
					pt.V_ID() = q[j].V_ID();
					flag = 1;
					break;
				}
			}
			if (!flag)
			{
				point_itr = Point_Set.find(pt);
				if (point_itr == Point_Set.end())
				{
					m++;
					pt.V_ID() = m;
					Point_List.push_back(pt);//插入新点
					Point_Set.insert(pt);
				}
				else pt.V_ID() = point_itr->V_ID();
			}
			p[i] = pt;
			Tri.P_ID(i) = pt.V_ID();
			read.get();
		}
		//---------------------------创建半边表--------------------------------------//
		HalfEdge he0(p[0].V_ID(), n, ++k); Half_Edge_List.push_back(he0); 
		Point_List[p[2].V_ID() - 1].H_E().push_back(k);//设置节点的出射半边表
		Tri.E(0) = k;//设置三角面片的入射半边表
		HalfEdge he1(p[1].V_ID(), n, ++k); Half_Edge_List.push_back(he1);
		Point_List[p[0].V_ID() - 1].H_E().push_back(k);
		Tri.E(1) = k;
		HalfEdge he2(p[2].V_ID(), n, ++k); Half_Edge_List.push_back(he2);
		Point_List[p[1].V_ID() - 1].H_E().push_back(k);
		Tri.E(2) = k;

		Tri_List.push_back(Tri);
		getline(read, line);
		getline(read, line);
	}
	read.close();
	std::cout << "文件读取结束...\n";
	vector<Point_3>(Point_List).swap(Point_List);//压缩剩余空间
	vector<HalfEdge>(Half_Edge_List).swap(Half_Edge_List);
	vector<Face_Tri>(Tri_List).swap(Tri_List);
	cout << "Number of Node: " << m << endl;
	cout << "Number of Triangle: " << n << endl;
	return 1;
}

//拓扑重建
Model Build_ToPo_Of_Model(vector<Point_3> &Point_List, vector<HalfEdge> &H_E_L, vector<Face_Tri> &Tri_List)
{//用节点链表和三角形链表创建模型拓扑结构
    std::cout << "拓扑重建...\n";
	//-----------------------遍历三角面片创建模型边链表及各节点入射半边表--------------------------------------//
	vector<Face_Tri>::iterator it_Tri;
	for (it_Tri = Tri_List.begin(); it_Tri != Tri_List.end(); it_Tri++)
	{//对每个三角面片
		//-------------------------------计算三角面片的单位法矢------------------------------------------------//
		double D_X, D_Y, D_Z;
		Point_3 p0 = Point_List[(*it_Tri).P_ID(0) - 1];
		Point_3 p1 = Point_List[(*it_Tri).P_ID(1) - 1];
		Point_3 p2 = Point_List[(*it_Tri).P_ID(2) - 1];
		/*std::cout << p0 << std::endl;
		std::cout << p1 << std::endl;
		std::cout << p2 << std::endl;*/
		D_X = (p2.Y()- p1.Y()) * (p0.Z() - p1.Z()) - (p2.Z() - p1.Z()) * (p0.Y() - p1.Y());
		D_Y = -(p2.X() - p1.X()) * (p0.Z() - p1.Z()) + (p2.Z() - p1.Z()) * (p0.X() - p1.X());
		D_Z = (p2.X() - p1.X()) * (p0.Y() - p1.Y()) - (p2.Y() - p1.Y()) * (p0.X() - p1.X());

		double Length = sqrt(D_X * D_X + D_Y * D_Y + D_Z * D_Z);
		D_X /= Length; D_Y /= Length; D_Z /= Length;
		
		it_Tri->Dir().push_back(D_X); it_Tri->Dir().push_back(D_Y); it_Tri->Dir().push_back(D_Z);
		//-------------------------------创建模型的边表------------------------------------------------//
		/*Edge e0(it_Tri->P[0], it_Tri->P[1]), e1(it_Tri->P[1], it_Tri->P[2]), e2(it_Tri->P[2], it_Tri->P[0]);
		e0.num = 0; e1.num = 1; e2.num = 2; 
		e0.flag = 0; e1.flag = 0; e2.flag = 0;
		e0.Tri_List.push_back(it_Tri->num); e1.Tri_List.push_back(it_Tri->num); e2.Tri_List.push_back(it_Tri->num);
		//-------------------------------设置点的入射半边------------------------------------------------//
		Point_List[it_Tri->P[0].num - 1].Edge_List.push_back(e2);
		Point_List[it_Tri->P[1].num - 1].Edge_List.push_back(e0);
		Point_List[it_Tri->P[2].num - 1].Edge_List.push_back(e1);*/
		//std::cout << it_Tri->Dir()[0] << "  " << it_Tri->Dir()[1] << "  " << it_Tri->Dir()[2] << std::endl; 
	}
	//--------------------------------------------------------------------------------------------------------//
	vector<Edge> Edge_List;
	Build_Edge_List(Point_List, H_E_L, Edge_List);//创建模型的边链表,同时建立每条边的邻接三角形链表
	/*Build_TriandEdge_List_Of_Face(Tri_List, Edge_List);//为每个三角形面片创建边链表和邻接面表*/
	Model MyModel(Point_List, H_E_L, Edge_List, Tri_List);
	return MyModel;
}

void Build_Edge_List(vector<Point_3> &Point_List, vector<HalfEdge> &H_E_L, vector<Edge> &Edge_List)
{
	vector<Point_3>::iterator itp = Point_List.begin();
	int edge_num = 0;//边编号
	//Edge_List.reserve(20000);
	while (itp != Point_List.end())
	{
		vector<int>::iterator ite_of_p = itp->H_E().begin();
		while (ite_of_p != itp->H_E().end())
		{
			if (H_E_L[*ite_of_p - 1].E_ID() != 0)
				ite_of_p++;
			else
			{//创建一条新边
				edge_num++;//边编号
				Edge e(H_E_L[*ite_of_p- 1].Half_P(), itp->V_ID(), edge_num);
				H_E_L[*ite_of_p- 1].E_ID() = edge_num;//设置半边所属的边的编号
				e.H_Edge().push_back(*ite_of_p);//设置组成边的半边(每条边最多由两条半边组成)
				//查找另一条半边
				Point_3 p_temp = Point_List[H_E_L[*ite_of_p- 1].Half_P() - 1];
				//在p_temp的出射半边表中查找入射点为itp的半边
				vector<int>::iterator ite_of_p_temp = p_temp.H_E().begin();
				while (ite_of_p_temp != p_temp.H_E().end())
				{
					if (H_E_L[*ite_of_p_temp - 1].Half_P() == itp->V_ID() && H_E_L[*ite_of_p_temp - 1].E_ID() == 0)
					{
						H_E_L[*ite_of_p_temp - 1].E_ID() = edge_num;
						e.H_Edge().push_back(*ite_of_p_temp);//设置组成边的另一条半边
						break;
					}
					else ite_of_p_temp++;
				}
				Point_3 p_temp0 = Point_List[H_E_L[*ite_of_p- 1].Half_P() - 1];
				Point_3 p_temp1 = Point_List[itp->V_ID() - 1];
				double length = (p_temp0.X() - p_temp1.X()) * (p_temp0.X() - p_temp1.X())
					           +(p_temp0.Y() - p_temp1.Y()) * (p_temp0.Y() - p_temp1.Y())
							   +(p_temp0.Z() - p_temp1.Z()) * (p_temp0.Z() - p_temp1.Z());
				e.Length() = length;//计算边的长度*/
				Edge_List.push_back(e);
				//cout << ite_of_p->num << endl;
			}
		}
		itp++;
	}
	std::cout << "拓扑重建完成...\n";
	cout << "Number of Edge: " << edge_num <<endl;
}

//网格重分
void Get_Feature_Line(Model& MyModel)
{//提取模型的特征边和特征点
	//----------------------------------------------基于二面角--------------------------------------------------------------------//
	double alpha1 = 0.0;//小于10度的最大值
	double alpha2 = 180.0;//大于10度的最小值
	double alpha_fa = 0.0;//二面角阀值
	struct Ratio_Of_Edge
	{
		int num;//边的编号
		double r_o_e;//周长比值
	};//存储边的周长比值
	vector<Ratio_Of_Edge> Ratio;
	//-----------------------------计算每条边的二面角值和周长比值----------------------------------------------------//
	vector<Edge>::iterator it_Edge;
	for (it_Edge = MyModel.edge().begin(); it_Edge != MyModel.edge().end(); it_Edge++)
	{//对模型的每条边
		if (it_Edge->H_Edge().size() != 2)
		{
			it_Edge->Angle() = 180.0;
			std::cout << *it_Edge << "  " << it_Edge->H_Edge().size() << std::endl;
		}
		else
		{
			//cout << "#" << it_Edge->P0() << ", " << it_Edge->P1() << endl;
			vector<int> F_L_E = Find_Tri_List_OfEdge(it_Edge, MyModel.h_edge());
			//-----------------------计算邻接三角形的周长比-------------------------------------------//
			int he0, he1, he2;
			he0 = MyModel.face()[F_L_E[0] - 1].E(0);
			he1 = MyModel.face()[F_L_E[0] - 1].E(1);
			he2 = MyModel.face()[F_L_E[0] - 1].E(2);
			int e0, e1, e2;
			e0 = MyModel.h_edge()[he0 - 1].E_ID();
			e1 = MyModel.h_edge()[he1 - 1].E_ID();
			e2 = MyModel.h_edge()[he2 - 1].E_ID();
			//std::cout << "!!!  " << e0 << "  " << e1 << " " << e2 << std::endl;
			double lengthA;
			lengthA = MyModel.edge()[e0 - 1].Length() + MyModel.edge()[e1 - 1].Length() + MyModel.edge()[e2 - 1].Length();
			he0 = MyModel.face()[F_L_E[1] - 1].E(0);
			he1 = MyModel.face()[F_L_E[1] - 1].E(1);
			he2 = MyModel.face()[F_L_E[1] - 1].E(2);

			e0 = MyModel.h_edge()[he0 - 1].E_ID();
			e1 = MyModel.h_edge()[he1 - 1].E_ID();
			e2 = MyModel.h_edge()[he2 - 1].E_ID();
			double lengthB = MyModel.edge()[e0 - 1].Length() + MyModel.edge()[e1 - 1].Length() + MyModel.edge()[e2 - 1].Length();
			if (lengthA < lengthB)
			{//保证lengthA > lengthB
				double temp;
				temp = lengthA;
				lengthA = lengthB;
				lengthB = temp;
			}
			Ratio_Of_Edge rae;
			rae.num = it_Edge->E_ID();  rae.r_o_e = sqrt(lengthA / lengthB);
			Ratio.push_back(rae);
			//---------------------------------计算二面角---------------------------------------------//
			double A = sqrt(MyModel.face()[F_L_E[0] - 1].Dir()[0] * MyModel.face()[F_L_E[0] - 1].Dir()[0] 
			               +MyModel.face()[F_L_E[0] - 1].Dir()[1] * MyModel.face()[F_L_E[0] - 1].Dir()[1]
			               +MyModel.face()[F_L_E[0] - 1].Dir()[2] * MyModel.face()[F_L_E[0] - 1].Dir()[2]
			               );
			double B = sqrt(MyModel.face()[F_L_E[1] - 1].Dir()[0] * MyModel.face()[F_L_E[1] - 1].Dir()[0] 
			               +MyModel.face()[F_L_E[1] - 1].Dir()[1] * MyModel.face()[F_L_E[1] - 1].Dir()[1]
			               +MyModel.face()[F_L_E[1] - 1].Dir()[2] * MyModel.face()[F_L_E[1] - 1].Dir()[2]
			               );
			double AB = MyModel.face()[F_L_E[0] - 1].Dir()[0] * MyModel.face()[F_L_E[1] - 1].Dir()[0]
			           +MyModel.face()[F_L_E[0] - 1].Dir()[1] * MyModel.face()[F_L_E[1] - 1].Dir()[1]
			           +MyModel.face()[F_L_E[0] - 1].Dir()[2] * MyModel.face()[F_L_E[1] - 1].Dir()[2];
			it_Edge->Angle() = acos(AB / (A * B)) * 180 / PI;
			if (it_Edge->Angle() > alpha1 && it_Edge->Angle() < 10.0) alpha1 = it_Edge->Angle();
			else if (it_Edge->Angle() < alpha2 && it_Edge->Angle() > 10.0) alpha2 = it_Edge->Angle();
			//cout << it_Edge->Angle() << endl;
		}
	}
	alpha_fa = (alpha1 + alpha2) / 2.0;
	cout << "alpha = " << alpha_fa << endl;
	for (it_Edge = MyModel.edge().begin(); it_Edge != MyModel.edge().end(); it_Edge++)
	{
		if (it_Edge->Angle() > alpha_fa) 
		{
			it_Edge->Feature() = 1;//设置该边为特征边
			MyModel.point()[it_Edge->P0() - 1].Feature() = 1;//设置边的两端点为特征点
			MyModel.point()[it_Edge->P1() - 1].Feature() = 1;
		}
	}
	//--------------------------------------------------------------------基于边长比值--------------------------------------------------------------//
	vector<Ratio_Of_Edge>::iterator it_roe;
	double sum = 0.0;
	for (it_roe = Ratio.begin(); it_roe != Ratio.end(); it_roe++)
	{
		sum += it_roe->r_o_e;
	}
	double Tao = 4 * sum / Ratio.size();
	for (it_roe = Ratio.begin(); it_roe != Ratio.end(); it_roe++)
	{
		if (it_roe->r_o_e >= Tao)
		{
			int ne = it_roe->num;
			if (MyModel.edge()[ne - 1].Feature() != 1)
			{
				MyModel.edge()[ne - 1].Feature() = 2;//设置该边为特征边
				MyModel.point()[MyModel.edge()[ne - 1].P0() - 1].Feature() = 1;//设置边的两端点为特征点
				MyModel.point()[MyModel.edge()[ne - 1].P1() - 1].Feature() = 1;
			}
		}
	}
}

vector<int> Find_Tri_List_OfEdge(vector<Edge>::iterator e, vector<HalfEdge>& H_L)
{
	vector<int> F_L;
	vector<int>::iterator it_H_E_L;
	for (it_H_E_L = e->H_Edge().begin(); it_H_E_L != e->H_Edge().end(); it_H_E_L++)
	{
		int fa = H_L[*it_H_E_L - 1].Half_F();
		F_L.push_back(fa);
	}
	return F_L;
}

vector<Sub_Set> Surface_div(Model& MyModel)
{//返回子域及子域边界(边界未整合)
	std::cout << "子域分解...\n";
	vector<Sub_Set> Div_Result;//子域分解结果
	vector<Face_Tri>& F_L = MyModel.face();
	vector<Edge>& E_L = MyModel.edge();
	vector<HalfEdge>& H_E_L = MyModel.h_edge();
	int temp = 0;
	int Num_Set = 0;//统计子域个数
	while (temp = Search_Face(F_L))
	{//如果有尚未分组的面片 
		Sub_Set Sub_Surface;
		Num_Set++;
		stack<int> fl_Temp;//临时栈 用于存储待处理的面片
		fl_Temp.push(temp);

		//以下通过广度优先搜索 对面片进行分组
		while (fl_Temp.size() != 0)
		{//如果栈不为空
			int f = fl_Temp.top();//弹出栈顶元素
			fl_Temp.pop();
			//cout << f << "...\n";
			Sub_Surface.Sub_Face.push_back(f);//将面片编号加入面片组
			F_L[f - 1].sub_face() = Num_Set;//设置面片所属的面片组序号
			for (int j = 0; j < 3; j++)
			{//对该面片的每条边
				int hei = F_L[f - 1].E(j);
				int ei = H_E_L[hei - 1].E_ID();
				set<int>::iterator it_s;
				it_s = Sub_Surface.Loop_Edge.find(ei);
				if (it_s != Sub_Surface.Loop_Edge.end())
					Sub_Surface.Loop_Edge.erase(ei);
				else Sub_Surface.Loop_Edge.insert(ei);
				if (E_L[ei - 1].Feature() == 0)
				{//若果该边不是特征边
					for (vector<int>::iterator it_H_E = E_L[ei - 1].H_Edge().begin(); it_H_E != E_L[ei - 1].H_Edge().end(); it_H_E++)
					{//对该边的每条半边
						if (*it_H_E != hei)
						{
							if (F_L[H_E_L[*it_H_E - 1].Half_F() - 1].sub_face() == 0)//且邻接面不在面片组中
							    fl_Temp.push(H_E_L[*it_H_E - 1].Half_F());//将该面片压入栈
							break;
						}
					}
				}
			}
		}
		//cout << "!!!\n";
		Div_Result.push_back(Sub_Surface);
	}
	std::cout << "子域分解结束...\n";
	std::cout << "Number of SubSet: " << Div_Result.size() << std::endl;
	return Div_Result;
}

int Search_Face(vector<Face_Tri> F_L)
{//查找有无未分组的面片，如果有 返回其编号  否则返回0
	for (int i = 0; i < F_L.size(); i++)
	{
		if (F_L[i].sub_face() == 0)
		{
			return F_L[i].F_ID();
			break;
		}
	}
	return 0;
}

//输出函数
void Show_Node_List(vector<Point_3>& Point_List)
{//输出模型的节点链表
	vector<Point_3>::iterator itp;
	cout << "模型节点列表:\n";
	for (itp = Point_List.begin(); itp != Point_List.end(); itp++)
	{
		if (itp->Feature() == 1)
		  cout << *itp << endl;
	}
}

void Show_Edge_List(Model& MyModel)
{//输出模型的边链表
	vector<Edge>& Edge_List = MyModel.edge(); 
	vector<Edge>::iterator ite;
	cout << "模型边列表:\n";
	for (ite = Edge_List.begin(); ite != Edge_List.end(); ite++)
	{
		//if (ite->Feature() == 2)
		//{
			cout << "# " << ite->E_ID() << endl;
			cout << "(" << ite->P0() << ", " << ite->P1() << ")" << "... " << ite->Feature() << endl;
		//}
	}
}

void Show_Triangle_List(Model& MyModel)
{
	vector<Face_Tri>& Tri_List = MyModel.face();
	cout << "模型的三角面片列表:\n";
	vector<Face_Tri>::iterator itt;
	for (itt = Tri_List.begin(); itt != Tri_List.end(); itt++)
	{
		cout << *itt << endl;
		/*cout << itt->Direction[0] << "    " << itt->Direction[1] << "     " << itt->Direction[2] << endl;*/
	}
}

void Show_Edge_List_OfNode(Model& MyModel)
{//输出节点的出射半边表
	vector<Point_3>& P_L = MyModel.point();
	vector<HalfEdge>& H_E_L = MyModel.h_edge();
	vector<Edge>& E_L = MyModel.edge();
	vector<Point_3>::iterator it_Point;
	cout << "各节点的出射半边列表:\n";
	for (it_Point = P_L.begin(); it_Point != P_L.end(); it_Point++)
	{//对每个节点
		//if (ite->flag < 0)
		cout << "#" << it_Point->V_ID() << endl;
		vector<int>::iterator it_Edge;
		for (it_Edge = it_Point->H_E().begin(); it_Edge != it_Point->H_E().end(); it_Edge++)
		{
			int h_e_temp = *it_Edge;//节点的出射半边编号
			int p_temp = H_E_L[h_e_temp - 1].Half_P();
			cout << "(" << it_Point->V_ID() << ", " << p_temp << ")" << endl;
		}
	}
}

void Show_Tri_List_OfNode(Model& MyModel)
{//输出各节点的邻接面表
	//索引方式:点-->半边-->面
	vector<Point_3> P_L = MyModel.point();
	vector<HalfEdge> H_E_L = MyModel.h_edge();
	vector<int>::iterator it_H_E;//半边索引
	vector<Point_3>::iterator it_P;
	cout << "各节点的邻接面表:\n";
	for (it_P = P_L.begin(); it_P != P_L.end(); it_P++)
	{//对每个顶点
		cout << "#" << it_P->V_ID() << endl;
		for (it_H_E = it_P->H_E().begin(); it_H_E != it_P->H_E().end(); it_H_E++)
		{
			cout << H_E_L[*it_H_E - 1].Half_F() << ", ";
		}
		cout << endl;
	}
	//vector<Point_3>
}

void Show_Tri_List_OfEdge(Model& MyModel)
{//输出模型每条边的邻接面表
	//索引方式:边-->半边-->面
	vector<Edge>& E_L = MyModel.edge();
	vector<HalfEdge>& H_L = MyModel.h_edge();
	vector<Edge>::iterator it_E;
	for (it_E = E_L.begin(); it_E != E_L.end(); it_E++)
	{//对每条边
		cout << "#" << it_E->E_ID() << "(" << it_E->P0() << ", " << it_E->P1() << ")" << endl;
		vector<int>::iterator it_H_E_L;//半边索引
		for (it_H_E_L = it_E->H_Edge().begin(); it_H_E_L != it_E->H_Edge().end(); it_H_E_L++)
		{//对该边所包含的每条半边
			cout << H_L[*it_H_E_L - 1].Half_F() << ", " << endl;
		}
	}
}

void Show_Tri_List_OfFace(Model& MyModel)
{//输出三角形的邻接面片
	//索引方式:面-->半边-->边-->半边-->面
	vector<Face_Tri>& F_T = MyModel.face();
	vector<HalfEdge>& H_E_L = MyModel.h_edge();
	vector<Edge>& E_L = MyModel.edge();
	vector<Face_Tri>::iterator it_F;
	for (it_F = F_T.begin(); it_F != F_T.end(); it_F++)
	{//对每个三角面片
		cout << "#" << it_F->F_ID() << endl;
		cout << "(" << it_F->P_ID(0) << ", " << it_F->P_ID(1) << ", " << it_F->P_ID(2) << ")" << endl;
		vector<int>::iterator it_H_E;
		for (int i = 0; i < 3; i++)
		{//对每条入射半边
			int hei = it_F->E(i);//半边编号
			int ei = H_E_L[hei - 1].E_ID();//半边所属的边的编号
			for (it_H_E = E_L[ei - 1].H_Edge().begin(); it_H_E != E_L[ei - 1].H_Edge().end(); it_H_E++)
			{//对该边的每条半边
				if (*it_H_E != hei)
				{
					cout << H_E_L[*it_H_E - 1].Half_F() << ", ";
					break;
				}
			}
		}
		cout << endl;
	}
}

void Show_Half_Edge_List(vector<HalfEdge> H_E_L)
{//输出模型的半边表
	vector<HalfEdge>::iterator it_h_e_l;
	for (it_h_e_l = H_E_L.begin(); it_h_e_l != H_E_L.end(); it_h_e_l++)
	{
		cout << "#" << it_h_e_l->HE_ID() << "  ";
		cout << it_h_e_l->Half_P() << endl;
	}
}

void Write_Data(Model& MyModel)
{
	vector<Point_3>& P_L = MyModel.point();
	vector<Face_Tri>& F_L = MyModel.face();
	ofstream outfile;
	outfile.open("123.face",ios::out);
	vector<Point_3>::iterator itp;
	//outfile << "Coordinates\n";
	outfile << "OFF\n";
	outfile << P_L.size() << "  " << P_L.size() << "  " << "0\n\n";
	for (itp = P_L.begin(); itp != P_L.end(); itp++)
	{
		outfile << *itp << endl;
	}
	//outfile << "end coordinates\n\n";
	vector<Face_Tri>::iterator itt;
	//outfile << "Elements\n";
	for (itt = F_L.begin(); itt != F_L.end(); itt++)
		outfile << *itt << endl;
	outfile << endl;
	//outfile << "end elements\n";
	outfile.close();
	std::cout << "数据写出完成...\n";
}

void Write_Data_F(Model& MyModel)
{
	vector<Point_3>& P_L = MyModel.point();
	vector<Edge>& E_L = MyModel.edge();
	ofstream outfile;
	outfile.open("123.edge",ios::out);
	vector<Edge>::iterator ite;
	for (ite = E_L.begin(); ite != E_L.end(); ite++)
	{
		if (ite->Feature() != 0)
		{//输出特征线的两个端点
			int NA = ite->P0();
			int NB = ite->P1();
			outfile << P_L[NA - 1] << endl;
			outfile << P_L[NB - 1] << endl;
		}
	}
	outfile.close();
}