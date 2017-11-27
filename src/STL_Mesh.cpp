#include <iostream>
#include <ctime>
#include "STL__Kernel.h"
#include "STL__Function.h"

int main()
{
	clock_t start,finish;
	double totaltime;
	double MySize = 0.0;
	start=clock();

	vector<Point_3> Point_List;
	vector<HalfEdge> Half_Edge_List;
	vector<Face_Tri> Tri_List;
	Model MyModel;

	Read_STL(Point_List,Half_Edge_List, Tri_List);//读取STL文件
	MyModel = Build_ToPo_Of_Model(Point_List, Half_Edge_List, Tri_List);//建立拓扑关系
	Get_Feature_Line(MyModel);//几何特征识别

	//vector<vector<int>> Div_Result;
	/*vector<Sub_Set> Div_Result = Surface_div(MyModel);;
	vector<vector<int>>::iterator it_tt;
	vector<int>::iterator it_t;
	for (int i = 0; i < Div_Result.size(); i++)
	{
		for (it_t = Div_Result[i].Sub_Face.begin(); it_t != Div_Result[i].Sub_Face.end(); it_t++)
			std::cout << *it_t << "  ";
		std::cout << std::endl;
		set<int>::iterator it_s;
		for (it_s = Div_Result[i].Loop_Edge.begin(); it_s != Div_Result[i].Loop_Edge.end(); it_s++)
		{
			cout << *it_s << " ";
		}
		std:cout << std::endl;
	}*/

	finish=clock();
	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	std::cout<<"\n此程序的运行时间为"<<totaltime<<"秒！"<< std::endl;

	//Show_Node_List(MyModel.point());
	//Show_Edge_List(MyModel);
	//Show_Edge_List_OfNode(MyModel);
	//Show_Tri_List_OfNode(MyModel);
	//Show_Triangle_List(MyModel);
	//Show_Tri_List_OfEdge(MyModel);
	//Show_Tri_List_OfFace(MyModel);
	Write_Data(MyModel);
	//Write_Data_F(MyModel);
	std::cin.get();
	return 0;
}